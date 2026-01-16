#!/usr/bin/env python
"""FEP window runner for Vast.ai GPU deployment."""
import os
import sys
import numpy as np
import time
from pathlib import Path
from openmm import *
from openmm.app import *
from openmm import unit

# Configuration
TIMESTEP = 2.0 * unit.femtoseconds
TEMPERATURE = 310.0 * unit.kelvin
FRICTION = 1.0 / unit.picoseconds
EQUIL_STEPS = 50000      # 100 ps equilibration
PROD_STEPS = 500000      # 1 ns production
ENERGY_INTERVAL = 500    # Save every 1 ps


def validate_inputs(base_dir, sys_name, window_idx):
    """Validate inputs before running. Returns (sys_dir, lambda_schedule, n_windows) or raises."""
    sys_dir = Path(base_dir) / sys_name

    # Check required files exist
    required_files = [
        sys_dir / "alchemical_system.xml",
        sys_dir / "topology.pdb",
        sys_dir / "positions.npy",
        sys_dir / "lambda_schedule.npy"
    ]
    for f in required_files:
        if not f.exists():
            raise FileNotFoundError(f"Required file missing: {f}")

    # Load and validate lambda schedule
    lambda_schedule = np.load(sys_dir / "lambda_schedule.npy")
    if lambda_schedule.ndim != 2 or lambda_schedule.shape[1] < 2:
        raise ValueError(f"Lambda schedule has wrong shape: {lambda_schedule.shape}, expected (N, 2+)")

    n_windows = len(lambda_schedule)

    # Validate window index (0-based)
    if window_idx < 0 or window_idx >= n_windows:
        raise ValueError(f"Window index {window_idx} out of range [0, {n_windows-1}]")

    return sys_dir, lambda_schedule, n_windows


def run_window(base_dir, sys_name, window_idx, gpu_id):
    """Run a single FEP window."""

    # Validate inputs first
    sys_dir, lambda_schedule, n_windows = validate_inputs(base_dir, sys_name, window_idx)

    # Now create output directory
    window_dir = sys_dir / f"window_{window_idx:02d}"
    window_dir.mkdir(parents=True, exist_ok=True)

    lambda_elec = float(lambda_schedule[window_idx, 0])
    lambda_sterics = float(lambda_schedule[window_idx, 1])

    print(f"\n{'='*60}")
    print(f"Window {window_idx}: {sys_name}")
    print(f"lambda_elec={lambda_elec:.4f}, lambda_sterics={lambda_sterics:.4f}")
    print(f"GPU: {gpu_id}")
    print(f"{'='*60}")

    # Load system
    print("Loading system...")
    with open(sys_dir / "alchemical_system.xml", 'r') as f:
        system = XmlSerializer.deserialize(f.read())

    # Load topology and positions
    pdb = PDBFile(str(sys_dir / "topology.pdb"))
    positions_nm = np.load(sys_dir / "positions.npy")
    positions = [Vec3(x, y, z) * unit.nanometers for x, y, z in positions_nm]

    # Set lambda values in system defaults
    found_elec = False
    found_sterics = False
    for i in range(system.getNumForces()):
        force = system.getForce(i)
        if hasattr(force, 'getNumGlobalParameters'):
            for j in range(force.getNumGlobalParameters()):
                name = force.getGlobalParameterName(j)
                if name == 'lambda_electrostatics':
                    force.setGlobalParameterDefaultValue(j, lambda_elec)
                    found_elec = True
                elif name == 'lambda_sterics':
                    force.setGlobalParameterDefaultValue(j, lambda_sterics)
                    found_sterics = True

    if not found_elec:
        raise RuntimeError("System missing 'lambda_electrostatics' parameter")
    if not found_sterics:
        raise RuntimeError("System missing 'lambda_sterics' parameter")

    # Create integrator
    integrator = LangevinMiddleIntegrator(TEMPERATURE, FRICTION, TIMESTEP)

    # Setup platform - try CUDA first, fallback to OpenCL
    simulation = None
    for platform_name in ['CUDA', 'OpenCL']:
        try:
            platform = Platform.getPlatformByName(platform_name)
            if platform_name == 'CUDA':
                properties = {'CudaDeviceIndex': str(gpu_id), 'CudaPrecision': 'mixed'}
            else:
                properties = {'OpenCLPlatformIndex': '0', 'OpenCLDeviceIndex': str(gpu_id), 'OpenCLPrecision': 'mixed'}

            simulation = Simulation(pdb.topology, system, integrator, platform, properties)
            print(f"  Using {platform_name} platform on GPU {gpu_id}")
            break
        except Exception as e:
            print(f"  {platform_name} failed: {e}")
            # Need new integrator for retry
            integrator = LangevinMiddleIntegrator(TEMPERATURE, FRICTION, TIMESTEP)

    if simulation is None:
        raise RuntimeError("No GPU platform available")

    simulation.context.setPositions(positions)

    # Set lambda in context (validated above that these exist)
    simulation.context.setParameter('lambda_electrostatics', lambda_elec)
    simulation.context.setParameter('lambda_sterics', lambda_sterics)

    # Minimize
    print("Minimizing...")
    simulation.minimizeEnergy(maxIterations=1000)

    state = simulation.context.getState(getEnergy=True)
    min_energy = state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
    print(f"  After minimize: {min_energy:.0f} kJ/mol")

    if np.isnan(min_energy) or abs(min_energy) > 1e10:
        raise RuntimeError(f"Bad energy after minimization: {min_energy}")

    # Set velocities
    simulation.context.setVelocitiesToTemperature(TEMPERATURE)

    # Equilibration
    print(f"Equilibrating ({EQUIL_STEPS} steps)...")
    t0 = time.time()
    simulation.step(EQUIL_STEPS)
    print(f"  Done in {time.time()-t0:.1f}s")

    # Production - collect u_nk data
    print(f"Production ({PROD_STEPS} steps)...")
    t0 = time.time()

    # Handle any remainder steps
    n_samples = PROD_STEPS // ENERGY_INTERVAL
    remainder = PROD_STEPS % ENERGY_INTERVAL
    if remainder > 0:
        print(f"  Note: {remainder} remainder steps will be run at end")

    u_nk = np.zeros((n_samples, n_windows))
    kT = unit.MOLAR_GAS_CONSTANT_R * TEMPERATURE

    for sample_idx in range(n_samples):
        simulation.step(ENERGY_INTERVAL)

        # Evaluate energy at all lambda states
        for k in range(n_windows):
            lam_e_k = float(lambda_schedule[k, 0])
            lam_s_k = float(lambda_schedule[k, 1])

            simulation.context.setParameter('lambda_electrostatics', lam_e_k)
            simulation.context.setParameter('lambda_sterics', lam_s_k)

            state = simulation.context.getState(getEnergy=True)
            energy_kT = state.getPotentialEnergy() / kT

            if np.isnan(energy_kT):
                raise RuntimeError(f"NaN energy at sample {sample_idx}, lambda state {k}")

            u_nk[sample_idx, k] = energy_kT

        # Reset to this window's lambda
        simulation.context.setParameter('lambda_electrostatics', lambda_elec)
        simulation.context.setParameter('lambda_sterics', lambda_sterics)

        if (sample_idx + 1) % 100 == 0:
            elapsed = time.time() - t0
            rate = (sample_idx + 1) * ENERGY_INTERVAL / elapsed
            print(f"  Sample {sample_idx+1}/{n_samples} ({rate:.0f} steps/s)")

    # Run remainder steps
    if remainder > 0:
        simulation.step(remainder)

    prod_time = time.time() - t0
    print(f"  Production done in {prod_time:.1f}s ({PROD_STEPS/prod_time:.0f} steps/s)")

    # Save results
    np.save(window_dir / "u_nk.npy", u_nk)

    state = simulation.context.getState(getPositions=True)
    pos = state.getPositions(asNumpy=True).value_in_unit(unit.nanometer)
    np.save(window_dir / "final_positions.npy", pos)

    box = state.getPeriodicBoxVectors()
    box_arr = np.array([[v.x, v.y, v.z] for v in box])
    np.save(window_dir / "final_box_vectors.npy", box_arr)

    print(f"[PASS] Window {window_idx} complete!")
    print(f"  u_nk shape: {u_nk.shape}")

    return True


def main():
    if len(sys.argv) < 5:
        print("Usage: python run_fep_gpu.py <base_dir> <sys_name> <window_start> <window_end> [gpu_id]")
        print("  window_start/end are 0-based indices (inclusive)")
        sys.exit(1)

    base_dir = sys.argv[1]
    sys_name = sys.argv[2]
    window_start = int(sys.argv[3])
    window_end = int(sys.argv[4])
    gpu_id = int(sys.argv[5]) if len(sys.argv) > 5 else 0

    print("="*60)
    print(f"FEP GPU Runner")
    print(f"System: {sys_name}")
    print(f"Windows: {window_start}-{window_end} (0-based, inclusive)")
    print(f"GPU: {gpu_id}")
    print("="*60)

    failures = []
    for window_idx in range(window_start, window_end + 1):
        try:
            run_window(base_dir, sys_name, window_idx, gpu_id)
        except Exception as e:
            print(f"[FAIL] Window {window_idx}: {e}")
            import traceback
            traceback.print_exc()
            failures.append(window_idx)

    print("\n" + "="*60)
    if failures:
        print(f"[FAILED] {len(failures)} windows failed: {failures}")
        sys.exit(1)
    else:
        print("[SUCCESS] All windows complete!")
    print("="*60)


if __name__ == "__main__":
    main()
