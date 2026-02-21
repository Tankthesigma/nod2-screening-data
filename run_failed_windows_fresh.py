#!/usr/bin/env python
"""
RUN FAILED WINDOWS FRESH - Run the 13 failed windows from scratch (not seeded)

Based on canary test results: softcore_alpha=0.5 is FINE when starting fresh.
The original failures were due to checkpoint/seeding issues, not softcore_alpha.

Failed windows to run:
  wt_complex: 15, 16, 17, 18, 19
  mut_complex: 15, 17, 18, 19
  solvent: 16, 17, 18, 19

Strategy: Start each window fresh from minimized topology positions,
NOT seeded from previous windows.
"""
import os
import sys
import numpy as np
import argparse

from openmm import XmlSerializer, LangevinMiddleIntegrator, Platform, MonteCarloBarostat, Vec3
from openmm.app import PDBFile, Simulation, StateDataReporter
from openmm import unit

# Production parameters (conservative for stability at extreme lambda)
TEMPERATURE = 310.0 * unit.kelvin
FRICTION = 2.0 / unit.picosecond
TIMESTEP = 1.0 * unit.femtoseconds  # Conservative 1fs for stability
EQUIL_STEPS = 500000   # 500 ps equilibration
PROD_STEPS = 2000000   # 2 ns production
SAMPLE_INTERVAL = 1000  # Save every 1 ps
ENERGY_INTERVAL = 100000  # Report energy every 100 ps
CHECKPOINT_INTERVAL = 500000  # Checkpoint every 500 ps

# Failed windows
FAILED_WINDOWS = {
    "wt_complex": [15, 16, 17, 18, 19],
    "mut_complex": [15, 17, 18, 19],
    "solvent": [16, 17, 18, 19]
}


def get_box_vectors_from_pdb(pdb_path):
    """Extract box vectors from PDB CRYST1 record."""
    with open(pdb_path, 'r') as f:
        for line in f:
            if line.startswith('CRYST1'):
                a = float(line[6:15]) / 10.0  # Angstroms to nm
                b = float(line[15:24]) / 10.0
                c = float(line[24:33]) / 10.0
                return (
                    Vec3(a, 0, 0) * unit.nanometer,
                    Vec3(0, b, 0) * unit.nanometer,
                    Vec3(0, 0, c) * unit.nanometer
                )
    return None


def run_window_fresh(system_name, window_idx, base_path):
    """Run a single window fresh from minimized positions."""

    system_dir = os.path.join(base_path, system_name)
    window_dir = os.path.join(system_dir, f"window_{window_idx:02d}")

    print("=" * 60)
    print(f"RUNNING FRESH: {system_name}/window_{window_idx:02d}")
    print("=" * 60)

    # Load lambda schedule
    lambda_schedule = np.load(os.path.join(system_dir, "lambda_schedule.npy"))
    n_windows = len(lambda_schedule)
    lam_e, lam_s, lam_r = lambda_schedule[window_idx]
    print(f"Lambda: elec={lam_e:.4f}, sterics={lam_s:.4f}, restraints={lam_r:.4f}")

    has_restraints = system_name != "solvent"

    # Load system
    with open(os.path.join(system_dir, "alchemical_system.xml"), "r") as f:
        system = XmlSerializer.deserialize(f.read())

    # Remove barostat for NVT production (standard FEP practice)
    # Equilibration is done with barostat implicitly through minimization
    for i in range(system.getNumForces()):
        force = system.getForce(i)
        if isinstance(force, MonteCarloBarostat):
            system.removeForce(i)
            print("Running NVT (barostat removed)")
            break

    # Load topology
    pdb_path = os.path.join(system_dir, "topology.pdb")
    pdb = PDBFile(pdb_path)
    topology = pdb.topology

    # Load positions (from original topology, not seed)
    positions = np.load(os.path.join(system_dir, "positions.npy"))
    positions = positions * unit.nanometer

    # Get box vectors from PDB
    box_vectors = get_box_vectors_from_pdb(pdb_path)
    if box_vectors is None:
        print("[FAIL] No CRYST1 record in topology.pdb")
        return False

    # Create integrator
    integrator = LangevinMiddleIntegrator(TEMPERATURE, FRICTION, TIMESTEP)

    # Select platform
    try:
        platform = Platform.getPlatformByName('CUDA')
        properties = {'Precision': 'mixed'}
        print("Using CUDA platform")
    except Exception:
        platform = Platform.getPlatformByName('CPU')
        properties = {}
        print("Using CPU platform")

    # Create simulation
    simulation = Simulation(topology, system, integrator, platform, properties)
    context = simulation.context

    # Set box vectors FIRST, then positions
    context.setPeriodicBoxVectors(*box_vectors)
    context.setPositions(positions)

    # Validate lambda parameters exist
    params = context.getParameters()
    if 'lambda_electrostatics' not in params:
        print("[FAIL] lambda_electrostatics parameter missing")
        return False
    if 'lambda_sterics' not in params:
        print("[FAIL] lambda_sterics parameter missing")
        return False

    # Set lambda parameters
    context.setParameter('lambda_electrostatics', float(lam_e))
    context.setParameter('lambda_sterics', float(lam_s))
    if has_restraints:
        if 'lambda_restraints' not in params:
            print("[FAIL] lambda_restraints missing for complex system")
            return False
        context.setParameter('lambda_restraints', float(lam_r))

    # Initialize velocities
    context.setVelocitiesToTemperature(TEMPERATURE)

    # Minimize
    print("\nMinimizing...")
    try:
        simulation.minimizeEnergy(maxIterations=2000)
        state = context.getState(getEnergy=True)
        energy = state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
        print(f"Energy after minimization: {energy:.1f} kJ/mol")
        if not np.isfinite(energy):
            print("[FAIL] NaN after minimization")
            return False
    except Exception as e:
        print(f"[FAIL] Minimization failed: {e}")
        return False

    # Create output directory
    os.makedirs(window_dir, exist_ok=True)

    # Equilibration
    equil_ps = EQUIL_STEPS / 1000
    print(f"\nEquilibrating ({EQUIL_STEPS} steps = {equil_ps:.0f} ps)...")
    try:
        chunk = 50000
        for i in range(0, EQUIL_STEPS, chunk):
            simulation.step(min(chunk, EQUIL_STEPS - i))
            state = context.getState(getEnergy=True)
            energy = state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
            if not np.isfinite(energy):
                print(f"[FAIL] NaN at step {i + chunk}")
                return False
            if (i + chunk) % 100000 == 0:
                print(f"  Step {i + chunk}: E = {energy:.1f} kJ/mol")
        print(f"[OK] Equilibration complete, E = {energy:.1f} kJ/mol")
    except Exception as e:
        print(f"[FAIL] Equilibration failed: {e}")
        return False

    # Production
    prod_ps = PROD_STEPS / 1000
    n_samples = PROD_STEPS // SAMPLE_INTERVAL
    print(f"\nProduction ({PROD_STEPS} steps = {prod_ps:.0f} ps, {n_samples} samples)...")

    # Prepare u_nk storage
    kT = (unit.MOLAR_GAS_CONSTANT_R * TEMPERATURE).value_in_unit(unit.kilojoules_per_mole)
    u_nk_data = []
    times = []

    try:
        for sample in range(n_samples):
            simulation.step(SAMPLE_INTERVAL)

            # Collect u_nk
            u_row = []
            for k in range(n_windows):
                lam_e_k, lam_s_k, lam_r_k = lambda_schedule[k]
                context.setParameter('lambda_electrostatics', float(lam_e_k))
                context.setParameter('lambda_sterics', float(lam_s_k))
                if has_restraints:
                    context.setParameter('lambda_restraints', float(lam_r_k))

                state = context.getState(getEnergy=True)
                e_kj = state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
                u = e_kj / kT
                u_row.append(u)

            # Restore lambda
            context.setParameter('lambda_electrostatics', float(lam_e))
            context.setParameter('lambda_sterics', float(lam_s))
            if has_restraints:
                context.setParameter('lambda_restraints', float(lam_r))

            # Check for NaN
            if not all(np.isfinite(u_row)):
                print(f"[FAIL] NaN in u_nk at sample {sample + 1}")
                return False

            u_nk_data.append(u_row)
            times.append((sample + 1) * SAMPLE_INTERVAL * TIMESTEP.value_in_unit(unit.picosecond))

            # Progress report every 10%
            report_interval = max(1, n_samples // 10)
            if (sample + 1) % report_interval == 0:
                state = context.getState(getEnergy=True)
                energy = state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
                print(f"  Sample {sample + 1}/{n_samples} ({times[-1]:.1f} ps): E = {energy:.1f} kJ/mol")

            # Checkpoint periodically
            checkpoint_samples = max(1, CHECKPOINT_INTERVAL // SAMPLE_INTERVAL)
            if (sample + 1) % checkpoint_samples == 0:
                simulation.saveCheckpoint(os.path.join(window_dir, "checkpoint.chk"))

        print("[OK] Production complete")
    except Exception as e:
        print(f"[FAIL] Production failed: {e}")
        return False

    # Save outputs
    print("\nSaving outputs...")

    # Save u_nk
    u_nk = np.array(u_nk_data)
    np.save(os.path.join(window_dir, "u_nk.npy"), u_nk)
    print(f"  u_nk.npy: {u_nk.shape}")

    # Save times
    np.save(os.path.join(window_dir, "times.npy"), np.array(times))

    # Save final checkpoint
    simulation.saveCheckpoint(os.path.join(window_dir, "checkpoint.chk"))
    print("  checkpoint.chk")

    # Save final positions
    state = context.getState(getPositions=True)
    pos = state.getPositions(asNumpy=True).value_in_unit(unit.nanometer)
    np.save(os.path.join(window_dir, "final_positions.npy"), pos)
    print("  final_positions.npy")

    # Save final box vectors
    box = state.getPeriodicBoxVectors()
    box_array = np.array([
        [box[0][0].value_in_unit(unit.nanometer), box[0][1].value_in_unit(unit.nanometer), box[0][2].value_in_unit(unit.nanometer)],
        [box[1][0].value_in_unit(unit.nanometer), box[1][1].value_in_unit(unit.nanometer), box[1][2].value_in_unit(unit.nanometer)],
        [box[2][0].value_in_unit(unit.nanometer), box[2][1].value_in_unit(unit.nanometer), box[2][2].value_in_unit(unit.nanometer)]
    ], dtype=np.float64)
    np.save(os.path.join(window_dir, "final_box_vectors.npy"), box_array)
    print("  final_box_vectors.npy")

    print("\n" + "=" * 60)
    print(f"[PASS] {system_name}/window_{window_idx:02d} COMPLETE")
    print("=" * 60)
    return True


def main():
    parser = argparse.ArgumentParser(description="Run failed windows fresh")
    parser.add_argument("--system", help="System name (wt_complex, mut_complex, solvent)")
    parser.add_argument("--window", type=int, help="Window index")
    parser.add_argument("--base", default="C:/Users/vasud/nod2-screening-data/fep_pmx",
                        help="Base FEP directory")
    parser.add_argument("--all", action="store_true", help="Run all failed windows")

    args = parser.parse_args()

    if args.all:
        # Run all failed windows
        results = {}
        for system_name, windows in FAILED_WINDOWS.items():
            for window_idx in windows:
                success = run_window_fresh(system_name, window_idx, args.base)
                results[(system_name, window_idx)] = success

        # Summary
        print("\n" + "=" * 60)
        print("RESULTS SUMMARY")
        print("=" * 60)
        passed = sum(1 for s in results.values() if s)
        failed = len(results) - passed
        for (system_name, window_idx), success in results.items():
            status = "PASS" if success else "FAIL"
            print(f"  {system_name}/window_{window_idx:02d}: [{status}]")
        print(f"\nTotal: {passed} passed, {failed} failed")
        sys.exit(0 if failed == 0 else 1)

    elif args.system and args.window is not None:
        # Run single window
        success = run_window_fresh(args.system, args.window, args.base)
        sys.exit(0 if success else 1)

    else:
        print("Error: Specify --system and --window, or use --all")
        sys.exit(1)


if __name__ == "__main__":
    main()
