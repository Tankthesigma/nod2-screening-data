#!/usr/bin/env python
"""
RESCUE TEST - Short run to verify seeding works before full production

Tests if windows 15+ can run stably when seeded from previous checkpoint.
Uses conservative parameters: 1fs timestep, 2/ps friction, 500ps equil, 200ps prod
"""
import os
import sys
import numpy as np
import argparse

from openmm import XmlSerializer, LangevinMiddleIntegrator, Platform, MonteCarloBarostat, Vec3
from openmm.app import PDBFile, Simulation
from openmm import unit

# CONSERVATIVE parameters
TEMPERATURE = 310.0 * unit.kelvin
FRICTION = 2.0 / unit.picosecond
TIMESTEP = 1.0 * unit.femtoseconds
EQUIL_STEPS = 500000   # 500 ps equilibration
PROD_STEPS = 200000    # 200 ps production (short test)
ENERGY_INTERVAL = 1000 # Sample every 1 ps


def run_rescue_test(system_name, window_idx, seed_idx, base_path):
    """Run rescue test for a single window."""

    system_dir = os.path.join(base_path, system_name)
    window_dir = os.path.join(system_dir, f"window_{window_idx:02d}")
    seed_dir = os.path.join(system_dir, f"window_{seed_idx:02d}")
    seed_checkpoint = os.path.join(seed_dir, "checkpoint.chk")

    print("=" * 60)
    print(f"RESCUE TEST: {system_name}/window_{window_idx:02d}")
    print("=" * 60)
    print(f"Seed: window_{seed_idx:02d}")

    # Verify seed exists
    if not os.path.exists(seed_checkpoint):
        print(f"[FAIL] Seed checkpoint not found: {seed_checkpoint}")
        return False

    # Load lambda schedule
    lambda_schedule = np.load(os.path.join(system_dir, "lambda_schedule.npy"))
    n_windows = len(lambda_schedule)
    lam_e, lam_s, lam_r = lambda_schedule[window_idx]
    print(f"Target lambda: elec={lam_e:.4f}, sterics={lam_s:.4f}, restraints={lam_r:.4f}")

    has_restraints = system_name != "solvent"

    # Load system
    with open(os.path.join(system_dir, "alchemical_system.xml"), "r") as f:
        system = XmlSerializer.deserialize(f.read())

    # Check for barostat
    has_barostat = False
    for i in range(system.getNumForces()):
        if isinstance(system.getForce(i), MonteCarloBarostat):
            has_barostat = True
            break
    print(f"System has barostat: {has_barostat}")

    # Load topology
    pdb = PDBFile(os.path.join(system_dir, "topology.pdb"))
    topology = pdb.topology

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

    # Load seed checkpoint
    print(f"\nLoading checkpoint from window_{seed_idx:02d}...")
    try:
        simulation.loadCheckpoint(seed_checkpoint)
        print("[OK] Checkpoint loaded")
    except Exception as e:
        print(f"[FAIL] Checkpoint load failed: {e}")
        return False

    # Check state after loading
    state = context.getState(getPositions=True, getEnergy=True)
    energy = state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
    box = state.getPeriodicBoxVectors()
    box_x = box[0][0].value_in_unit(unit.nanometer)
    box_y = box[1][1].value_in_unit(unit.nanometer)
    box_z = box[2][2].value_in_unit(unit.nanometer)
    print(f"Energy after load: {energy:.1f} kJ/mol")
    print(f"Box: {box_x:.3f} x {box_y:.3f} x {box_z:.3f} nm")

    if not np.isfinite(energy):
        print("[FAIL] Energy is NaN/Inf after checkpoint load")
        return False

    # Set target lambda parameters
    print(f"\nSetting target lambda (sterics={lam_s:.2f})...")
    params = context.getParameters()
    if 'lambda_electrostatics' in params:
        context.setParameter('lambda_electrostatics', float(lam_e))
    if 'lambda_sterics' in params:
        context.setParameter('lambda_sterics', float(lam_s))
    if has_restraints and 'lambda_restraints' in params:
        context.setParameter('lambda_restraints', float(lam_r))

    # Check energy at new lambda
    state = context.getState(getEnergy=True)
    energy = state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
    print(f"Energy at target lambda: {energy:.1f} kJ/mol")

    if not np.isfinite(energy):
        print("[FAIL] Energy is NaN/Inf at target lambda")
        return False

    # Minimize
    print("\nMinimizing (1000 steps)...")
    try:
        simulation.minimizeEnergy(maxIterations=1000)
        state = context.getState(getEnergy=True)
        energy = state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
        print(f"Energy after minimization: {energy:.1f} kJ/mol")
        if not np.isfinite(energy):
            print("[FAIL] Energy is NaN/Inf after minimization")
            return False
    except Exception as e:
        print(f"[FAIL] Minimization failed: {e}")
        return False

    # Equilibration
    print(f"\nRunning equilibration ({EQUIL_STEPS} steps = 500 ps)...")
    try:
        # Run in chunks to detect failure early
        chunk_size = 50000
        for i in range(0, EQUIL_STEPS, chunk_size):
            simulation.step(min(chunk_size, EQUIL_STEPS - i))
            state = context.getState(getEnergy=True)
            energy = state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
            if not np.isfinite(energy):
                print(f"[FAIL] Energy became NaN at step {i + chunk_size}")
                return False
            if (i + chunk_size) % 100000 == 0:
                print(f"  Step {i + chunk_size}: E = {energy:.1f} kJ/mol")
        print(f"[OK] Equilibration complete, E = {energy:.1f} kJ/mol")
    except Exception as e:
        print(f"[FAIL] Equilibration failed: {e}")
        return False

    # Test u_nk evaluation at all lambda states
    print("\nTesting u_nk evaluation across all 20 lambda states...")
    test_energies = []
    nan_count = 0
    kT = (unit.MOLAR_GAS_CONSTANT_R * TEMPERATURE).value_in_unit(unit.kilojoules_per_mole)

    for k in range(n_windows):
        lam_e_k, lam_s_k, lam_r_k = lambda_schedule[k]
        if 'lambda_electrostatics' in params:
            context.setParameter('lambda_electrostatics', float(lam_e_k))
        if 'lambda_sterics' in params:
            context.setParameter('lambda_sterics', float(lam_s_k))
        if has_restraints and 'lambda_restraints' in params:
            context.setParameter('lambda_restraints', float(lam_r_k))

        state = context.getState(getEnergy=True)
        e_kj = state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
        u = e_kj / kT
        test_energies.append(u)

        if not np.isfinite(u):
            nan_count += 1
            print(f"  lambda[{k}] (sterics={lam_s_k:.2f}): NaN/Inf")

    # Restore target lambda
    if 'lambda_electrostatics' in params:
        context.setParameter('lambda_electrostatics', float(lam_e))
    if 'lambda_sterics' in params:
        context.setParameter('lambda_sterics', float(lam_s))
    if has_restraints and 'lambda_restraints' in params:
        context.setParameter('lambda_restraints', float(lam_r))

    if nan_count > 0:
        print(f"[FAIL] {nan_count}/{n_windows} lambda states have NaN/Inf energies")
        return False
    print(f"[OK] All 20 lambda states evaluated successfully")

    # Short production (200 ps)
    print(f"\nRunning short production ({PROD_STEPS} steps = 200 ps)...")
    n_samples = PROD_STEPS // ENERGY_INTERVAL
    try:
        for sample in range(n_samples):
            simulation.step(ENERGY_INTERVAL)
            if (sample + 1) % 50 == 0:
                state = context.getState(getEnergy=True)
                energy = state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
                if not np.isfinite(energy):
                    print(f"[FAIL] Energy became NaN at sample {sample + 1}")
                    return False
                print(f"  Sample {sample + 1}/{n_samples}: E = {energy:.1f} kJ/mol")
        print(f"[OK] Production complete")
    except Exception as e:
        print(f"[FAIL] Production failed: {e}")
        return False

    # Save outputs
    os.makedirs(window_dir, exist_ok=True)

    # Save checkpoint
    checkpoint_path = os.path.join(window_dir, "checkpoint.chk")
    simulation.saveCheckpoint(checkpoint_path)
    print(f"[OK] Saved checkpoint.chk")

    # Save final positions
    state = context.getState(getPositions=True)
    pos = state.getPositions(asNumpy=True).value_in_unit(unit.nanometer)
    np.save(os.path.join(window_dir, "final_positions.npy"), pos)
    print(f"[OK] Saved final_positions.npy")

    # Save box vectors
    box = state.getPeriodicBoxVectors()
    box_array = np.array([
        [box[0][0].value_in_unit(unit.nanometer), box[0][1].value_in_unit(unit.nanometer), box[0][2].value_in_unit(unit.nanometer)],
        [box[1][0].value_in_unit(unit.nanometer), box[1][1].value_in_unit(unit.nanometer), box[1][2].value_in_unit(unit.nanometer)],
        [box[2][0].value_in_unit(unit.nanometer), box[2][1].value_in_unit(unit.nanometer), box[2][2].value_in_unit(unit.nanometer)]
    ], dtype=np.float64)
    np.save(os.path.join(window_dir, "final_box_vectors.npy"), box_array)
    print(f"[OK] Saved final_box_vectors.npy")

    print("\n" + "=" * 60)
    print(f"[PASS] RESCUE TEST PASSED for {system_name}/window_{window_idx:02d}")
    print("=" * 60)
    return True


def main():
    parser = argparse.ArgumentParser(description="Run rescue test")
    parser.add_argument("--system", required=True)
    parser.add_argument("--window", type=int, required=True)
    parser.add_argument("--seed", type=int, required=True)
    parser.add_argument("--base", default="C:/Users/vasud/nod2-screening-data/fep_pmx")

    args = parser.parse_args()
    success = run_rescue_test(args.system, args.window, args.seed, args.base)
    sys.exit(0 if success else 1)


if __name__ == "__main__":
    main()
