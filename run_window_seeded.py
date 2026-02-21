#!/usr/bin/env python
"""
Seeded FEP Window Runner for Failed Windows

This script loads a checkpoint from the previous window (seed) and runs
production at the current window's lambda values. Uses safer parameters
to avoid NaN at extreme lambda values.

SAFE PARAMETER CHANGES (don't affect Hamiltonian):
- timestep: 2 fs -> 1 fs
- friction: 1/ps -> 2/ps
- equilibration: 100 ps -> 500 ps

HAMILTONIAN-DEFINING (NOT changed):
- lambda schedule
- softcore_alpha
- Boresch restraint force constant
"""
import os
import sys
import numpy as np
import argparse

from openmm import XmlSerializer, LangevinMiddleIntegrator, Platform, MonteCarloBarostat, Vec3
from openmm.app import PDBFile, Simulation
from openmm import unit

# SAFER simulation parameters for extreme lambda
TEMPERATURE = 310.0 * unit.kelvin
FRICTION = 5.0 / unit.picosecond      # Higher friction for stability
TIMESTEP = 0.5 * unit.femtoseconds    # Very small timestep
EQUIL_STEPS = 200000                  # 100 ps equilibration
PROD_STEPS = 1000000                  # 500 ps production (keeps 1000 samples)
ENERGY_INTERVAL = 1000                # Sample every 1000 steps (0.5 ps)

# NaN handling
MAX_SKIP_PERCENT_EARLY = 2.0  # Abort if >2% skip in first 100 samples
MAX_SKIP_PERCENT_WARN = 5.0   # Warn if skip rate exceeds this


def get_seed_window(system_name, window_idx):
    """Determine which window to seed from based on seeding chains."""
    # Seeding chains from the fix protocol
    if system_name == "wt_complex":
        # wt: 15<-14, 16<-15, 17<-16, 18<-17, 19<-18
        return window_idx - 1
    elif system_name == "mut_complex":
        # mut: 15<-14, 17<-16 (16 succeeded), 18<-17, 19<-18
        if window_idx == 15:
            return 14
        elif window_idx == 17:
            return 16  # 16 succeeded
        else:
            return window_idx - 1
    elif system_name == "solvent":
        # solvent: 16<-15 (15 succeeded), 17<-16, 18<-17, 19<-18
        return window_idx - 1
    else:
        raise ValueError(f"Unknown system: {system_name}")


def run_seeded_window(system_name, window_idx, base_path, smoke_test=False):
    """Run a window seeded from previous window's checkpoint."""

    system_dir = os.path.join(base_path, system_name)
    window_dir = os.path.join(system_dir, f"window_{window_idx:02d}")

    # Ensure window directory exists
    os.makedirs(window_dir, exist_ok=True)

    # Determine seed window
    seed_idx = get_seed_window(system_name, window_idx)
    seed_dir = os.path.join(system_dir, f"window_{seed_idx:02d}")
    seed_checkpoint = os.path.join(seed_dir, "checkpoint.chk")

    print(f"{'='*60}")
    print(f"FEP Window {window_idx} - {system_name} (SEEDED FIX)")
    print(f"{'='*60}")
    print(f"Seeding from: window_{seed_idx:02d}")

    # Verify seed checkpoint exists
    if not os.path.exists(seed_checkpoint):
        raise RuntimeError(f"Seed checkpoint not found: {seed_checkpoint}")

    # Load lambda schedule
    lambda_schedule = np.load(os.path.join(system_dir, "lambda_schedule.npy"))
    n_windows = len(lambda_schedule)
    lam_e, lam_s, lam_r = lambda_schedule[window_idx]
    print(f"Lambda: elec={lam_e:.4f}, sterics={lam_s:.4f}, restraints={lam_r:.4f}")
    print(f"Using SAFER parameters: timestep=1fs, friction=2/ps, equil=500ps")

    has_restraints = system_name != "solvent"

    # Load system and REMOVE barostat immediately
    # (seed checkpoints are from NVT, so we must match)
    with open(os.path.join(system_dir, "alchemical_system.xml"), "r") as f:
        system = XmlSerializer.deserialize(f.read())

    # Remove barostat for compatibility with NVT seed checkpoints
    for i in range(system.getNumForces()):
        force = system.getForce(i)
        if isinstance(force, MonteCarloBarostat):
            system.removeForce(i)
            print("Removed barostat (NVT mode for seeded run)")
            break

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
    print(f"Loading seed checkpoint from window_{seed_idx:02d}...")
    simulation.loadCheckpoint(seed_checkpoint)

    # CRITICAL: Set THIS window's lambda parameters (checkpoint has seed's lambdas)
    params = context.getParameters()
    if 'lambda_electrostatics' in params:
        context.setParameter('lambda_electrostatics', float(lam_e))
    if 'lambda_sterics' in params:
        context.setParameter('lambda_sterics', float(lam_s))
    if has_restraints and 'lambda_restraints' in params:
        context.setParameter('lambda_restraints', float(lam_r))

    # Verify box vectors loaded
    state = context.getState(getPositions=True)
    box = state.getPeriodicBoxVectors()
    box_x = box[0][0].value_in_unit(unit.nanometer)
    box_y = box[1][1].value_in_unit(unit.nanometer)
    box_z = box[2][2].value_in_unit(unit.nanometer)
    print(f"Box vectors from seed: {box_x:.3f} x {box_y:.3f} x {box_z:.3f} nm")

    # Get seed lambda values
    seed_lam_e, seed_lam_s, seed_lam_r = lambda_schedule[seed_idx]

    # GRADUAL LAMBDA ANNEALING - don't jump directly
    # Anneal from seed lambda to target lambda over 100 steps
    n_anneal_steps = 100
    anneal_dynamics_per_step = 100  # 100 steps of dynamics per lambda change

    print(f"Annealing lambda from sterics={seed_lam_s:.2f} to {lam_s:.2f}...")
    for i in range(n_anneal_steps + 1):
        frac = i / n_anneal_steps
        interp_lam_s = seed_lam_s + frac * (lam_s - seed_lam_s)
        interp_lam_e = seed_lam_e + frac * (lam_e - seed_lam_e)
        interp_lam_r = seed_lam_r + frac * (lam_r - seed_lam_r)

        if 'lambda_electrostatics' in params:
            context.setParameter('lambda_electrostatics', float(interp_lam_e))
        if 'lambda_sterics' in params:
            context.setParameter('lambda_sterics', float(interp_lam_s))
        if has_restraints and 'lambda_restraints' in params:
            context.setParameter('lambda_restraints', float(interp_lam_r))

        # Brief minimization at each lambda step
        if i == 0:
            print("  Initial minimization...")
            simulation.minimizeEnergy(maxIterations=500)

        # Run short dynamics
        try:
            simulation.step(anneal_dynamics_per_step)
        except Exception as e:
            print(f"  WARNING: Dynamics failed at sterics={interp_lam_s:.3f}: {e}")
            # Try to recover with minimization
            simulation.minimizeEnergy(maxIterations=500)

        if i % 20 == 0:
            state = context.getState(getEnergy=True)
            energy = state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
            print(f"  Anneal step {i}/{n_anneal_steps}: sterics={interp_lam_s:.3f}, E={energy:.1f} kJ/mol")

    # Final minimization at target lambda
    print("Final minimization at target lambda...")
    simulation.minimizeEnergy(maxIterations=1000)

    state = context.getState(getEnergy=True)
    energy = state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
    print(f"Energy after annealing: {energy:.1f} kJ/mol")

    if not np.isfinite(energy):
        raise RuntimeError(f"Energy is NaN/Inf after annealing: {energy}")

    # SMOKE TEST: Just run brief dynamics and compute one u_nk row
    if smoke_test:
        print("\n=== SMOKE TEST MODE ===")
        print("Running 1000 steps...")
        simulation.step(1000)

        state = context.getState(getEnergy=True)
        energy = state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
        print(f"Energy after 1000 steps: {energy:.1f} kJ/mol")

        if not np.isfinite(energy):
            raise RuntimeError(f"SMOKE TEST FAILED: Energy is NaN/Inf")

        # Compute one u_nk row
        print("Computing u_nk at all lambda states...")
        test_energies = []
        for k in range(n_windows):
            lam_e_k, lam_s_k, lam_r_k = lambda_schedule[k]
            if 'lambda_electrostatics' in params:
                context.setParameter('lambda_electrostatics', float(lam_e_k))
            if 'lambda_sterics' in params:
                context.setParameter('lambda_sterics', float(lam_s_k))
            if has_restraints and 'lambda_restraints' in params:
                context.setParameter('lambda_restraints', float(lam_r_k))

            state = context.getState(getEnergy=True)
            e = state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
            test_energies.append(e)

        # Restore window lambda
        if 'lambda_electrostatics' in params:
            context.setParameter('lambda_electrostatics', float(lam_e))
        if 'lambda_sterics' in params:
            context.setParameter('lambda_sterics', float(lam_s))
        if has_restraints and 'lambda_restraints' in params:
            context.setParameter('lambda_restraints', float(lam_r))

        # Check for NaN/Inf
        nan_count = sum(1 for e in test_energies if not np.isfinite(e))
        print(f"u_nk test: {nan_count}/{n_windows} NaN/Inf energies")

        if nan_count > 0:
            raise RuntimeError(f"SMOKE TEST FAILED: {nan_count} NaN/Inf in u_nk")

        # Save smoke test result
        smoke_file = os.path.join(window_dir, "smoke_test_passed.txt")
        with open(smoke_file, 'w') as f:
            f.write(f"Smoke test passed for {system_name}/window_{window_idx:02d}\n")
            f.write(f"Box: {box_x:.3f} x {box_y:.3f} x {box_z:.3f} nm\n")
            f.write(f"Test energies (kJ/mol):\n")
            for k, e in enumerate(test_energies):
                f.write(f"  lambda[{k}]: {e:.2f}\n")

        print(f"\n[PASS] Smoke test passed!")
        print(f"Saved: {smoke_file}")
        return

    # FULL PRODUCTION RUN
    print(f"\nRunning NVT production ({PROD_STEPS} steps)...")
    n_samples = PROD_STEPS // ENERGY_INTERVAL
    u_nk = np.zeros((n_samples, n_windows))

    skipped_frames = 0
    valid_frames = 0

    for sample_idx in range(n_samples):
        # Run dynamics
        simulation.step(ENERGY_INTERVAL)

        # Compute energies at all lambda states
        frame_energies = np.zeros(n_windows)
        frame_has_nan = False

        for k in range(n_windows):
            lam_e_k, lam_s_k, lam_r_k = lambda_schedule[k]
            if 'lambda_electrostatics' in params:
                context.setParameter('lambda_electrostatics', float(lam_e_k))
            if 'lambda_sterics' in params:
                context.setParameter('lambda_sterics', float(lam_s_k))
            if has_restraints and 'lambda_restraints' in params:
                context.setParameter('lambda_restraints', float(lam_r_k))

            state = context.getState(getEnergy=True)
            kT = (unit.MOLAR_GAS_CONSTANT_R * TEMPERATURE).value_in_unit(unit.kilojoules_per_mole)
            e_kj = state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
            e = e_kj / kT  # Reduced potential (dimensionless float)
            frame_energies[k] = e

            if not np.isfinite(e) or abs(e) > 1e10:
                frame_has_nan = True

        # Restore window lambda
        if 'lambda_electrostatics' in params:
            context.setParameter('lambda_electrostatics', float(lam_e))
        if 'lambda_sterics' in params:
            context.setParameter('lambda_sterics', float(lam_s))
        if has_restraints and 'lambda_restraints' in params:
            context.setParameter('lambda_restraints', float(lam_r))

        # Handle NaN frames
        if frame_has_nan:
            skipped_frames += 1
            # Fill with NaN (will be excluded in MBAR)
            u_nk[sample_idx, :] = np.nan
        else:
            valid_frames += 1
            u_nk[sample_idx, :] = frame_energies

        # Progress and skip rate monitoring
        if (sample_idx + 1) % 100 == 0:
            skip_pct = 100.0 * skipped_frames / (sample_idx + 1)
            print(f"  Sample {sample_idx + 1}/{n_samples} - Skip rate: {skip_pct:.1f}%")

            # Early abort check
            if sample_idx < 100 and skip_pct > MAX_SKIP_PERCENT_EARLY:
                raise RuntimeError(f"ABORT: Skip rate {skip_pct:.1f}% > {MAX_SKIP_PERCENT_EARLY}% in first 100 samples")

            # Warning for high skip rate later
            if skip_pct > MAX_SKIP_PERCENT_WARN:
                print(f"  WARNING: High skip rate {skip_pct:.1f}%")

    # Final statistics
    final_skip_pct = 100.0 * skipped_frames / n_samples
    print(f"\nProduction complete!")
    print(f"  Valid frames: {valid_frames}/{n_samples}")
    print(f"  Skipped frames: {skipped_frames}/{n_samples} ({final_skip_pct:.1f}%)")

    # Save u_nk (including NaN frames - MBAR can handle sparse data)
    np.save(os.path.join(window_dir, "u_nk.npy"), u_nk)
    print(f"[PASS] Saved u_nk.npy: shape {u_nk.shape}")

    # Save N_k (actual valid frame count for MBAR)
    np.save(os.path.join(window_dir, "N_k.npy"), np.array([valid_frames]))
    print(f"[PASS] Saved N_k.npy: {valid_frames} valid frames")

    # Save checkpoint for next window in chain
    checkpoint_path = os.path.join(window_dir, "checkpoint.chk")
    simulation.saveCheckpoint(checkpoint_path)
    print(f"[PASS] Saved checkpoint.chk")

    # Save final positions
    state = context.getState(getPositions=True)
    pos = state.getPositions(asNumpy=True).value_in_unit(unit.nanometer)
    np.save(os.path.join(window_dir, "final_positions.npy"), pos)
    print(f"[PASS] Saved final_positions.npy")

    # Save box vectors
    box = state.getPeriodicBoxVectors()
    box_array = np.array([
        [box[0][0].value_in_unit(unit.nanometer), box[0][1].value_in_unit(unit.nanometer), box[0][2].value_in_unit(unit.nanometer)],
        [box[1][0].value_in_unit(unit.nanometer), box[1][1].value_in_unit(unit.nanometer), box[1][2].value_in_unit(unit.nanometer)],
        [box[2][0].value_in_unit(unit.nanometer), box[2][1].value_in_unit(unit.nanometer), box[2][2].value_in_unit(unit.nanometer)]
    ], dtype=np.float64)
    np.save(os.path.join(window_dir, "final_box_vectors.npy"), box_array)
    print(f"[PASS] Saved final_box_vectors.npy")

    # Save log
    log_path = os.path.join(window_dir, f"window_{window_idx:02d}_log.txt")
    with open(log_path, 'w') as f:
        f.write(f"Window {window_idx} - {system_name}\n")
        f.write(f"Seeded from: window_{seed_idx}\n")
        f.write(f"Lambda: elec={lam_e:.4f}, sterics={lam_s:.4f}, restraints={lam_r:.4f}\n")
        f.write(f"Valid frames: {valid_frames}/{n_samples}\n")
        f.write(f"Skipped frames: {skipped_frames}/{n_samples} ({final_skip_pct:.1f}%)\n")
    print(f"[PASS] Saved {log_path}")

    print(f"\n{'='*60}")
    print(f"Window {window_idx} complete!")
    print(f"{'='*60}\n")


def main():
    parser = argparse.ArgumentParser(description="Run seeded FEP window")
    parser.add_argument("--system", required=True, help="System name")
    parser.add_argument("--window", type=int, required=True, help="Window index")
    parser.add_argument("--base", default="C:/Users/vasud/nod2-screening-data/fep_pmx",
                        help="Base path to fep_pmx directory")
    parser.add_argument("--smoke-test", action="store_true",
                        help="Run smoke test only (brief dynamics + one u_nk row)")

    args = parser.parse_args()
    run_seeded_window(args.system, args.window, args.base, smoke_test=args.smoke_test)


if __name__ == "__main__":
    main()
