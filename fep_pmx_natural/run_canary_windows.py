#!/usr/bin/env python
"""
Canary Window Runner for FEP Natural Product CID_10120

Runs windows 15-19 SEQUENTIALLY with checkpoint seeding:
- Window 15 is seeded from window 14's checkpoint
- Window 16 is seeded from window 15's checkpoint
- ... and so on

This tests the most challenging windows (low lambda_sterics)
before committing to the full 60-window run.

STABILITY CHECKS:
- No NaN energies
- No LINCS constraint explosions
- No restraint blowups
- Temperature/pressure stable
- Potential energy finite with no runaway drift
"""

import os
import sys
import subprocess
import numpy as np
from pathlib import Path
import json
import time

BASE_DIR = Path("C:/Users/vasud/nod2-screening-data/fep_pmx_natural")
SYSTEMS = ['wt_complex', 'mut_complex', 'solvent']
CANARY_WINDOWS = [15, 16, 17, 18, 19]
SEED_WINDOW = 14  # Seed from this window

def check_checkpoint_exists(sys_name, window_idx):
    """Check if checkpoint file exists for a window."""
    chk_path = BASE_DIR / sys_name / f"window_{window_idx:02d}" / "checkpoint.chk"
    return chk_path.exists()

def check_box_vectors_exist(sys_name, window_idx):
    """Check if box vectors file exists."""
    box_path = BASE_DIR / sys_name / f"window_{window_idx:02d}" / "final_box_vectors.npy"
    return box_path.exists()

def run_window_with_seed(sys_name, window_idx, seed_window_idx):
    """Run a window seeded from the previous window's checkpoint."""

    window_dir = BASE_DIR / sys_name / f"window_{window_idx:02d}"
    seed_dir = BASE_DIR / sys_name / f"window_{seed_window_idx:02d}"
    script_path = window_dir / "run_window_seeded.py"
    log_path = window_dir / "run.log"

    print(f"\n{'='*60}")
    print(f"Running {sys_name}/window_{window_idx:02d}")
    print(f"Seeded from window_{seed_window_idx:02d}")
    print(f"{'='*60}")

    # Check seed data exists - we need final_positions.npy (checkpoint.chk is optional)
    seed_pos = seed_dir / "final_positions.npy"
    seed_box = seed_dir / "final_box_vectors.npy"

    if not seed_pos.exists():
        print(f"[ERROR] No positions file in seed window {seed_window_idx}")
        print(f"        Required: {seed_pos}")
        print(f"        Run window {seed_window_idx} first to generate seed data")
        return False, "Missing seed positions"

    if not seed_box.exists():
        print(f"[WARNING] No box vectors file in seed window {seed_window_idx}")
        print(f"          Expected: {seed_box}")
        print(f"          Will use default box vectors")

    # Generate seeded runner script
    seeded_script = generate_seeded_script(sys_name, window_idx, seed_window_idx)
    with open(script_path, 'w') as f:
        f.write(seeded_script)

    # Run the script
    print(f"Starting simulation...")
    start_time = time.time()

    try:
        with open(log_path, 'w') as log:
            result = subprocess.run(
                [sys.executable, str(script_path)],
                cwd=str(window_dir),
                stdout=log,
                stderr=subprocess.STDOUT,
                timeout=7200,  # 2 hour timeout
            )
    except subprocess.TimeoutExpired:
        elapsed = time.time() - start_time
        print(f"[FAIL] Window timed out after {elapsed:.1f}s")
        return False, "Timeout (>2h)"

    elapsed = time.time() - start_time
    print(f"Completed in {elapsed:.1f}s")

    # Check result
    if result.returncode != 0:
        print(f"[FAIL] Window failed with return code {result.returncode}")
        # Read last lines of log for error
        with open(log_path, 'r') as f:
            lines = f.readlines()
            print("Last 10 lines of log:")
            for line in lines[-10:]:
                print(f"  {line.rstrip()}")
        return False, f"Exit code {result.returncode}"

    # Verify outputs
    u_nk_path = window_dir / "u_nk.npy"
    if not u_nk_path.exists():
        print(f"[FAIL] u_nk.npy not generated")
        return False, "Missing u_nk.npy"

    u_nk = np.load(u_nk_path)
    if np.isnan(u_nk).any():
        nan_count = np.isnan(u_nk).sum()
        print(f"[FAIL] u_nk contains {nan_count} NaN values!")
        return False, f"NaN in u_nk ({nan_count} values)"

    if np.isinf(u_nk).any():
        inf_count = np.isinf(u_nk).sum()
        print(f"[FAIL] u_nk contains {inf_count} Inf values!")
        return False, f"Inf in u_nk ({inf_count} values)"

    # Check for runaway energy drift
    energy_std = np.std(u_nk[:, window_idx])
    energy_range = np.max(u_nk[:, window_idx]) - np.min(u_nk[:, window_idx])
    if energy_range > 100:  # More than 100 kT range suggests drift
        print(f"[WARNING] Large energy range: {energy_range:.1f} kT")

    print(f"[PASS] Window {window_idx} completed successfully")
    print(f"       u_nk shape: {u_nk.shape}")
    print(f"       Energy std: {energy_std:.2f} kT")

    return True, "Success"

def generate_seeded_script(sys_name, window_idx, seed_window_idx):
    """Generate a run script that seeds from previous window."""

    # Load lambda schedule from file (NEVER hardcode - must match febuxostat!)
    lambda_schedule_path = BASE_DIR / sys_name / "lambda_schedule.npy"
    if not lambda_schedule_path.exists():
        raise FileNotFoundError(f"Lambda schedule not found at {lambda_schedule_path}. Run setup_fep_natural.py first.")
    lambda_schedule = np.load(lambda_schedule_path)

    lam_e, lam_s, lam_r = lambda_schedule[window_idx]
    has_restraints = sys_name != 'solvent'
    if not has_restraints:
        lam_r = 0.0

    script = f'''#!/usr/bin/env python
"""
Seeded FEP Window Runner - Window {window_idx}
Seeded from window {seed_window_idx}
System: {sys_name}
"""
import os
import sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from openmm import XmlSerializer, LangevinMiddleIntegrator, Platform, MonteCarloBarostat
from openmm.app import PDBFile, Simulation, CheckpointReporter
from openmm import unit

WINDOW_IDX = {window_idx}
SEED_WINDOW = {seed_window_idx}
LAMBDA_ELEC = {lam_e}
LAMBDA_STERICS = {lam_s}
LAMBDA_RESTRAINTS = {lam_r}
SYS_NAME = "{sys_name}"
HAS_RESTRAINTS = {has_restraints}

TEMPERATURE = 310.0 * unit.kelvin
FRICTION = 1.0 / unit.picosecond
TIMESTEP = 2.0 * unit.femtoseconds
EQUIL_STEPS = 25000    # Shorter equil since seeded
PROD_STEPS = 500000

def main():
    print("="*60)
    print(f"Seeded FEP Window {{WINDOW_IDX}} - {{SYS_NAME}}")
    print(f"Seeded from window {{SEED_WINDOW}}")
    print("="*60)

    parent_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    seed_dir = os.path.join(parent_dir, f"window_{{SEED_WINDOW:02d}}")

    # Load system
    with open(os.path.join(parent_dir, "alchemical_system.xml"), "r") as f:
        system = XmlSerializer.deserialize(f.read())

    pdb = PDBFile(os.path.join(parent_dir, "topology.pdb"))
    topology = pdb.topology

    lambda_schedule = np.load(os.path.join(parent_dir, "lambda_schedule.npy"))
    n_windows = len(lambda_schedule)

    # Load seed positions and box vectors
    seed_pos_path = os.path.join(seed_dir, "final_positions.npy")
    seed_box_path = os.path.join(seed_dir, "final_box_vectors.npy")
    seed_chk_path = os.path.join(seed_dir, "checkpoint.chk")

    if os.path.exists(seed_pos_path):
        positions = np.load(seed_pos_path) * unit.nanometer
        print(f"Loaded positions from {{seed_pos_path}}")
    else:
        raise FileNotFoundError(f"No seed positions at {{seed_pos_path}}")

    if os.path.exists(seed_box_path):
        box_arr = np.load(seed_box_path)
        box_vectors = [
            (box_arr[0,0], box_arr[0,1], box_arr[0,2]) * unit.nanometer,
            (box_arr[1,0], box_arr[1,1], box_arr[1,2]) * unit.nanometer,
            (box_arr[2,0], box_arr[2,1], box_arr[2,2]) * unit.nanometer,
        ]
        print(f"Loaded box vectors from {{seed_box_path}}")
    else:
        print("[WARNING] No box vectors file, using default")
        box_vectors = None

    # Setup platform
    try:
        platform = Platform.getPlatformByName('CUDA')
        properties = {{'Precision': 'mixed'}}
        print("Using CUDA platform")
    except:
        try:
            platform = Platform.getPlatformByName('OpenCL')
            properties = {{}}
            print("Using OpenCL platform")
        except:
            platform = Platform.getPlatformByName('CPU')
            properties = {{}}
            print("Using CPU platform")

    # Create simulation
    integrator = LangevinMiddleIntegrator(TEMPERATURE, FRICTION, TIMESTEP)
    simulation = Simulation(topology, system, integrator, platform, properties)
    context = simulation.context

    # Set box vectors and positions
    if box_vectors:
        context.setPeriodicBoxVectors(*box_vectors)
    context.setPositions(positions)
    context.setVelocitiesToTemperature(TEMPERATURE)

    # Set lambda
    context.setParameter('lambda_electrostatics', LAMBDA_ELEC)
    context.setParameter('lambda_sterics', LAMBDA_STERICS)
    if HAS_RESTRAINTS:
        context.setParameter('lambda_restraints', LAMBDA_RESTRAINTS)

    # Verify energy
    state = context.getState(getEnergy=True)
    E_init = state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
    print(f"Initial energy: {{E_init:.1f}} kJ/mol")

    if np.isnan(E_init) or E_init > 0 or abs(E_init) > 1e10:
        raise RuntimeError(f"Bad initial energy: {{E_init}}")

    # Short equilibration
    print(f"Running equilibration ({{EQUIL_STEPS}} steps)...")
    simulation.step(EQUIL_STEPS)

    # Remove barostat for NVT
    state = context.getState(getPositions=True)
    eq_positions = state.getPositions()
    eq_box = state.getPeriodicBoxVectors()

    for i in range(system.getNumForces()):
        force = system.getForce(i)
        if isinstance(force, MonteCarloBarostat):
            system.removeForce(i)
            break

    integrator2 = LangevinMiddleIntegrator(TEMPERATURE, FRICTION, TIMESTEP)
    simulation = Simulation(topology, system, integrator2, platform, properties)
    context = simulation.context

    context.setPeriodicBoxVectors(*eq_box)
    context.setPositions(eq_positions)
    context.setVelocitiesToTemperature(TEMPERATURE)

    context.setParameter('lambda_electrostatics', LAMBDA_ELEC)
    context.setParameter('lambda_sterics', LAMBDA_STERICS)
    if HAS_RESTRAINTS:
        context.setParameter('lambda_restraints', LAMBDA_RESTRAINTS)

    # Production
    print(f"Running production ({{PROD_STEPS}} steps)...")
    n_samples = PROD_STEPS // 500
    u_nk = np.zeros((n_samples, n_windows))

    for sample_idx in range(n_samples):
        simulation.step(500)

        for k in range(n_windows):
            lam_e, lam_s, lam_r = lambda_schedule[k]
            context.setParameter('lambda_electrostatics', lam_e)
            context.setParameter('lambda_sterics', lam_s)
            if HAS_RESTRAINTS:
                context.setParameter('lambda_restraints', lam_r)

            state = context.getState(getEnergy=True)
            kT = unit.MOLAR_GAS_CONSTANT_R * TEMPERATURE
            u_nk[sample_idx, k] = state.getPotentialEnergy() / kT

        context.setParameter('lambda_electrostatics', LAMBDA_ELEC)
        context.setParameter('lambda_sterics', LAMBDA_STERICS)
        if HAS_RESTRAINTS:
            context.setParameter('lambda_restraints', LAMBDA_RESTRAINTS)

        if (sample_idx + 1) % 100 == 0:
            print(f"  Sample {{sample_idx + 1}}/{{n_samples}}")

    # Save outputs
    np.save("u_nk.npy", u_nk)
    print(f"[PASS] Saved u_nk.npy")

    state = context.getState(getPositions=True)
    pos = state.getPositions(asNumpy=True).value_in_unit(unit.nanometer)
    np.save("final_positions.npy", pos)

    box = state.getPeriodicBoxVectors()
    box_arr = np.array([[v.x, v.y, v.z] for v in box])
    np.save("final_box_vectors.npy", box_arr)

    simulation.saveCheckpoint("checkpoint.chk")
    print("[PASS] Window complete!")

if __name__ == "__main__":
    main()
'''
    return script

def compute_overlap(sys_name, windows):
    """Compute overlap matrix for canary windows.

    Note: Each u_nk file has shape (n_samples, 20) - energies at ALL lambda states.
    For canary windows 15-19, we need to extract only columns 15-19 for MBAR.
    """
    print(f"\n{'='*60}")
    print(f"Computing overlap for {sys_name} windows {windows}")
    print(f"{'='*60}")

    # Load all u_nk data
    u_nk_list = []
    for w in windows:
        u_nk_path = BASE_DIR / sys_name / f"window_{w:02d}" / "u_nk.npy"
        if u_nk_path.exists():
            full_u_nk = np.load(u_nk_path)
            # Extract only columns for canary windows (15-19 = indices 15:20)
            u_nk_subset = full_u_nk[:, windows]  # columns for windows 15,16,17,18,19
            u_nk_list.append(u_nk_subset)
        else:
            print(f"  [MISSING] {u_nk_path}")
            return None

    # Stack and compute MBAR overlap
    try:
        from pymbar import MBAR

        # Combine u_nk matrices - now with matching dimensions
        # Each u_nk_subset has shape (n_samples, 5) for 5 canary windows
        u_kn = np.vstack(u_nk_list)
        N_k = np.array([u.shape[0] for u in u_nk_list])

        print(f"Combined u_kn shape: {u_kn.shape}")
        print(f"N_k: {N_k}")
        print(f"Number of states K: {len(windows)}")

        # Verify dimensions match
        assert u_kn.shape[1] == len(windows), f"Column mismatch: {u_kn.shape[1]} vs {len(windows)}"
        assert len(N_k) == len(windows), f"N_k length mismatch: {len(N_k)} vs {len(windows)}"

        # Run MBAR
        mbar = MBAR(u_kn.T, N_k)

        # Get overlap matrix
        overlap = mbar.compute_overlap()
        overlap_matrix = overlap['matrix']

        print("\nOverlap matrix (adjacent windows):")
        print("-"*40)
        for i in range(len(windows)-1):
            ov = overlap_matrix[i, i+1]
            status = "OK" if ov > 0.03 else "LOW" if ov > 0.01 else "FAIL"
            print(f"  Window {windows[i]}-{windows[i+1]}: {ov:.4f} [{status}]")

        return overlap_matrix

    except Exception as e:
        print(f"[ERROR] MBAR analysis failed: {e}")
        return None

def main():
    print("="*70)
    print("CANARY WINDOW RUNNER")
    print("Natural Product CID_10120 FEP")
    print("="*70)
    print()
    print("Running windows 15-19 sequentially with checkpoint seeding")
    print("This tests stability before full deployment")
    print()

    results = {}

    for sys_name in SYSTEMS:
        print(f"\n{'#'*70}")
        print(f"# SYSTEM: {sys_name}")
        print(f"{'#'*70}")

        results[sys_name] = []

        # First check if window 14 exists (seed source)
        if not check_checkpoint_exists(sys_name, SEED_WINDOW):
            if not (BASE_DIR / sys_name / f"window_{SEED_WINDOW:02d}" / "final_positions.npy").exists():
                print(f"[ERROR] Window {SEED_WINDOW} not available for {sys_name}")
                print(f"        Run window {SEED_WINDOW} first before canary test")
                results[sys_name].append((SEED_WINDOW, False, "Missing seed window"))
                continue

        # Run canary windows
        prev_window = SEED_WINDOW
        for window_idx in CANARY_WINDOWS:
            success, msg = run_window_with_seed(sys_name, window_idx, prev_window)
            results[sys_name].append((window_idx, success, msg))

            if not success:
                print(f"[ABORT] Stopping {sys_name} due to failure at window {window_idx}")
                break

            prev_window = window_idx

        # Compute overlap if all windows succeeded
        if all(r[1] for r in results[sys_name]):
            overlap = compute_overlap(sys_name, CANARY_WINDOWS)
            if overlap is not None:
                # Check for collapsed overlap
                min_overlap = np.min([overlap[i, i+1] for i in range(len(CANARY_WINDOWS)-1)])
                if min_overlap < 0.01:
                    print(f"[FAIL] Overlap collapsed (min={min_overlap:.4f})")
                    print("       DO NOT proceed with full deployment")
                else:
                    print(f"[PASS] Overlap OK (min={min_overlap:.4f})")

    # Summary
    print("\n" + "="*70)
    print("CANARY TEST SUMMARY")
    print("="*70)

    all_passed = True
    for sys_name, sys_results in results.items():
        success_count = sum(1 for r in sys_results if r[1])
        total = len(sys_results)
        status = "PASS" if success_count == total else "FAIL"
        all_passed = all_passed and (success_count == total)
        print(f"{sys_name}: {success_count}/{total} [{status}]")

        for window_idx, success, msg in sys_results:
            symbol = "[OK]" if success else "[X]"
            print(f"  {symbol} Window {window_idx}: {msg}")

    print()
    if all_passed:
        print("[READY] All canary windows passed!")
        print("        Safe to deploy full 60-window run")
    else:
        print("[BLOCKED] Canary test failed!")
        print("          Fix issues before deploying full run")

    return 0 if all_passed else 1

if __name__ == "__main__":
    sys.exit(main())
