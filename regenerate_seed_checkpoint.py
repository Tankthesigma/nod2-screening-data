#!/usr/bin/env python
"""
Regenerate checkpoint.chk for seed windows that lack box vectors.
This re-runs the window simulation and saves checkpoint + box vectors.
Does NOT overwrite existing u_nk.npy files.

Seed windows needed:
- wt_complex/window_14
- mut_complex/window_14
- mut_complex/window_16
- solvent/window_15
"""
import os
import sys
import numpy as np
import argparse

from openmm import XmlSerializer, LangevinMiddleIntegrator, Platform, MonteCarloBarostat
from openmm.app import PDBFile, Simulation
from openmm import unit

# Simulation parameters (MUST match original)
TEMPERATURE = 310.0 * unit.kelvin
FRICTION = 1.0 / unit.picosecond
TIMESTEP = 2.0 * unit.femtoseconds
EQUIL_STEPS = 50000    # 100 ps NPT equilibration
PROD_STEPS = 500000    # 1 ns production
ENERGY_INTERVAL = 500


def run_seed_window(system_name, window_idx, base_path):
    """Re-run a seed window to generate checkpoint.chk"""

    system_dir = os.path.join(base_path, system_name)
    window_dir = os.path.join(system_dir, f"window_{window_idx:02d}")

    print(f"{'='*60}")
    print(f"Regenerating checkpoint for {system_name}/window_{window_idx:02d}")
    print(f"{'='*60}")

    # Check if u_nk.npy exists (seed must have succeeded before)
    unk_path = os.path.join(window_dir, "u_nk.npy")
    if not os.path.exists(unk_path):
        raise RuntimeError(f"u_nk.npy not found - this window didn't succeed: {unk_path}")

    # Backup existing u_nk.npy
    unk_backup = os.path.join(window_dir, "u_nk_backup.npy")
    if not os.path.exists(unk_backup):
        import shutil
        shutil.copy(unk_path, unk_backup)
        print(f"Backed up u_nk.npy to u_nk_backup.npy")

    # Load lambda schedule
    lambda_schedule = np.load(os.path.join(system_dir, "lambda_schedule.npy"))
    lam_e, lam_s, lam_r = lambda_schedule[window_idx]
    print(f"Lambda: elec={lam_e:.4f}, sterics={lam_s:.4f}, restraints={lam_r:.4f}")

    has_restraints = system_name != "solvent"

    # Load system
    with open(os.path.join(system_dir, "alchemical_system.xml"), "r") as f:
        system = XmlSerializer.deserialize(f.read())

    # Load initial positions
    positions = np.load(os.path.join(system_dir, "positions.npy"))
    positions = positions * unit.nanometer

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
    simulation.context.setPositions(positions)

    # Initialize velocities before NPT (avoid cold start)
    simulation.context.setVelocitiesToTemperature(TEMPERATURE)

    # Set lambda parameters
    context = simulation.context
    context.setParameter('lambda_electrostatics', float(lam_e))
    context.setParameter('lambda_sterics', float(lam_s))
    if has_restraints:
        context.setParameter('lambda_restraints', float(lam_r))

    # Minimize
    print("Minimizing energy...")
    simulation.minimizeEnergy(maxIterations=1000)

    # NPT Equilibration
    print(f"Running NPT equilibration ({EQUIL_STEPS} steps)...")
    simulation.step(EQUIL_STEPS)

    # Get state after equilibration
    state = context.getState(getPositions=True)
    positions_after_equil = state.getPositions()
    box_vectors = state.getPeriodicBoxVectors()

    # Remove barostat for NVT
    print("Removing barostat for NVT production...")
    for i in range(system.getNumForces()):
        force = system.getForce(i)
        if isinstance(force, MonteCarloBarostat):
            system.removeForce(i)
            break

    # Create fresh simulation
    integrator2 = LangevinMiddleIntegrator(TEMPERATURE, FRICTION, TIMESTEP)
    simulation = Simulation(topology, system, integrator2, platform, properties)
    context = simulation.context

    # Restore state
    context.setPeriodicBoxVectors(*box_vectors)
    context.setPositions(positions_after_equil)
    context.setVelocitiesToTemperature(TEMPERATURE)

    # Set lambda parameters
    context.setParameter('lambda_electrostatics', float(lam_e))
    context.setParameter('lambda_sterics', float(lam_s))
    if has_restraints:
        context.setParameter('lambda_restraints', float(lam_r))

    # Verify energy
    state = context.getState(getEnergy=True)
    energy = state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
    print(f"NVT initial energy: {energy:.1f} kJ/mol")

    if not np.isfinite(energy) or energy > 0 or abs(energy) > 1e10:
        raise RuntimeError(f"Energy looks wrong: {energy}")

    # Run production (we don't need to recompute u_nk, just get to final state)
    print(f"Running NVT production ({PROD_STEPS} steps)...")
    n_samples = PROD_STEPS // ENERGY_INTERVAL

    for sample_idx in range(n_samples):
        simulation.step(ENERGY_INTERVAL)
        if (sample_idx + 1) % 200 == 0:
            print(f"  Progress: {sample_idx + 1}/{n_samples}")

    # Save checkpoint (contains positions, velocities, box vectors)
    checkpoint_path = os.path.join(window_dir, "checkpoint.chk")
    simulation.saveCheckpoint(checkpoint_path)
    print(f"[PASS] Saved checkpoint.chk")

    # Also save box vectors separately for verification
    state = context.getState(getPositions=True)
    box = state.getPeriodicBoxVectors()
    # Convert to nm (plain floats) for np.save
    box_array = np.array([
        [box[0].x / unit.nanometer, box[0].y / unit.nanometer, box[0].z / unit.nanometer],
        [box[1].x / unit.nanometer, box[1].y / unit.nanometer, box[1].z / unit.nanometer],
        [box[2].x / unit.nanometer, box[2].y / unit.nanometer, box[2].z / unit.nanometer]
    ], dtype=np.float64)
    box_path = os.path.join(window_dir, "final_box_vectors.npy")
    np.save(box_path, box_array)
    print(f"[PASS] Saved final_box_vectors.npy: {box_array[0,0]:.3f} x {box_array[1,1]:.3f} x {box_array[2,2]:.3f} nm")

    # Save final positions (overwrite with consistent state)
    pos = state.getPositions(asNumpy=True).value_in_unit(unit.nanometer)
    np.save(os.path.join(window_dir, "final_positions.npy"), pos)
    print(f"[PASS] Saved final_positions.npy")

    # Restore original u_nk.npy from backup
    import shutil
    shutil.copy(unk_backup, unk_path)
    print(f"[PASS] Restored original u_nk.npy from backup")

    print(f"\n{'='*60}")
    print(f"Seed window {system_name}/window_{window_idx:02d} complete!")
    print(f"{'='*60}\n")


def main():
    parser = argparse.ArgumentParser(description="Regenerate seed window checkpoints")
    parser.add_argument("--system", required=True, help="System name (wt_complex, mut_complex, solvent)")
    parser.add_argument("--window", type=int, required=True, help="Window index")
    parser.add_argument("--base", default="C:/Users/vasud/nod2-screening-data/fep_pmx",
                        help="Base path to fep_pmx directory")

    args = parser.parse_args()
    run_seed_window(args.system, args.window, args.base)


if __name__ == "__main__":
    main()
