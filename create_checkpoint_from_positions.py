#!/usr/bin/env python
"""
Create checkpoint.chk from existing final_positions.npy WITHOUT re-running simulation.

The seed windows already succeeded and have final_positions.npy.
We just need to create a checkpoint file for seeding the next window.

Box vectors are extracted from topology.pdb CRYST1 record (NPT should have
relaxed close to this for a well-equilibrated system).
"""
import os
import sys
import numpy as np
import argparse
import re

from openmm import XmlSerializer, LangevinMiddleIntegrator, Platform, MonteCarloBarostat, Vec3
from openmm.app import PDBFile, Simulation
from openmm import unit

# Simulation parameters (MUST match original)
TEMPERATURE = 310.0 * unit.kelvin
FRICTION = 1.0 / unit.picosecond
TIMESTEP = 2.0 * unit.femtoseconds


def get_box_from_pdb(pdb_path):
    """Extract box vectors from CRYST1 record in PDB file."""
    with open(pdb_path, 'r') as f:
        for line in f:
            if line.startswith('CRYST1'):
                # CRYST1    a      b      c    alpha  beta  gamma
                parts = line.split()
                a = float(parts[1]) / 10.0  # Angstrom to nm
                b = float(parts[2]) / 10.0
                c = float(parts[3]) / 10.0
                # Assuming orthorhombic box (90 degree angles)
                return np.array([[a, 0, 0], [0, b, 0], [0, 0, c]], dtype=np.float64)
    raise RuntimeError(f"No CRYST1 record found in {pdb_path}")


def create_checkpoint(system_name, window_idx, base_path):
    """Create checkpoint from existing final_positions.npy."""

    system_dir = os.path.join(base_path, system_name)
    window_dir = os.path.join(system_dir, f"window_{window_idx:02d}")

    print(f"{'='*60}")
    print(f"Creating checkpoint for {system_name}/window_{window_idx:02d}")
    print(f"{'='*60}")

    # Check final_positions.npy exists
    pos_path = os.path.join(window_dir, "final_positions.npy")
    if not os.path.exists(pos_path):
        raise RuntimeError(f"final_positions.npy not found: {pos_path}")

    # Load positions (already in nm as plain floats)
    positions = np.load(pos_path)
    positions = positions * unit.nanometer
    print(f"Loaded positions: {len(positions)} atoms")

    # Get box vectors from topology.pdb
    pdb_path = os.path.join(system_dir, "topology.pdb")
    box_nm = get_box_from_pdb(pdb_path)
    print(f"Box from PDB: {box_nm[0,0]:.3f} x {box_nm[1,1]:.3f} x {box_nm[2,2]:.3f} nm")

    # Load lambda schedule
    lambda_schedule = np.load(os.path.join(system_dir, "lambda_schedule.npy"))
    lam_e, lam_s, lam_r = lambda_schedule[window_idx]
    print(f"Lambda: elec={lam_e:.4f}, sterics={lam_s:.4f}, restraints={lam_r:.4f}")

    has_restraints = system_name != "solvent"

    # Load system (need to remove barostat for NVT checkpoint)
    with open(os.path.join(system_dir, "alchemical_system.xml"), "r") as f:
        system = XmlSerializer.deserialize(f.read())

    # Remove barostat (we're creating NVT checkpoint)
    for i in range(system.getNumForces()):
        force = system.getForce(i)
        if isinstance(force, MonteCarloBarostat):
            system.removeForce(i)
            print("Removed barostat for NVT checkpoint")
            break

    # Load topology
    pdb = PDBFile(pdb_path)
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

    # Set box vectors FIRST (use Vec3 with units)
    box_vectors = (
        Vec3(box_nm[0, 0], box_nm[0, 1], box_nm[0, 2]) * unit.nanometer,
        Vec3(box_nm[1, 0], box_nm[1, 1], box_nm[1, 2]) * unit.nanometer,
        Vec3(box_nm[2, 0], box_nm[2, 1], box_nm[2, 2]) * unit.nanometer
    )
    context.setPeriodicBoxVectors(*box_vectors)

    # Set positions
    context.setPositions(positions)

    # Set lambda parameters (with guard for parameter existence)
    params = context.getParameters()
    if 'lambda_electrostatics' in params:
        context.setParameter('lambda_electrostatics', float(lam_e))
    if 'lambda_sterics' in params:
        context.setParameter('lambda_sterics', float(lam_s))
    if has_restraints and 'lambda_restraints' in params:
        context.setParameter('lambda_restraints', float(lam_r))

    # Initialize velocities
    context.setVelocitiesToTemperature(TEMPERATURE)

    # Aggressive minimization to fix box mismatch clashes
    print("Minimizing energy (up to 5000 steps)...")
    simulation.minimizeEnergy(maxIterations=5000)

    # Check energy
    state = context.getState(getEnergy=True)
    energy = state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
    print(f"Energy after minimization: {energy:.1f} kJ/mol")

    if not np.isfinite(energy):
        raise RuntimeError(f"Energy is NaN/Inf: {energy}")

    # Skip dynamics - just save checkpoint with minimized positions
    # This is safer for windows with box vector mismatch
    print("Skipping dynamics - saving checkpoint with minimized positions")

    # Save checkpoint
    checkpoint_path = os.path.join(window_dir, "checkpoint.chk")
    simulation.saveCheckpoint(checkpoint_path)
    print(f"[PASS] Saved checkpoint.chk")

    # Save box vectors
    state = context.getState(getPositions=True)
    box = state.getPeriodicBoxVectors()
    # Extract values properly - box vectors are Vec3 with units
    box_array = np.array([
        [box[0][0].value_in_unit(unit.nanometer), box[0][1].value_in_unit(unit.nanometer), box[0][2].value_in_unit(unit.nanometer)],
        [box[1][0].value_in_unit(unit.nanometer), box[1][1].value_in_unit(unit.nanometer), box[1][2].value_in_unit(unit.nanometer)],
        [box[2][0].value_in_unit(unit.nanometer), box[2][1].value_in_unit(unit.nanometer), box[2][2].value_in_unit(unit.nanometer)]
    ], dtype=np.float64)
    box_path = os.path.join(window_dir, "final_box_vectors.npy")
    np.save(box_path, box_array)
    print(f"[PASS] Saved final_box_vectors.npy: {box_array[0,0]:.3f} x {box_array[1,1]:.3f} x {box_array[2,2]:.3f} nm")

    print(f"\n{'='*60}")
    print(f"Checkpoint created for {system_name}/window_{window_idx:02d}!")
    print(f"{'='*60}\n")


def main():
    parser = argparse.ArgumentParser(description="Create checkpoint from existing positions")
    parser.add_argument("--system", required=True, help="System name")
    parser.add_argument("--window", type=int, required=True, help="Window index")
    parser.add_argument("--base", default="C:/Users/vasud/nod2-screening-data/fep_pmx",
                        help="Base path to fep_pmx directory")

    args = parser.parse_args()
    create_checkpoint(args.system, args.window, args.base)


if __name__ == "__main__":
    main()
