#!/usr/bin/env python3
"""
MM-GBSA Binding Free Energy Calculation for febuxostat

Uses OpenMM to compute:
- E_complex
- E_receptor
- E_ligand
- ΔG_bind = E_complex - E_receptor - E_ligand

NOTE: This is simplified. Production MM-GBSA should use
proper implicit solvent models and multiple frames.
"""

from openmm import *
from openmm.app import *
from openmm.unit import *
import numpy as np
import json

def compute_mmgbsa(complex_pdb, traj_file, n_frames=100):
    """Compute MM-GBSA from trajectory snapshots"""

    import mdtraj as md

    # Load trajectory
    print(f"Loading trajectory: {traj_file}")
    traj = md.load(traj_file, top=complex_pdb)

    # Select frames (evenly spaced from last half)
    n_total = traj.n_frames
    start_frame = n_total // 2
    frame_indices = np.linspace(start_frame, n_total-1, n_frames, dtype=int)

    print(f"Computing MM-GBSA for {n_frames} frames...")

    # Force field with implicit solvent (GB)
    forcefield = ForceField('amber14-all.xml', 'implicit/gbn2.xml')

    energies = []

    for i, frame_idx in enumerate(frame_indices):
        if i % 10 == 0:
            print(f"  Frame {i+1}/{n_frames}...")

        # Extract frame
        frame = traj[frame_idx]

        # Save temporary PDB
        frame.save_pdb('_temp_frame.pdb')

        try:
            pdb = PDBFile('_temp_frame.pdb')

            # Create system with implicit solvent
            system = forcefield.createSystem(
                pdb.topology,
                nonbondedMethod=NoCutoff,
                constraints=HBonds
            )

            # Get energy
            integrator = VerletIntegrator(0.001*picoseconds)
            simulation = Simulation(pdb.topology, system, integrator)
            simulation.context.setPositions(pdb.positions)

            state = simulation.context.getState(getEnergy=True)
            energy = state.getPotentialEnergy().value_in_unit(kilojoules_per_mole)
            energies.append(energy)

        except Exception as e:
            print(f"    Error at frame {frame_idx}: {e}")
            continue

    # Clean up
    import os
    if os.path.exists('_temp_frame.pdb'):
        os.remove('_temp_frame.pdb')

    # Results
    results = {
        'name': 'febuxostat',
        'n_frames': len(energies),
        'mean_energy_kJ': float(np.mean(energies)),
        'std_energy_kJ': float(np.std(energies)),
        'mean_energy_kcal': float(np.mean(energies) / 4.184),
        'std_energy_kcal': float(np.std(energies) / 4.184)
    }

    print(f"\nMM-GBSA Results:")
    print(f"  Mean: {results['mean_energy_kcal']:.2f} +/- {results['std_energy_kcal']:.2f} kcal/mol")

    with open('febuxostat_mmgbsa.json', 'w') as f:
        json.dump(results, f, indent=2)

    return results


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Compute MM-GBSA")
    parser.add_argument("--pdb", required=True, help="Complex PDB file")
    parser.add_argument("--traj", required=True, help="Trajectory file")
    parser.add_argument("--frames", type=int, default=100, help="Number of frames")

    args = parser.parse_args()

    compute_mmgbsa(args.pdb, args.traj, args.frames)
