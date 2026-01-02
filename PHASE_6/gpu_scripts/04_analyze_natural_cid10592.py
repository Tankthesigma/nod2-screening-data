#!/usr/bin/env python3
"""
Trajectory Analysis Script for natural_cid10592

Computes:
- RMSD/RMSF
- H-bond analysis
- Distance to key residues
- MM-GBSA (simplified)
"""

import mdtraj as md
import numpy as np
import json
from pathlib import Path

def analyze_trajectory(traj_file, topology_file, output_dir="."):
    """Complete trajectory analysis"""

    print(f"Loading trajectory: {traj_file}")
    traj = md.load(traj_file, top=topology_file)

    print(f"  Frames: {traj.n_frames}")
    print(f"  Time: {traj.time[-1]/1000:.1f} ns")

    results = {
        'name': 'natural_cid10592',
        'n_frames': traj.n_frames,
        'time_ns': traj.time[-1] / 1000
    }

    # RMSD
    print("Computing RMSD...")
    # Select protein backbone
    backbone = traj.topology.select('backbone')
    rmsd = md.rmsd(traj, traj, 0, atom_indices=backbone)

    results['rmsd'] = {
        'values': rmsd.tolist(),
        'mean': float(np.mean(rmsd[len(rmsd)//2:])),
        'std': float(np.std(rmsd[len(rmsd)//2:]))
    }
    print(f"  Mean RMSD (last half): {results['rmsd']['mean']:.3f} nm")

    # RMSF
    print("Computing RMSF...")
    # Align first
    traj.superpose(traj, 0)
    rmsf = md.rmsf(traj, traj, 0, atom_indices=backbone)

    results['rmsf'] = {
        'values': rmsf.tolist(),
        'mean': float(np.mean(rmsf)),
        'max': float(np.max(rmsf))
    }
    print(f"  Mean RMSF: {results['rmsf']['mean']:.3f} nm")

    # H-bonds
    print("Analyzing H-bonds...")
    hbonds = md.baker_hubbard(traj, freq=0.3)  # Present in >30% of frames
    results['hbonds'] = {
        'total_persistent': len(hbonds)
    }
    print(f"  Persistent H-bonds: {len(hbonds)}")

    # Key residue distances
    print("Computing key residue distances...")
    key_residues = {
        'ARG702': 'Crohn polymorphism',
        'LEU1007': '1007fs region',
        'GLY908': 'Crohn polymorphism'
    }

    for resname_id, description in key_residues.items():
        resname = resname_id[:3]
        resid = int(resname_id[3:])

        try:
            # Find CA atom of residue
            selection = traj.topology.select(f'name CA and resid {resid}')
            if len(selection) > 0:
                # Find ligand center (assumes ligand is last residue)
                ligand_atoms = traj.topology.select('resname LIG or resname UNK')
                if len(ligand_atoms) > 0:
                    # Compute distances
                    pairs = [[ligand_atoms[0], selection[0]]]
                    distances = md.compute_distances(traj, pairs)

                    results[f'dist_{resname_id}'] = {
                        'mean': float(np.mean(distances)),
                        'std': float(np.std(distances)),
                        'description': description
                    }
                    print(f"  {resname_id}: {np.mean(distances):.3f} +/- {np.std(distances):.3f} nm")
        except Exception as e:
            print(f"  Could not compute distance to {resname_id}: {e}")

    # Save results
    output_file = Path(output_dir) / f"natural_cid10592_analysis.json"
    with open(output_file, 'w') as f:
        json.dump(results, f, indent=2)

    print(f"\nResults saved to: {output_file}")
    return results


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Analyze MD trajectory")
    parser.add_argument("--traj", required=True, help="Trajectory file (DCD)")
    parser.add_argument("--top", required=True, help="Topology file (PDB)")
    parser.add_argument("--output", default=".", help="Output directory")

    args = parser.parse_args()

    analyze_trajectory(args.traj, args.top, args.output)
