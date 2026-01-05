#!/usr/bin/env python3
"""Debug RMSD calculation."""

import MDAnalysis as mda
from MDAnalysis.analysis import rms
from pathlib import Path
import numpy as np

TRAJ_DIR = Path(r'C:\Users\vasud\nod2-screening-data\PHASE_6\trajectories')

print("Loading budesonide...")
u = mda.Universe(str(TRAJ_DIR / 'budesonide_rep1_solvated.pdb'), str(TRAJ_DIR / 'budesonide_rep1.dcd'))

# Ligand selection
lig_sel = "resid 1 and resname UNK"
lig_heavy_sel = f"({lig_sel}) and not name H*"

ligand = u.select_atoms(lig_sel)
ligand_heavy = u.select_atoms(lig_heavy_sel)

print(f"Ligand atoms: {len(ligand)}")
print(f"Ligand heavy atoms: {len(ligand_heavy)}")

if len(ligand_heavy) == 0:
    print("WARNING: No heavy atoms selected! Trying without H filter...")
    lig_heavy_sel = lig_sel
    ligand_heavy = u.select_atoms(lig_heavy_sel)
    print(f"Ligand atoms (no filter): {len(ligand_heavy)}")

# Print atom names to see if H filter is wrong
print(f"\nFirst 10 ligand atom names: {[a.name for a in ligand.atoms[:10]]}")

# Test RMSD calculation
print("\n=== Testing rms.RMSD ===")
try:
    R = rms.RMSD(
        u, u,
        select="protein and backbone",
        groupselections=[lig_heavy_sel],
        ref_frame=0
    )
    R.run(step=100, verbose=True)  # Only 20 frames for testing

    print(f"\nRMSD results shape: {R.results.rmsd.shape}")
    print(f"Columns: Frame, Time, Backbone RMSD, Ligand RMSD")
    print(f"First 5 rows:\n{R.results.rmsd[:5]}")
    print(f"Last 5 rows:\n{R.results.rmsd[-5:]}")

    print(f"\nBackbone RMSD: mean={R.results.rmsd[:,2].mean():.2f}, max={R.results.rmsd[:,2].max():.2f}")
    print(f"Ligand RMSD: mean={R.results.rmsd[:,3].mean():.2f}, max={R.results.rmsd[:,3].max():.2f}")

except Exception as e:
    print(f"ERROR: {e}")
    import traceback
    traceback.print_exc()

# Manual RMSD calculation for comparison
print("\n=== Manual RMSD check ===")
u.trajectory[0]
ref_lig_pos = ligand_heavy.positions.copy()
ref_lig_com = ref_lig_pos.mean(axis=0)
print(f"Frame 0 ligand COM: {ref_lig_com}")

u.trajectory[-1]
cur_lig_pos = ligand_heavy.positions.copy()
cur_lig_com = cur_lig_pos.mean(axis=0)
print(f"Last frame ligand COM: {cur_lig_com}")

# Simple RMSD without alignment
diff = cur_lig_pos - ref_lig_pos
simple_rmsd = np.sqrt((diff**2).sum() / len(diff))
print(f"Simple RMSD (no alignment): {simple_rmsd:.2f}")
