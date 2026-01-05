#!/usr/bin/env python3
"""Debug RMSD calculation - full stride."""

import MDAnalysis as mda
from MDAnalysis.analysis import rms
from pathlib import Path
import numpy as np

TRAJ_DIR = Path(r'C:\Users\vasud\nod2-screening-data\PHASE_6\trajectories')

print("Loading budesonide rep1...")
u = mda.Universe(str(TRAJ_DIR / 'budesonide_rep1_solvated.pdb'), str(TRAJ_DIR / 'budesonide_rep1.dcd'))

# Ligand selection - same as main script
lig_sel = "resid 1 and resname UNK"
lig_heavy_sel = f"({lig_sel}) and not name H*"

ligand_heavy = u.select_atoms(lig_heavy_sel)
print(f"Ligand heavy atoms: {len(ligand_heavy)}")

# Run with stride=10 like the main script
print("\n=== RMSD with stride=10 (like main script) ===")
R = rms.RMSD(
    u, u,
    select="protein and backbone",
    groupselections=[lig_heavy_sel],
    ref_frame=0
)
R.run(step=10, verbose=True)

print(f"\nTotal frames analyzed: {len(R.results.rmsd)}")
print(f"Backbone RMSD: mean={R.results.rmsd[:,2].mean():.2f}, std={R.results.rmsd[:,2].std():.2f}")
print(f"Ligand RMSD: mean={R.results.rmsd[:,3].mean():.2f}, std={R.results.rmsd[:,3].std():.2f}")

# Check if there are huge jumps
lig_rmsd = R.results.rmsd[:,3]
print(f"\nLigand RMSD range: {lig_rmsd.min():.2f} - {lig_rmsd.max():.2f}")
print(f"Frames with RMSD > 20Å: {np.sum(lig_rmsd > 20)}")
print(f"Frames with RMSD > 50Å: {np.sum(lig_rmsd > 50)}")

# Plot distribution
import matplotlib.pyplot as plt
plt.figure(figsize=(10,4))
plt.subplot(1,2,1)
plt.plot(R.results.rmsd[:,1]/1000, R.results.rmsd[:,3])
plt.xlabel('Time (ns)')
plt.ylabel('Ligand RMSD (Å)')
plt.title('Budesonide Ligand RMSD over time')

plt.subplot(1,2,2)
plt.hist(lig_rmsd, bins=50)
plt.xlabel('Ligand RMSD (Å)')
plt.ylabel('Count')
plt.title('Distribution')

plt.tight_layout()
plt.savefig(r'C:\Users\vasud\nod2-screening-data\PHASE_7_analysis\debug_rmsd_budesonide.png', dpi=150)
print("\nSaved debug plot to debug_rmsd_budesonide.png")
