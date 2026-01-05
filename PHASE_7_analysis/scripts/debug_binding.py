#!/usr/bin/env python3
"""Debug ligand binding - check if ligand stays bound to protein."""

import MDAnalysis as mda
from MDAnalysis.analysis import rms, align
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

TRAJ_DIR = Path(r'C:\Users\vasud\nod2-screening-data\PHASE_6\trajectories')

# Load budesonide
print("Loading budesonide rep1...")
u = mda.Universe(str(TRAJ_DIR / 'budesonide_rep1_solvated.pdb'),
                 str(TRAJ_DIR / 'budesonide_rep1.dcd'))

protein = u.select_atoms("protein")
protein_bb = u.select_atoms("protein and backbone")
ligand = u.select_atoms("resid 1 and resname UNK")
ligand_heavy = u.select_atoms("(resid 1 and resname UNK) and not name H*")

print(f"Protein atoms: {len(protein)}")
print(f"Ligand atoms: {len(ligand)}, heavy: {len(ligand_heavy)}")

# Track positions over trajectory
times = []
prot_com = []
lig_com = []
lig_prot_dist = []  # Distance from ligand to protein COM
box_sizes = []

for ts in u.trajectory[::10]:  # Every 10th frame
    times.append(ts.time / 1000)  # ps to ns

    # COM positions
    pc = protein.center_of_mass()
    lc = ligand.center_of_mass()

    prot_com.append(pc)
    lig_com.append(lc)

    # Distance (raw, without PBC correction)
    dist = np.linalg.norm(lc - pc)
    lig_prot_dist.append(dist)

    box_sizes.append(ts.dimensions[:3] if ts.dimensions is not None else [100,100,100])

times = np.array(times)
prot_com = np.array(prot_com)
lig_com = np.array(lig_com)
lig_prot_dist = np.array(lig_prot_dist)
box_sizes = np.array(box_sizes)

# Plot results
fig, axes = plt.subplots(3, 1, figsize=(12, 10))

# Plot 1: Protein and Ligand COM x-coordinate
ax1 = axes[0]
ax1.plot(times, prot_com[:, 0], label='Protein COM X')
ax1.plot(times, lig_com[:, 0], label='Ligand COM X')
ax1.axhline(y=box_sizes[0, 0], color='gray', linestyle='--', alpha=0.5, label='Box X')
ax1.axhline(y=0, color='gray', linestyle='--', alpha=0.5)
ax1.set_xlabel('Time (ns)')
ax1.set_ylabel('X coordinate (Ang)')
ax1.legend()
ax1.set_title('COM X-coordinate Over Time')

# Plot 2: Distance between ligand and protein COM
ax2 = axes[1]
ax2.plot(times, lig_prot_dist)
ax2.axhline(y=box_sizes[0, 0]/2, color='red', linestyle='--', alpha=0.5, label='Half box')
ax2.set_xlabel('Time (ns)')
ax2.set_ylabel('Distance (Ang)')
ax2.legend()
ax2.set_title('Ligand-Protein COM Distance')

# Plot 3: Detect jumps in ligand position
lig_displacement = np.diff(lig_com, axis=0)
lig_jump_mag = np.linalg.norm(lig_displacement, axis=1)
ax3 = axes[2]
ax3.plot(times[1:], lig_jump_mag)
ax3.axhline(y=20, color='red', linestyle='--', alpha=0.5, label='20 Ang threshold')
ax3.set_xlabel('Time (ns)')
ax3.set_ylabel('Frame-to-frame displacement (Ang)')
ax3.legend()
ax3.set_title('Ligand Frame-to-Frame Displacement (PBC jumps)')

plt.tight_layout()
plt.savefig(Path(__file__).parent.parent / 'debug_binding.png', dpi=150)
print(f"Saved debug plot")

# Summary stats
print(f"\n=== SUMMARY ===")
print(f"Protein COM range: X={prot_com[:,0].min():.1f}-{prot_com[:,0].max():.1f}")
print(f"Ligand COM range: X={lig_com[:,0].min():.1f}-{lig_com[:,0].max():.1f}")
print(f"Lig-Prot distance: {lig_prot_dist.mean():.1f} +/- {lig_prot_dist.std():.1f} Ang")
print(f"Max frame-to-frame ligand jump: {lig_jump_mag.max():.1f} Ang")
print(f"Frames with jump > 20 Ang: {np.sum(lig_jump_mag > 20)}")
print(f"Box size: {box_sizes[0]}")
