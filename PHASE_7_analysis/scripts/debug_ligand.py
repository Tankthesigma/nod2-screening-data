#!/usr/bin/env python3
"""Debug ligand selection and RMSD issues."""

import MDAnalysis as mda
from pathlib import Path
import numpy as np

TRAJ_DIR = Path(r'C:\Users\vasud\nod2-screening-data\PHASE_6\trajectories')

# Check budesonide
print("Loading budesonide...")
u = mda.Universe(str(TRAJ_DIR / 'budesonide_rep1_solvated.pdb'), str(TRAJ_DIR / 'budesonide_rep1.dcd'))

print('=== BUDESONIDE ===')
print(f'Total atoms: {len(u.atoms)}')

protein = u.select_atoms("protein")
print(f'Protein atoms: {len(protein)}')

# Check non-protein residues
nonprot = u.select_atoms('not protein')
print(f'Non-protein atoms: {len(nonprot)}')

# Group by resname
from collections import Counter
resnames = Counter(r.resname for r in nonprot.residues)
print(f'Non-protein residues by type: {dict(resnames)}')

# Find the ligand candidate
WATER = {'HOH', 'WAT', 'SOL', 'TIP3', 'TIP3P'}
IONS = {'NA', 'CL', 'NA+', 'CL-'}

candidates = []
for res in nonprot.residues:
    rname = res.resname.upper()
    if rname not in WATER and rname not in IONS:
        if len(res.atoms) >= 5:
            candidates.append((res.resname, res.resid, len(res.atoms)))

print(f'Ligand candidates (>=5 atoms): {candidates}')

# Check ligand positions at frame 0 vs last frame
if candidates:
    lig_resname = candidates[0][0]
    lig = u.select_atoms(f'resname {lig_resname}')
    print(f'\nLigand resname: {lig_resname}, atoms: {len(lig)}')

    u.trajectory[0]
    com0 = lig.center_of_mass()
    print(f'Frame 0 ligand COM: {com0}')

    u.trajectory[-1]
    com_last = lig.center_of_mass()
    print(f'Last frame ligand COM: {com_last}')

    # Distance moved
    dist = np.linalg.norm(com_last - com0)
    print(f'Distance COM moved: {dist:.2f} Angstrom')

    # Box dimensions
    print(f'\nBox dimensions: {u.trajectory[0].dimensions[:3]}')

    # Check if ligand crossed PBC
    box = u.trajectory[0].dimensions[:3]
    if any(abs(com_last[i] - com0[i]) > box[i]/2 for i in range(3)):
        print("WARNING: Ligand may have crossed periodic boundary!")

# Check protein COM movement
u.trajectory[0]
prot_com0 = protein.center_of_mass()
u.trajectory[-1]
prot_com_last = protein.center_of_mass()
prot_dist = np.linalg.norm(prot_com_last - prot_com0)
print(f'\nProtein COM moved: {prot_dist:.2f} Angstrom')
