#!/usr/bin/env python
"""
Select Boresch restraint anchor atoms for Natural Product CID_10120.

Boresch restraints require 6 atoms:
- 3 protein atoms (P1, P2, P3) - typically backbone CA atoms
- 3 ligand atoms (L1, L2, L3) - typically heavy atoms

The restraints are:
- 1 distance: r(P1-L1)
- 2 angles: theta(P2-P1-L1), theta(P1-L1-L2)
- 3 dihedrals: phi(P3-P2-P1-L1), phi(P2-P1-L1-L2), phi(P1-L1-L2-L3)

Force constants (MUST MATCH FEBUXOSTAT):
- k_distance = 4184.0 kJ/(mol*nm^2)  = 10 kcal/(mol*A^2)
- k_angle = 41.84 kJ/(mol*rad^2)     = 10 kcal/(mol*rad^2)
- k_dihedral = 41.84 kJ/(mol*rad^2)  = 10 kcal/(mol*rad^2)
"""

import numpy as np
from pathlib import Path
import json

try:
    from openmm.app import PDBFile
    from openmm import unit
    OPENMM_AVAILABLE = True
except ImportError:
    OPENMM_AVAILABLE = False
    print("[WARNING] OpenMM not available, using fallback parser")

BASE_DIR = Path("C:/Users/vasud/nod2-screening-data/fep_pmx_natural")
INPUT_PDB = Path("C:/Users/vasud/nod2-screening-data/PHASE_6/structures/complex_natural_cid10592.pdb")

# Force constants (MUST MATCH FEBUXOSTAT)
K_DISTANCE = 4184.0      # kJ/(mol*nm^2)
K_ANGLE = 41.84          # kJ/(mol*rad^2)
K_DIHEDRAL = 41.84       # kJ/(mol*rad^2)

def parse_pdb_simple(pdb_path):
    """Parse PDB file without OpenMM."""
    if not Path(pdb_path).exists():
        raise FileNotFoundError(f"PDB file not found: {pdb_path}")
    atoms = []
    with open(pdb_path, 'r') as f:
        for line in f:
            if line.startswith('ATOM') or line.startswith('HETATM'):
                # Handle PDB insertion codes (e.g., "702A") by stripping non-digits
                resid_str = line[22:26].strip()
                try:
                    resid = int(resid_str)
                except ValueError:
                    # Strip insertion code (last char if it's a letter)
                    resid = int(''.join(c for c in resid_str if c.isdigit()))
                atom = {
                    'serial': int(line[6:11]),
                    'name': line[12:16].strip(),
                    'resname': line[17:20].strip(),
                    'chain': line[21],
                    'resid': resid,
                    'x': float(line[30:38]),
                    'y': float(line[38:46]),
                    'z': float(line[46:54]),
                }
                atoms.append(atom)
    return atoms

def distance(a1, a2):
    """Calculate distance between two atoms."""
    return np.sqrt((a1['x']-a2['x'])**2 + (a1['y']-a2['y'])**2 + (a1['z']-a2['z'])**2)

def find_ligand_atoms(atoms):
    """Find ligand atoms (non-standard residues, excluding water/ions)."""
    standard_residues = {
        'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE',
        'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL',
        'HOH', 'WAT', 'NA', 'CL', 'K', 'MG', 'CA', 'ZN', 'FE',
        'HIE', 'HID', 'HIP', 'CYX', 'ASH', 'GLH',
    }
    ligand_atoms = [a for a in atoms if a['resname'] not in standard_residues]
    return ligand_atoms

def find_protein_ca_atoms(atoms):
    """Find protein CA (alpha carbon) atoms, excluding calcium ions."""
    # Standard protein residues that have CA atoms
    protein_residues = {
        'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE',
        'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL',
        'HIE', 'HID', 'HIP', 'CYX', 'ASH', 'GLH',  # protonation variants
    }
    ca_atoms = [a for a in atoms if a['name'] == 'CA' and a['resname'] in protein_residues]
    return ca_atoms

def select_boresch_anchors(atoms):
    """Select optimal Boresch anchor atoms."""
    print("="*70)
    print("SELECTING BORESCH ANCHOR ATOMS")
    print("="*70)

    # Find ligand and protein atoms
    ligand_atoms = find_ligand_atoms(atoms)
    ca_atoms = find_protein_ca_atoms(atoms)

    print(f"Found {len(ligand_atoms)} ligand atoms")
    print(f"Found {len(ca_atoms)} protein CA atoms")

    if len(ligand_atoms) < 3:
        raise ValueError("Need at least 3 ligand atoms for Boresch restraints")

    # Find ligand center of mass
    lig_com = np.array([
        np.mean([a['x'] for a in ligand_atoms]),
        np.mean([a['y'] for a in ligand_atoms]),
        np.mean([a['z'] for a in ligand_atoms]),
    ])
    print(f"Ligand COM: ({lig_com[0]:.2f}, {lig_com[1]:.2f}, {lig_com[2]:.2f}) A")

    # Find heavy ligand atoms (exclude H)
    heavy_ligand = [a for a in ligand_atoms if not a['name'].startswith('H')]
    print(f"Heavy ligand atoms: {len(heavy_ligand)}")

    # Select L1: heavy atom closest to protein
    min_dist = float('inf')
    l1 = None
    for lig_atom in heavy_ligand:
        for ca in ca_atoms:
            d = distance(lig_atom, ca)
            if d < min_dist:
                min_dist = d
                l1 = lig_atom

    if l1 is None:
        raise ValueError("Could not find L1 anchor")

    print(f"\nL1 (closest to protein): {l1['name']} ({l1['resname']} {l1['resid']})")
    print(f"  Position: ({l1['x']:.2f}, {l1['y']:.2f}, {l1['z']:.2f}) A")

    # Select L2: heavy atom bonded to L1 (within ~1.6 A)
    l2_candidates = []
    for a in heavy_ligand:
        if a['serial'] != l1['serial']:
            d = distance(a, l1)
            if 1.2 < d < 1.8:  # typical bond length
                l2_candidates.append((d, a))

    l2_candidates.sort(key=lambda x: x[0])
    if not l2_candidates:
        # Fallback: just pick closest heavy atom
        l2_candidates = [(distance(a, l1), a) for a in heavy_ligand if a['serial'] != l1['serial']]
        l2_candidates.sort(key=lambda x: x[0])

    if not l2_candidates:
        raise ValueError("Cannot find L2: need at least 2 heavy atoms in ligand")
    l2 = l2_candidates[0][1]
    print(f"L2 (bonded to L1): {l2['name']} ({l2['resname']} {l2['resid']})")

    # Select L3: heavy atom bonded to L2 but not L1
    l3_candidates = []
    for a in heavy_ligand:
        if a['serial'] not in [l1['serial'], l2['serial']]:
            d = distance(a, l2)
            if 1.2 < d < 1.8:
                l3_candidates.append((d, a))

    l3_candidates.sort(key=lambda x: x[0])
    if not l3_candidates:
        l3_candidates = [(distance(a, l2), a) for a in heavy_ligand
                         if a['serial'] not in [l1['serial'], l2['serial']]]
        l3_candidates.sort(key=lambda x: x[0])

    if not l3_candidates:
        raise ValueError("Cannot find L3: need at least 3 heavy atoms in ligand")
    l3 = l3_candidates[0][1]
    print(f"L3 (bonded to L2): {l3['name']} ({l3['resname']} {l3['resid']})")

    # Select P1: CA atom closest to L1 (within reasonable distance)
    if not ca_atoms:
        raise ValueError("No protein CA atoms found - cannot select Boresch anchors")
    p1_candidates = [(distance(ca, l1), ca) for ca in ca_atoms]
    p1_candidates.sort(key=lambda x: x[0])
    p1 = p1_candidates[0][1]
    r0 = p1_candidates[0][0] / 10.0  # Convert to nm

    print(f"\nP1 (closest CA to L1): {p1['name']} ({p1['resname']} {p1['resid']})")
    print(f"  Distance P1-L1: {r0*10:.2f} A = {r0:.4f} nm")

    # Select P2: adjacent CA in sequence
    p1_resid = p1['resid']
    p2_candidates = [ca for ca in ca_atoms if abs(ca['resid'] - p1_resid) == 1]
    if not p2_candidates:
        p2_candidates = [ca for ca in ca_atoms if ca['resid'] != p1_resid]
    if not p2_candidates:
        raise ValueError("Cannot find P2: need at least 2 protein CA atoms")
    p2_candidates.sort(key=lambda x: distance(x, p1))
    p2 = p2_candidates[0]
    print(f"P2 (adjacent to P1): {p2['name']} ({p2['resname']} {p2['resid']})")

    # Select P3: another adjacent CA
    p3_candidates = [ca for ca in ca_atoms
                     if ca['resid'] not in [p1['resid'], p2['resid']]
                     and abs(ca['resid'] - p2['resid']) <= 2]
    if not p3_candidates:
        p3_candidates = [ca for ca in ca_atoms
                         if ca['resid'] not in [p1['resid'], p2['resid']]]
    if not p3_candidates:
        raise ValueError("Cannot find P3: need at least 3 protein CA atoms")
    p3_candidates.sort(key=lambda x: distance(x, p2))
    p3 = p3_candidates[0]
    print(f"P3 (adjacent to P2): {p3['name']} ({p3['resname']} {p3['resid']})")

    # Calculate equilibrium values
    print("\n" + "="*70)
    print("EQUILIBRIUM VALUES")
    print("="*70)

    # Distance r0 (P1-L1)
    print(f"r0 (P1-L1): {r0:.6f} nm")

    # Angles (in radians)
    def angle_3atoms(a1, a2, a3):
        """Calculate angle a1-a2-a3."""
        v1 = np.array([a1['x']-a2['x'], a1['y']-a2['y'], a1['z']-a2['z']])
        v2 = np.array([a3['x']-a2['x'], a3['y']-a2['y'], a3['z']-a2['z']])
        norm1 = np.linalg.norm(v1)
        norm2 = np.linalg.norm(v2)
        if norm1 < 1e-10 or norm2 < 1e-10:
            raise ValueError(f"Degenerate angle: atoms are coincident (norm1={norm1}, norm2={norm2})")
        cos_angle = np.dot(v1, v2) / (norm1 * norm2)
        return np.arccos(np.clip(cos_angle, -1, 1))

    theta1 = angle_3atoms(p2, p1, l1)  # P2-P1-L1
    theta2 = angle_3atoms(p1, l1, l2)  # P1-L1-L2

    print(f"theta1 (P2-P1-L1): {theta1:.6f} rad = {np.degrees(theta1):.2f} deg")
    print(f"theta2 (P1-L1-L2): {theta2:.6f} rad = {np.degrees(theta2):.2f} deg")

    # Dihedrals (in radians)
    def dihedral_4atoms(a1, a2, a3, a4):
        """Calculate dihedral angle a1-a2-a3-a4."""
        b1 = np.array([a2['x']-a1['x'], a2['y']-a1['y'], a2['z']-a1['z']])
        b2 = np.array([a3['x']-a2['x'], a3['y']-a2['y'], a3['z']-a2['z']])
        b3 = np.array([a4['x']-a3['x'], a4['y']-a3['y'], a4['z']-a3['z']])

        norm_b2 = np.linalg.norm(b2)
        if norm_b2 < 1e-10:
            raise ValueError("Degenerate dihedral: central bond has zero length")

        n1 = np.cross(b1, b2)
        n2 = np.cross(b2, b3)

        m1 = np.cross(n1, b2/norm_b2)

        x = np.dot(n1, n2)
        y = np.dot(m1, n2)

        return np.arctan2(y, x)

    phi1 = dihedral_4atoms(p3, p2, p1, l1)  # P3-P2-P1-L1
    phi2 = dihedral_4atoms(p2, p1, l1, l2)  # P2-P1-L1-L2
    phi3 = dihedral_4atoms(p1, l1, l2, l3)  # P1-L1-L2-L3

    print(f"phi1 (P3-P2-P1-L1): {phi1:.6f} rad = {np.degrees(phi1):.2f} deg")
    print(f"phi2 (P2-P1-L1-L2): {phi2:.6f} rad = {np.degrees(phi2):.2f} deg")
    print(f"phi3 (P1-L1-L2-L3): {phi3:.6f} rad = {np.degrees(phi3):.2f} deg")

    # Build anchor dict
    anchors = {
        'P1': {'serial': p1['serial'], 'name': p1['name'], 'resname': p1['resname'], 'resid': p1['resid']},
        'P2': {'serial': p2['serial'], 'name': p2['name'], 'resname': p2['resname'], 'resid': p2['resid']},
        'P3': {'serial': p3['serial'], 'name': p3['name'], 'resname': p3['resname'], 'resid': p3['resid']},
        'L1': {'serial': l1['serial'], 'name': l1['name'], 'resname': l1['resname'], 'resid': l1['resid']},
        'L2': {'serial': l2['serial'], 'name': l2['name'], 'resname': l2['resname'], 'resid': l2['resid']},
        'L3': {'serial': l3['serial'], 'name': l3['name'], 'resname': l3['resname'], 'resid': l3['resid']},
        'equilibrium': {
            'r0_nm': r0,
            'theta1_rad': theta1,
            'theta2_rad': theta2,
            'phi1_rad': phi1,
            'phi2_rad': phi2,
            'phi3_rad': phi3,
        },
        'force_constants': {
            'k_distance': K_DISTANCE,
            'k_angle': K_ANGLE,
            'k_dihedral': K_DIHEDRAL,
        }
    }

    return anchors

def save_anchors(anchors, sys_name):
    """Save anchor atoms to file."""
    out_path = BASE_DIR / sys_name / "boresch_anchors.json"
    out_path.parent.mkdir(parents=True, exist_ok=True)

    with open(out_path, 'w') as f:
        json.dump(anchors, f, indent=2)
    print(f"\n[SAVED] {out_path}")

    # Also save as numpy for compatibility
    np_path = BASE_DIR / sys_name / "boresch_anchors.npy"
    np.save(np_path, anchors)
    print(f"[SAVED] {np_path}")

def main():
    print("="*70)
    print("BORESCH RESTRAINT ANCHOR SELECTION")
    print("Natural Product CID_10120")
    print("="*70)
    print()

    # Parse input PDB
    print(f"Loading: {INPUT_PDB}")
    atoms = parse_pdb_simple(INPUT_PDB)
    print(f"Total atoms: {len(atoms)}")

    # Select anchors
    anchors = select_boresch_anchors(atoms)

    # Save for each system
    for sys_name in ['wt_complex', 'mut_complex']:
        save_anchors(anchors, sys_name)

    print("\n" + "="*70)
    print("SUMMARY")
    print("="*70)
    print(f"Protein anchors: P1={anchors['P1']['resid']}, P2={anchors['P2']['resid']}, P3={anchors['P3']['resid']}")
    print(f"Ligand anchors: L1={anchors['L1']['name']}, L2={anchors['L2']['name']}, L3={anchors['L3']['name']}")
    print(f"Distance r0: {anchors['equilibrium']['r0_nm']:.4f} nm")
    print()
    print("[DONE] Boresch anchors selected and saved")

if __name__ == "__main__":
    main()
