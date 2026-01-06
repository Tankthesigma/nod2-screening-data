#!/usr/bin/env python3
"""
PHASE A1: Verify Ligand Poses for Mutant MD Simulations

Verifies that ligand coordinates from Phase 6 are correctly transferred
and positioned near the binding pocket.
"""

import os
import sys
import math

# Paths
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
LIGANDS_DIR = os.path.join(SCRIPT_DIR, "ligands")
STRUCTURES_DIR = os.path.join(SCRIPT_DIR, "structures")

# Ligand files
FEBUXOSTAT_PDB = os.path.join(LIGANDS_DIR, "ligand_febuxostat.pdb")
NATURAL_PDB = os.path.join(LIGANDS_DIR, "ligand_natural_cid10592.pdb")

# Reference protein (for pocket distance check)
WT_PDB = os.path.join(os.path.dirname(SCRIPT_DIR), "PHASE_6", "structures", "NOD2_LRR_clean.pdb")

# Binding pocket residues
POCKET_RESIDUES = [1008, 1011, 1037]


def parse_pdb_atoms(pdb_path):
    """Extract atom coordinates from PDB file."""
    atoms = []
    with open(pdb_path, 'r') as f:
        for line in f:
            if line.startswith('ATOM') or line.startswith('HETATM'):
                try:
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    atom_name = line[12:16].strip()
                    res_name = line[17:20].strip()
                    res_num = int(line[22:26].strip()) if line[22:26].strip().isdigit() else 0
                    atoms.append({
                        'name': atom_name,
                        'res_name': res_name,
                        'res_num': res_num,
                        'x': x, 'y': y, 'z': z
                    })
                except (ValueError, IndexError):
                    continue
    return atoms


def get_centroid(atoms):
    """Calculate centroid of atom coordinates."""
    if not atoms:
        return None
    x = sum(a['x'] for a in atoms) / len(atoms)
    y = sum(a['y'] for a in atoms) / len(atoms)
    z = sum(a['z'] for a in atoms) / len(atoms)
    return (x, y, z)


def distance(p1, p2):
    """Euclidean distance between two points."""
    return math.sqrt((p1[0]-p2[0])**2 + (p1[1]-p2[1])**2 + (p1[2]-p2[2])**2)


def get_pocket_centroid(protein_atoms):
    """Get centroid of binding pocket residues."""
    pocket_atoms = [a for a in protein_atoms if a['res_num'] in POCKET_RESIDUES]
    if not pocket_atoms:
        return None
    return get_centroid(pocket_atoms)


def calculate_rmsd(atoms1, atoms2):
    """Calculate RMSD between two sets of atoms (by atom name matching)."""
    if len(atoms1) != len(atoms2):
        return None

    # Sort by atom name for matching
    atoms1_sorted = sorted(atoms1, key=lambda a: a['name'])
    atoms2_sorted = sorted(atoms2, key=lambda a: a['name'])

    sum_sq = 0
    for a1, a2 in zip(atoms1_sorted, atoms2_sorted):
        sum_sq += (a1['x']-a2['x'])**2 + (a1['y']-a2['y'])**2 + (a1['z']-a2['z'])**2

    return math.sqrt(sum_sq / len(atoms1))


def verify_ligand(name, pdb_path, pocket_centroid):
    """Verify a ligand pose."""
    result = {
        'name': name,
        'file': pdb_path,
        'exists': os.path.exists(pdb_path),
        'n_atoms': 0,
        'centroid': None,
        'pocket_distance': None,
        'in_pocket': False,
        'status': 'FAIL'
    }

    if not result['exists']:
        return result

    atoms = parse_pdb_atoms(pdb_path)
    result['n_atoms'] = len(atoms)

    if atoms:
        centroid = get_centroid(atoms)
        result['centroid'] = centroid

        if pocket_centroid:
            dist = distance(centroid, pocket_centroid)
            result['pocket_distance'] = dist
            # Ligand should be within 15 Angstroms of pocket centroid
            result['in_pocket'] = dist < 15.0

    if result['exists'] and result['n_atoms'] > 0 and result['in_pocket']:
        result['status'] = 'PASS'

    return result


def main():
    print("=" * 60)
    print("PHASE A1: Ligand Pose Verification")
    print("=" * 60)

    # Load protein and get pocket centroid
    print("\n1. Loading protein structure...")
    if not os.path.exists(WT_PDB):
        print(f"  ERROR: Protein PDB not found: {WT_PDB}")
        sys.exit(1)

    protein_atoms = parse_pdb_atoms(WT_PDB)
    pocket_centroid = get_pocket_centroid(protein_atoms)

    if pocket_centroid:
        print(f"  Pocket centroid (residues {POCKET_RESIDUES}):")
        print(f"    ({pocket_centroid[0]:.2f}, {pocket_centroid[1]:.2f}, {pocket_centroid[2]:.2f})")
    else:
        print("  WARNING: Could not calculate pocket centroid")

    # Verify ligands
    print("\n2. Verifying ligand poses...")
    print("-" * 60)

    ligands = [
        ("Febuxostat", FEBUXOSTAT_PDB),
        ("20a-Dihydrocortisol", NATURAL_PDB)
    ]

    results = []
    all_pass = True

    for name, pdb_path in ligands:
        result = verify_ligand(name, pdb_path, pocket_centroid)
        results.append(result)

        print(f"\n  {name}:")
        print(f"    File: {os.path.basename(pdb_path)}")
        print(f"    Atoms: {result['n_atoms']}")

        if result['centroid']:
            print(f"    Centroid: ({result['centroid'][0]:.2f}, {result['centroid'][1]:.2f}, {result['centroid'][2]:.2f})")

        if result['pocket_distance'] is not None:
            print(f"    Distance to pocket: {result['pocket_distance']:.2f} A")
            print(f"    In pocket (<15 A): {'Yes' if result['in_pocket'] else 'No'}")

        status_str = "[OK]" if result['status'] == 'PASS' else "[FAIL]"
        print(f"    Status: {status_str}")

        if result['status'] != 'PASS':
            all_pass = False

    # Write verification log
    print("\n" + "=" * 60)
    print("3. Writing verification log...")
    print("=" * 60)

    log_path = os.path.join(LIGANDS_DIR, "ligand_verification_log.txt")
    with open(log_path, 'w') as f:
        f.write("=" * 60 + "\n")
        f.write("PHASE A1: Ligand Pose Verification Log\n")
        f.write("=" * 60 + "\n\n")

        f.write("BINDING POCKET REFERENCE\n")
        f.write("-" * 40 + "\n")
        f.write(f"Pocket residues: {POCKET_RESIDUES}\n")
        if pocket_centroid:
            f.write(f"Pocket centroid: ({pocket_centroid[0]:.2f}, {pocket_centroid[1]:.2f}, {pocket_centroid[2]:.2f})\n")
        f.write("\n")

        f.write("LIGAND VERIFICATION\n")
        f.write("-" * 40 + "\n\n")

        for r in results:
            f.write(f"{r['name']}:\n")
            f.write(f"  File: {os.path.basename(r['file'])}\n")
            f.write(f"  Atoms: {r['n_atoms']}\n")
            if r['centroid']:
                f.write(f"  Centroid: ({r['centroid'][0]:.2f}, {r['centroid'][1]:.2f}, {r['centroid'][2]:.2f})\n")
            if r['pocket_distance'] is not None:
                f.write(f"  Distance to pocket: {r['pocket_distance']:.2f} A\n")
            f.write(f"  Status: {r['status']}\n\n")

        f.write("=" * 60 + "\n")
        f.write(f"OVERALL: {'PASSED' if all_pass else 'FAILED'}\n")
        f.write("=" * 60 + "\n")

    print(f"\nLog saved: {log_path}")

    # Summary
    print("\n" + "=" * 60)
    print("CHECKPOINT 2 COMPLETE")
    print("=" * 60)

    print(f"\nLigand poses transferred from Phase 6:")
    for r in results:
        status_str = "[OK]" if r['status'] == 'PASS' else "[FAIL]"
        print(f"  {status_str} {r['name']}: {r['pocket_distance']:.1f} A from pocket")

    if all_pass:
        print("\n[OK] All ligand poses verified - ready for MD setup")
    else:
        print("\n[FAIL] Some ligand poses failed verification!")
        sys.exit(1)


if __name__ == "__main__":
    main()
