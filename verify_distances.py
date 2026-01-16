#!/usr/bin/env python
"""
Check A: Verify 80 Angstrom distance and protein dimensions
"""
import numpy as np
import os

BASE = "C:/Users/vasud/nod2-screening-data/fep_complete/fep_pmx"

def parse_pdb_coords(pdb_file):
    """Parse PDB file and extract atom information."""
    atoms = []
    with open(pdb_file, 'r') as f:
        for line in f:
            if line.startswith('ATOM') or line.startswith('HETATM'):
                try:
                    atom_num = int(line[6:11].strip())
                    atom_name = line[12:16].strip()
                    res_name = line[17:20].strip()
                    chain = line[21]
                    res_num = int(line[22:26].strip())
                    x = float(line[30:38].strip())
                    y = float(line[38:46].strip())
                    z = float(line[46:54].strip())
                    atoms.append({
                        'atom_num': atom_num,
                        'atom_name': atom_name,
                        'res_name': res_name,
                        'chain': chain,
                        'res_num': res_num,
                        'x': x, 'y': y, 'z': z
                    })
                except (ValueError, IndexError):
                    continue
    return atoms

def get_protein_atoms(atoms):
    """Filter to get only protein atoms (chain A, not water/ions)."""
    protein_res = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY',
                   'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER',
                   'THR', 'TRP', 'TYR', 'VAL', 'HIE', 'HID', 'HIP']
    return [a for a in atoms if a['res_name'] in protein_res and a['chain'] == 'A']

def get_ligand_atoms(atoms):
    """Get ligand atoms (UNK residue)."""
    return [a for a in atoms if a['res_name'] == 'UNK']

def get_residue_atoms(atoms, res_num, chain='A'):
    """Get atoms for a specific residue."""
    return [a for a in atoms if a['res_num'] == res_num and a['chain'] == chain]

def calc_distance(a1, a2):
    """Calculate distance between two atoms."""
    return np.sqrt((a1['x']-a2['x'])**2 + (a1['y']-a2['y'])**2 + (a1['z']-a2['z'])**2)

def main():
    print("=" * 70)
    print("CHECK A: 80 ANGSTROM DISTANCE VERIFICATION")
    print("=" * 70)

    for sys_name in ['wt_complex', 'mut_complex']:
        print(f"\n{'='*70}")
        print(f"SYSTEM: {sys_name}")
        print(f"{'='*70}")

        pdb_file = f"{BASE}/{sys_name}/topology.pdb"
        if not os.path.exists(pdb_file):
            print(f"  ERROR: {pdb_file} not found")
            continue

        atoms = parse_pdb_coords(pdb_file)
        protein_atoms = get_protein_atoms(atoms)
        ligand_atoms = get_ligand_atoms(atoms)

        print(f"\n1. ATOM COUNTS:")
        print(f"   Total atoms: {len(atoms)}")
        print(f"   Protein atoms: {len(protein_atoms)}")
        print(f"   Ligand atoms: {len(ligand_atoms)}")

        # Protein dimensions
        if protein_atoms:
            xs = [a['x'] for a in protein_atoms]
            ys = [a['y'] for a in protein_atoms]
            zs = [a['z'] for a in protein_atoms]

            print(f"\n2. PROTEIN DIMENSIONS:")
            print(f"   X range: {min(xs):.1f} to {max(xs):.1f} = {max(xs)-min(xs):.1f} A")
            print(f"   Y range: {min(ys):.1f} to {max(ys):.1f} = {max(ys)-min(ys):.1f} A")
            print(f"   Z range: {min(zs):.1f} to {max(zs):.1f} = {max(zs)-min(zs):.1f} A")
            print(f"   Protein extent: {max(xs)-min(xs):.1f} x {max(ys)-min(ys):.1f} x {max(zs)-min(zs):.1f} A")

        # Residue 702 location
        res702 = get_residue_atoms(atoms, 702, 'A')
        res702_ca = [a for a in res702 if a['atom_name'] == 'CA']

        if res702:
            res_name = res702[0]['res_name']
            print(f"\n3. RESIDUE 702:")
            print(f"   Residue type: {res_name}")
            print(f"   Number of atoms: {len(res702)}")
            if res702_ca:
                ca = res702_ca[0]
                print(f"   CA position: ({ca['x']:.2f}, {ca['y']:.2f}, {ca['z']:.2f})")
        else:
            print(f"\n3. RESIDUE 702: NOT FOUND!")
            continue

        # Ligand center and distances
        if ligand_atoms:
            lig_xs = [a['x'] for a in ligand_atoms]
            lig_ys = [a['y'] for a in ligand_atoms]
            lig_zs = [a['z'] for a in ligand_atoms]
            lig_com = (np.mean(lig_xs), np.mean(lig_ys), np.mean(lig_zs))

            print(f"\n4. LIGAND:")
            print(f"   Number of atoms: {len(ligand_atoms)}")
            print(f"   Center of mass: ({lig_com[0]:.2f}, {lig_com[1]:.2f}, {lig_com[2]:.2f})")

            # Distance from R702 CA to ligand COM
            if res702_ca:
                ca = res702_ca[0]
                dist_com = np.sqrt((ca['x']-lig_com[0])**2 + (ca['y']-lig_com[1])**2 + (ca['z']-lig_com[2])**2)
                print(f"\n5. DISTANCES:")
                print(f"   R702 CA to ligand COM: {dist_com:.2f} A")

                # Distance to closest ligand atom
                min_dist = float('inf')
                closest_atom = None
                for lig_atom in ligand_atoms:
                    d = calc_distance(ca, lig_atom)
                    if d < min_dist:
                        min_dist = d
                        closest_atom = lig_atom

                print(f"   R702 CA to closest ligand atom: {min_dist:.2f} A ({closest_atom['atom_name']})")

                # Distance from R702 guanidinium (NH1, NH2) to ligand
                nh_atoms = [a for a in res702 if a['atom_name'] in ['NH1', 'NH2', 'NE1', 'NE']]
                if nh_atoms:
                    for nh in nh_atoms:
                        min_nh_dist = min(calc_distance(nh, la) for la in ligand_atoms)
                        print(f"   R702 {nh['atom_name']} to closest ligand: {min_nh_dist:.2f} A")

        # Binding site residues
        print(f"\n6. LRR BINDING SITE (residues 1007-1011):")
        for res_num in [1007, 1008, 1009, 1010, 1011]:
            res_atoms = get_residue_atoms(atoms, res_num, 'A')
            if res_atoms:
                ca = [a for a in res_atoms if a['atom_name'] == 'CA']
                if ca:
                    ca = ca[0]
                    # Distance from binding site to ligand
                    if ligand_atoms:
                        min_dist = min(calc_distance(ca, la) for la in ligand_atoms)
                        print(f"   {res_atoms[0]['res_name']}{res_num} CA to ligand: {min_dist:.2f} A")

    # Box dimensions check
    print(f"\n{'='*70}")
    print("BOX DIMENSION CHECK")
    print(f"{'='*70}")
    pdb_file = f"{BASE}/wt_complex/topology.pdb"
    with open(pdb_file, 'r') as f:
        for line in f:
            if line.startswith('CRYST1'):
                a = float(line[6:15])
                b = float(line[15:24])
                c = float(line[24:33])
                print(f"Box dimensions: {a:.2f} x {b:.2f} x {c:.2f} A")
                print(f"Can contain 80 A separation? {'YES' if min(a,b,c) > 80 else 'MARGINAL'}")
                break

if __name__ == "__main__":
    main()
