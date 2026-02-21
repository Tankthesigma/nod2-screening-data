#!/usr/bin/env python
"""
Check A: Verify 80 Angstrom distance and protein dimensions
No numpy required - pure Python
"""
import os
import math

BASE = "C:/Users/vasud/nod2-screening-data/fep_complete/fep_pmx"

def parse_pdb_coords(pdb_file):
    """Parse PDB file and extract atom information."""
    atoms = []
    with open(pdb_file, 'r') as f:
        for line in f:
            if line.startswith('ATOM') or line.startswith('HETATM'):
                try:
                    atom_name = line[12:16].strip()
                    res_name = line[17:20].strip()
                    chain = line[21]
                    res_num = int(line[22:26].strip())
                    x = float(line[30:38].strip())
                    y = float(line[38:46].strip())
                    z = float(line[46:54].strip())
                    atoms.append({
                        'atom_name': atom_name,
                        'res_name': res_name,
                        'chain': chain,
                        'res_num': res_num,
                        'x': x, 'y': y, 'z': z
                    })
                except (ValueError, IndexError):
                    continue
    return atoms

def calc_distance(a1, a2):
    """Calculate distance between two atoms."""
    return math.sqrt((a1['x']-a2['x'])**2 + (a1['y']-a2['y'])**2 + (a1['z']-a2['z'])**2)

def mean(values):
    return sum(values) / len(values)

def main():
    print("=" * 70)
    print("CHECK A: 80 ANGSTROM DISTANCE VERIFICATION")
    print("=" * 70)

    protein_res = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY',
                   'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER',
                   'THR', 'TRP', 'TYR', 'VAL', 'HIE', 'HID', 'HIP']

    for sys_name in ['wt_complex', 'mut_complex']:
        print(f"\n{'='*70}")
        print(f"SYSTEM: {sys_name}")
        print(f"{'='*70}")

        pdb_file = f"{BASE}/{sys_name}/topology.pdb"
        if not os.path.exists(pdb_file):
            print(f"  ERROR: {pdb_file} not found")
            continue

        atoms = parse_pdb_coords(pdb_file)
        protein_atoms = [a for a in atoms if a['res_name'] in protein_res and a['chain'] == 'A']
        ligand_atoms = [a for a in atoms if a['res_name'] == 'UNK']

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
            print(f"   X: {min(xs):.1f} to {max(xs):.1f} = {max(xs)-min(xs):.1f} A")
            print(f"   Y: {min(ys):.1f} to {max(ys):.1f} = {max(ys)-min(ys):.1f} A")
            print(f"   Z: {min(zs):.1f} to {max(zs):.1f} = {max(zs)-min(zs):.1f} A")

        # Residue 702 location
        res702 = [a for a in atoms if a['res_num'] == 702 and a['chain'] == 'A']
        res702_ca = [a for a in res702 if a['atom_name'] == 'CA']

        if res702:
            res_name = res702[0]['res_name']
            print(f"\n3. RESIDUE 702:")
            print(f"   Type: {res_name}")
            print(f"   Atoms: {len(res702)}")
            if res702_ca:
                ca = res702_ca[0]
                print(f"   CA: ({ca['x']:.2f}, {ca['y']:.2f}, {ca['z']:.2f})")

        # Ligand distances
        if ligand_atoms and res702_ca:
            lig_com = (mean([a['x'] for a in ligand_atoms]),
                       mean([a['y'] for a in ligand_atoms]),
                       mean([a['z'] for a in ligand_atoms]))

            print(f"\n4. LIGAND:")
            print(f"   Atoms: {len(ligand_atoms)}")
            print(f"   COM: ({lig_com[0]:.2f}, {lig_com[1]:.2f}, {lig_com[2]:.2f})")

            ca = res702_ca[0]
            dist_com = math.sqrt((ca['x']-lig_com[0])**2 + (ca['y']-lig_com[1])**2 + (ca['z']-lig_com[2])**2)

            print(f"\n5. DISTANCES FROM RESIDUE 702:")
            print(f"   CA to ligand COM: {dist_com:.2f} A")

            # Closest ligand atom
            min_dist = min(calc_distance(ca, la) for la in ligand_atoms)
            print(f"   CA to closest ligand atom: {min_dist:.2f} A")

        # Binding site distances
        print(f"\n6. LRR BINDING SITE DISTANCES:")
        for res_num in [1007, 1008, 1009, 1010, 1011]:
            res = [a for a in atoms if a['res_num'] == res_num and a['chain'] == 'A']
            if res and ligand_atoms:
                ca = [a for a in res if a['atom_name'] == 'CA']
                if ca:
                    min_dist = min(calc_distance(ca[0], la) for la in ligand_atoms)
                    print(f"   {res[0]['res_name']}{res_num} CA to ligand: {min_dist:.2f} A")

    # Box dimensions
    print(f"\n{'='*70}")
    print("BOX DIMENSIONS")
    print(f"{'='*70}")
    with open(f"{BASE}/wt_complex/topology.pdb", 'r') as f:
        for line in f:
            if line.startswith('CRYST1'):
                a = float(line[6:15])
                b = float(line[15:24])
                c = float(line[24:33])
                print(f"Box: {a:.2f} x {b:.2f} x {c:.2f} A")
                print(f"80 A separation possible? {'YES' if min(a,b,c) > 80 else 'MARGINAL'}")
                break

if __name__ == "__main__":
    main()
