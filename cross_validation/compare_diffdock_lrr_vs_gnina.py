#!/usr/bin/env python3
"""
Cross-validation: Compare DiffDock-LRR vs GNINA docking poses

DiffDock used full NOD2 structure with original numbering (1-1040)
GNINA also used same structure

Key binding residues:
  - GLU1008, ASN1010, ASP1011
  - ARG1034, ARG1037
"""

import numpy as np
import os
import re

# Paths
DIFFDOCK_LRR_DIR = r"C:\Users\vasud\nod2-screening-data\cross_validation\diffdock_lrr"
GNINA_SDF = r"C:\Users\vasud\nod2-screening-data\PHASE_6\structures\febuxostat_docked.sdf"
OUTPUT_DIR = r"C:\Users\vasud\nod2-screening-data\cross_validation\comparison_results"

# Key binding residues (original numbering)
KEY_RESIDUES = {
    'GLU1008': 1008,
    'ASN1010': 1010,
    'ASP1011': 1011,
    'ARG1034': 1034,
    'ARG1037': 1037
}

# Pocket residues for contact analysis
POCKET_RESIDUES = [1007, 1008, 1009, 1010, 1011, 1034, 1035, 1036, 1037]


def load_sdf_coords(sdf_file, pose_idx=0):
    """Load ligand heavy atom coordinates from SDF file (specific pose)."""
    coords = []
    with open(sdf_file, 'r') as f:
        content = f.read()

    # Split by molecule separator
    molecules = content.split('$$$$')
    if pose_idx >= len(molecules):
        return None

    lines = molecules[pose_idx].strip().split('\n')
    if len(lines) < 4:
        return None

    # Find counts line (contains V2000 or V3000)
    counts_idx = None
    for i, line in enumerate(lines):
        if 'V2000' in line or 'V3000' in line:
            counts_idx = i
            break

    if counts_idx is None:
        return None

    counts_line = lines[counts_idx]
    try:
        n_atoms = int(counts_line[:3].strip())
    except ValueError:
        return None

    # Atoms start right after counts line
    for i in range(counts_idx + 1, counts_idx + 1 + n_atoms):
        if i >= len(lines):
            break
        line = lines[i]
        try:
            x = float(line[0:10].strip())
            y = float(line[10:20].strip())
            z = float(line[20:30].strip())
            atom_symbol = line[31:34].strip()
            if atom_symbol != 'H':  # Skip hydrogens
                coords.append([x, y, z])
        except (ValueError, IndexError):
            continue

    return np.array(coords) if coords else None


def parse_pdb_resnum(line):
    """Parse residue number from PDB ATOM line, handling 4-digit numbers."""
    # Columns 23-26 for resSeq, but for 4-digit can overflow into column 22
    resnum_str = line[22:27].strip()
    # Remove any chain identifier that might be there
    resnum_str = re.sub(r'[A-Z]', '', resnum_str)
    try:
        return int(resnum_str)
    except ValueError:
        return None


def load_pdb_residue_coords(pdb_file, residue_nums, atom_name='CA'):
    """Load CA coordinates for specific residues from PDB."""
    coords = {}
    with open(pdb_file, 'r') as f:
        for line in f:
            if line.startswith('ATOM'):
                atom = line[12:16].strip()
                if atom != atom_name:
                    continue
                resnum = parse_pdb_resnum(line)
                if resnum and resnum in residue_nums:
                    try:
                        x = float(line[30:38].strip())
                        y = float(line[38:46].strip())
                        z = float(line[46:54].strip())
                        resname = line[17:20].strip()
                        coords[resnum] = {'xyz': np.array([x, y, z]), 'name': resname}
                    except (ValueError, IndexError):
                        continue
    return coords


def load_all_pocket_atoms(pdb_file, residue_nums):
    """Load all heavy atom coordinates for pocket residues."""
    atoms = []
    with open(pdb_file, 'r') as f:
        for line in f:
            if line.startswith('ATOM'):
                resnum = parse_pdb_resnum(line)
                if resnum and resnum in residue_nums:
                    atom_name = line[12:16].strip()
                    if atom_name.startswith('H') or len(atom_name) > 1 and atom_name[1] == 'H':
                        continue
                    try:
                        x = float(line[30:38].strip())
                        y = float(line[38:46].strip())
                        z = float(line[46:54].strip())
                        resname = line[17:20].strip()
                        atoms.append({
                            'xyz': np.array([x, y, z]),
                            'resnum': resnum,
                            'resname': resname,
                            'atom': atom_name
                        })
                    except (ValueError, IndexError):
                        continue
    return atoms


def calculate_centroid(coords):
    """Calculate centroid of coordinates."""
    return np.mean(coords, axis=0)


def calculate_rmsd_no_alignment(coords1, coords2):
    """Calculate RMSD between two coordinate sets (centroids only - different sizes ok)."""
    c1 = calculate_centroid(coords1)
    c2 = calculate_centroid(coords2)
    return np.linalg.norm(c1 - c2)


def find_contacts(ligand_coords, pocket_atoms, cutoff=4.0):
    """Find residues in contact with ligand (any heavy atom within cutoff)."""
    contacts = {}
    for lig_coord in ligand_coords:
        for atom in pocket_atoms:
            dist = np.linalg.norm(lig_coord - atom['xyz'])
            if dist < cutoff:
                key = f"{atom['resname']}{atom['resnum']}"
                if key not in contacts or dist < contacts[key]['min_dist']:
                    contacts[key] = {
                        'resnum': atom['resnum'],
                        'resname': atom['resname'],
                        'min_dist': dist,
                        'atom': atom['atom']
                    }
    return contacts


def find_hbonds(ligand_coords, pocket_atoms, cutoff=3.5):
    """Simple H-bond detection based on N/O distance."""
    hbond_atoms = ['N', 'O', 'OE1', 'OE2', 'OD1', 'OD2', 'ND2', 'NH1', 'NH2', 'NE']
    hbonds = []

    for lig_coord in ligand_coords:
        for atom in pocket_atoms:
            if atom['atom'] in hbond_atoms:
                dist = np.linalg.norm(lig_coord - atom['xyz'])
                if dist < cutoff:
                    hbonds.append({
                        'residue': f"{atom['resname']}{atom['resnum']}",
                        'atom': atom['atom'],
                        'distance': dist
                    })

    return hbonds


def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    print("=" * 70)
    print("CROSS-VALIDATION: DiffDock-LRR vs GNINA Docking Comparison")
    print("=" * 70)

    # Load protein
    diffdock_protein = os.path.join(DIFFDOCK_LRR_DIR, "protein_no_ligands.pdb")
    diffdock_rank1 = os.path.join(DIFFDOCK_LRR_DIR, "rank1_confidence-0.30.sdf")

    print("\n1. LOADING STRUCTURES")
    print("-" * 40)

    # Load DiffDock rank1 pose
    dd_coords = load_sdf_coords(diffdock_rank1)
    print(f"DiffDock rank1: {len(dd_coords)} heavy atoms")
    dd_centroid = calculate_centroid(dd_coords)
    print(f"  Centroid: ({dd_centroid[0]:.2f}, {dd_centroid[1]:.2f}, {dd_centroid[2]:.2f})")

    # Load GNINA pose (first molecule = best)
    gnina_coords = load_sdf_coords(GNINA_SDF, pose_idx=0)
    print(f"GNINA pose 1: {len(gnina_coords)} heavy atoms")
    gnina_centroid = calculate_centroid(gnina_coords)
    print(f"  Centroid: ({gnina_centroid[0]:.2f}, {gnina_centroid[1]:.2f}, {gnina_centroid[2]:.2f})")

    # Distance between ligand centroids
    centroid_dist = np.linalg.norm(dd_centroid - gnina_centroid)
    print(f"\n  Distance between centroids: {centroid_dist:.2f} A")

    # Load pocket residues
    print(f"\n2. POCKET ANALYSIS")
    print("-" * 40)

    pocket_ca = load_pdb_residue_coords(diffdock_protein, POCKET_RESIDUES)
    print(f"Found {len(pocket_ca)} pocket residues with CA atoms:")
    for resnum in sorted(pocket_ca.keys()):
        res = pocket_ca[resnum]
        print(f"  {res['name']}{resnum}: ({res['xyz'][0]:.2f}, {res['xyz'][1]:.2f}, {res['xyz'][2]:.2f})")

    # Calculate pocket center
    if pocket_ca:
        pocket_center = np.mean([res['xyz'] for res in pocket_ca.values()], axis=0)
        print(f"\nPocket center: ({pocket_center[0]:.2f}, {pocket_center[1]:.2f}, {pocket_center[2]:.2f})")

        # Distance from pocket center
        dd_to_pocket = np.linalg.norm(dd_centroid - pocket_center)
        gnina_to_pocket = np.linalg.norm(gnina_centroid - pocket_center)
        print(f"\nDistance to pocket center:")
        print(f"  DiffDock rank1: {dd_to_pocket:.2f} A")
        print(f"  GNINA pose 1:   {gnina_to_pocket:.2f} A")
    else:
        print("ERROR: No pocket residues found!")
        return

    # Load all pocket atoms
    pocket_atoms = load_all_pocket_atoms(diffdock_protein, POCKET_RESIDUES)
    print(f"\nLoaded {len(pocket_atoms)} heavy atoms from pocket residues")

    print(f"\n3. CONTACT ANALYSIS")
    print("-" * 40)

    # DiffDock contacts
    dd_contacts = find_contacts(dd_coords, pocket_atoms, cutoff=4.0)
    print(f"\nDiffDock rank1 contacts (< 4.0 A): {len(dd_contacts)}")
    for key, info in sorted(dd_contacts.items(), key=lambda x: x[1]['min_dist']):
        print(f"  {info['resname']}{info['resnum']}: {info['min_dist']:.2f} A ({info['atom']})")

    # GNINA contacts
    gnina_contacts = find_contacts(gnina_coords, pocket_atoms, cutoff=4.0)
    print(f"\nGNINA pose 1 contacts (< 4.0 A): {len(gnina_contacts)}")
    for key, info in sorted(gnina_contacts.items(), key=lambda x: x[1]['min_dist']):
        print(f"  {info['resname']}{info['resnum']}: {info['min_dist']:.2f} A ({info['atom']})")

    print(f"\n4. KEY RESIDUE CONTACT CHECK")
    print("-" * 40)

    print("\nDiffDock rank1:")
    dd_key_count = 0
    for name, resnum in KEY_RESIDUES.items():
        found = False
        for key, info in dd_contacts.items():
            if info['resnum'] == resnum:
                print(f"  {name}: YES ({info['min_dist']:.2f} A)")
                found = True
                dd_key_count += 1
                break
        if not found:
            print(f"  {name}: NO")

    print("\nGNINA pose 1:")
    gnina_key_count = 0
    for name, resnum in KEY_RESIDUES.items():
        found = False
        for key, info in gnina_contacts.items():
            if info['resnum'] == resnum:
                print(f"  {name}: YES ({info['min_dist']:.2f} A)")
                found = True
                gnina_key_count += 1
                break
        if not found:
            print(f"  {name}: NO")

    # H-bond analysis
    print(f"\n5. H-BOND ANALYSIS")
    print("-" * 40)

    dd_hbonds = find_hbonds(dd_coords, pocket_atoms)
    print(f"\nDiffDock rank1 potential H-bonds: {len(dd_hbonds)}")
    seen = set()
    for hb in sorted(dd_hbonds, key=lambda x: x['distance']):
        key = (hb['residue'], hb['atom'])
        if key not in seen:
            seen.add(key)
            print(f"  {hb['residue']} {hb['atom']}: {hb['distance']:.2f} A")

    gnina_hbonds = find_hbonds(gnina_coords, pocket_atoms)
    print(f"\nGNINA pose 1 potential H-bonds: {len(gnina_hbonds)}")
    seen = set()
    for hb in sorted(gnina_hbonds, key=lambda x: x['distance']):
        key = (hb['residue'], hb['atom'])
        if key not in seen:
            seen.add(key)
            print(f"  {hb['residue']} {hb['atom']}: {hb['distance']:.2f} A")

    # Analyze top 10 DiffDock poses
    print(f"\n6. TOP 10 DIFFDOCK POSES ANALYSIS")
    print("-" * 40)

    best_pose = None
    best_dist = float('inf')

    for rank in range(1, 11):
        sdf_files = [f for f in os.listdir(DIFFDOCK_LRR_DIR)
                     if f.startswith(f'rank{rank}_') and f.endswith('.sdf')]
        if sdf_files:
            sdf_path = os.path.join(DIFFDOCK_LRR_DIR, sdf_files[0])
            coords = load_sdf_coords(sdf_path)
            if coords is not None:
                centroid = calculate_centroid(coords)
                dist_to_pocket = np.linalg.norm(centroid - pocket_center)
                contacts = find_contacts(coords, pocket_atoms, cutoff=4.0)

                conf = sdf_files[0].split('confidence')[1].replace('.sdf', '')

                key_hits = []
                for name, resnum in KEY_RESIDUES.items():
                    for c in contacts.values():
                        if c['resnum'] == resnum:
                            key_hits.append(name.replace('1008', '').replace('1010', '').replace('1011', '').replace('1034', '').replace('1037', ''))
                            break

                print(f"\n  Rank {rank} (conf={conf}):")
                print(f"    Distance to pocket: {dist_to_pocket:.2f} A")
                print(f"    Total contacts: {len(contacts)}")
                print(f"    Key residues: {len(key_hits)}/5 ({', '.join(key_hits) if key_hits else 'None'})")

                if dist_to_pocket < best_dist:
                    best_dist = dist_to_pocket
                    best_pose = rank

    # Summary
    print(f"\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)

    agreement = "YES" if centroid_dist < 5.0 else "NO" if centroid_dist > 20.0 else "PARTIAL"

    print(f"""
DiffDock-LRR (Rank 1):
  - Distance to pocket center: {dd_to_pocket:.2f} A
  - Key residues contacted: {dd_key_count}/5
  - Potential H-bonds: {len(set((h['residue'], h['atom']) for h in dd_hbonds))}

GNINA (Best Pose):
  - Distance to pocket center: {gnina_to_pocket:.2f} A
  - Key residues contacted: {gnina_key_count}/5
  - Potential H-bonds: {len(set((h['residue'], h['atom']) for h in gnina_hbonds))}

Comparison:
  - Distance between ligand centroids: {centroid_dist:.2f} A
  - Binding site agreement: {agreement}
  - Best DiffDock pose by proximity: Rank {best_pose} ({best_dist:.2f} A)
""")

    if centroid_dist < 5.0:
        print("  CONCLUSION: DiffDock and GNINA agree on binding site!")
    elif dd_to_pocket < 10.0 and gnina_to_pocket < 10.0:
        print("  CONCLUSION: Both methods place ligand in same pocket region")
    else:
        print(f"  CONCLUSION: Methods disagree - DiffDock found different site")

    # Save results
    results_file = os.path.join(OUTPUT_DIR, "diffdock_lrr_vs_gnina.txt")
    with open(results_file, 'w') as f:
        f.write("CROSS-VALIDATION: DiffDock-LRR vs GNINA\n")
        f.write("=" * 50 + "\n\n")
        f.write(f"DiffDock rank1 distance to pocket: {dd_to_pocket:.2f} A\n")
        f.write(f"GNINA pose1 distance to pocket: {gnina_to_pocket:.2f} A\n")
        f.write(f"Distance between centroids: {centroid_dist:.2f} A\n\n")
        f.write(f"DiffDock key residues: {dd_key_count}/5\n")
        f.write(f"GNINA key residues: {gnina_key_count}/5\n\n")
        f.write(f"Agreement: {agreement}\n")

    print(f"\nResults saved to: {results_file}")


if __name__ == "__main__":
    main()
