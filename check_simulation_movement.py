#!/usr/bin/env python
"""
Check B: Verify protein movement during FEP simulations
Compare initial positions to final positions at various windows
"""
import os
import math
import struct

BASE = "C:/Users/vasud/nod2-screening-data/fep_complete/fep_pmx"

def load_npy_positions(filepath):
    """Load positions from .npy file without numpy."""
    with open(filepath, 'rb') as f:
        # Read magic string
        magic = f.read(6)
        if magic[:5] != b'\x93NUMPY':
            raise ValueError("Not a valid numpy file")

        version = (magic[5], f.read(1)[0]) if len(magic) > 5 else (1, 0)

        # Read header length
        if version[0] == 1:
            header_len = struct.unpack('<H', f.read(2))[0]
        else:
            header_len = struct.unpack('<I', f.read(4))[0]

        # Read header
        header = f.read(header_len).decode('utf-8')

        # Parse shape from header
        import re
        shape_match = re.search(r"'shape':\s*\(([^)]+)\)", header)
        if shape_match:
            shape = tuple(int(x.strip()) for x in shape_match.group(1).split(',') if x.strip())
        else:
            raise ValueError("Could not parse shape")

        # Parse dtype
        dtype_match = re.search(r"'descr':\s*'([^']+)'", header)
        dtype_str = dtype_match.group(1) if dtype_match else '<f8'

        # Read data
        n_atoms = shape[0]
        coords = []
        for i in range(n_atoms):
            x = struct.unpack('<d', f.read(8))[0]
            y = struct.unpack('<d', f.read(8))[0]
            z = struct.unpack('<d', f.read(8))[0]
            coords.append((x, y, z))

        return coords

def get_ca_indices(pdb_file):
    """Get indices of CA atoms from PDB file."""
    ca_indices = []
    atom_idx = 0
    with open(pdb_file, 'r') as f:
        for line in f:
            if line.startswith('ATOM') or line.startswith('HETATM'):
                atom_name = line[12:16].strip()
                res_name = line[17:20].strip()
                chain = line[21]
                if atom_name == 'CA' and chain == 'A':
                    # Check if protein residue
                    protein_res = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY',
                                   'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER',
                                   'THR', 'TRP', 'TYR', 'VAL', 'HIE', 'HID', 'HIP']
                    if res_name in protein_res:
                        ca_indices.append(atom_idx)
                atom_idx += 1
    return ca_indices

def calculate_rmsd(coords1, coords2, indices):
    """Calculate RMSD for specified atom indices."""
    sum_sq = 0.0
    n = 0
    for idx in indices:
        if idx < len(coords1) and idx < len(coords2):
            dx = coords1[idx][0] - coords2[idx][0]
            dy = coords1[idx][1] - coords2[idx][1]
            dz = coords1[idx][2] - coords2[idx][2]
            sum_sq += dx*dx + dy*dy + dz*dz
            n += 1
    if n == 0:
        return 0.0
    return math.sqrt(sum_sq / n) * 10  # nm to Angstrom

def main():
    print("=" * 70)
    print("CHECK B: PROTEIN MOVEMENT VERIFICATION")
    print("=" * 70)
    print("\nComparing initial positions to final positions at each window")
    print("RMSD > 1 A indicates protein is NOT frozen\n")

    for sys_name in ['wt_complex', 'mut_complex', 'solvent']:
        print(f"\n{'='*70}")
        print(f"SYSTEM: {sys_name}")
        print(f"{'='*70}")

        base_path = f"{BASE}/{sys_name}"
        initial_pos_file = f"{base_path}/positions.npy"
        pdb_file = f"{base_path}/topology.pdb"

        if not os.path.exists(initial_pos_file):
            print(f"  Initial positions not found: {initial_pos_file}")
            continue

        if not os.path.exists(pdb_file):
            print(f"  Topology not found: {pdb_file}")
            continue

        # Get CA indices
        if sys_name != 'solvent':
            ca_indices = get_ca_indices(pdb_file)
            print(f"  CA atoms for RMSD: {len(ca_indices)}")
        else:
            # For solvent, use all atoms
            ca_indices = list(range(100))  # First 100 atoms
            print(f"  Using first 100 atoms for RMSD")

        # Load initial positions
        try:
            initial_coords = load_npy_positions(initial_pos_file)
            print(f"  Initial positions loaded: {len(initial_coords)} atoms")
        except Exception as e:
            print(f"  Error loading initial positions: {e}")
            continue

        # Check each window
        print(f"\n  Window | Final RMSD (A) | Status")
        print(f"  " + "-" * 40)

        for window_idx in [0, 5, 10, 15, 19]:
            final_pos_file = f"{base_path}/window_{window_idx:02d}/final_positions.npy"

            if os.path.exists(final_pos_file):
                try:
                    final_coords = load_npy_positions(final_pos_file)
                    rmsd = calculate_rmsd(initial_coords, final_coords, ca_indices)
                    status = "MOVED" if rmsd > 1.0 else "FROZEN?" if rmsd < 0.5 else "MINIMAL"
                    print(f"  {window_idx:6d} | {rmsd:14.2f} | {status}")
                except Exception as e:
                    print(f"  {window_idx:6d} | ERROR: {e}")
            else:
                print(f"  {window_idx:6d} | FILE NOT FOUND")

    print(f"\n{'='*70}")
    print("INTERPRETATION")
    print(f"{'='*70}")
    print("  RMSD > 1-2 A: Normal protein fluctuation (NOT frozen)")
    print("  RMSD > 3-5 A: Significant conformational change")
    print("  RMSD < 0.5 A: Potentially frozen or minimized only")

if __name__ == "__main__":
    main()
