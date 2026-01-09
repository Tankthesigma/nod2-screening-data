#!/usr/bin/env python3
"""
PHASE A1 ANALYSIS - CHECKPOINT 1: File Discovery
Finds all 12 trajectory/topology pairs and verifies them.
"""

import os
import sys
from pathlib import Path
import struct
import platform

# Detect if running in WSL or Windows
if 'microsoft' in platform.uname().release.lower() or os.path.exists('/mnt/c'):
    # WSL - use Unix-style path
    BASE_DIR = Path("/mnt/c/Users/vasud/nod2-screening-data/PHASE_A1_mutant_MD")
else:
    # Native Windows
    BASE_DIR = Path(r"C:\Users\vasud\nod2-screening-data\PHASE_A1_mutant_MD")

# All directories containing trajectory data
TRAJECTORY_DIRS = [
    BASE_DIR / "trajectories",
    BASE_DIR / "vast_downloads" / "5080",
    BASE_DIR / "vast_downloads" / "5090",
    BASE_DIR / "vast_downloads" / "new_5090",
]

def get_dcd_atom_count(dcd_path):
    """Read atom count from DCD file header."""
    try:
        with open(dcd_path, 'rb') as f:
            # Skip first block size (4 bytes)
            f.read(4)
            # Read header signature
            header = f.read(4)
            # Skip to NATOM position (after header block)
            # DCD format: 4 bytes block size, 4 bytes 'CORD', then header info
            f.seek(0)
            f.read(4)  # block size
            f.read(4)  # CORD
            # Read NSET (number of frames)
            nset = struct.unpack('i', f.read(4))[0]
            # Read ISTART
            f.read(4)
            # Read NSAVC
            f.read(4)
            # Skip 5 unused integers
            f.read(20)
            # Read NAMNF (number of fixed atoms, usually 0)
            f.read(4)
            # Read DELTA (timestep)
            f.read(4)
            # Skip 9 unused values
            f.read(36)
            # Read CHARMM version
            f.read(4)
            # End of first block
            f.read(4)  # block size end

            # Second block: title
            block_size = struct.unpack('i', f.read(4))[0]
            f.read(block_size)
            f.read(4)  # block size end

            # Third block: NATOM
            f.read(4)  # block size
            natom = struct.unpack('i', f.read(4))[0]

            return natom, nset
    except Exception as e:
        return None, str(e)

def get_pdb_atom_count(pdb_path):
    """Count ATOM/HETATM records in PDB file."""
    count = 0
    try:
        with open(pdb_path, 'r') as f:
            for line in f:
                if line.startswith('ATOM') or line.startswith('HETATM'):
                    count += 1
        return count
    except Exception as e:
        return None

def get_unique_resnames(pdb_path):
    """Get all unique residue names from PDB, excluding standard ones."""
    standard_residues = {
        # Amino acids
        'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE',
        'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL',
        # Common variants
        'HIE', 'HID', 'HIP', 'CYX', 'ASH', 'GLH',
        # Water
        'HOH', 'WAT', 'TIP3', 'TIP', 'SOL', 'H2O',
        # Ions
        'NA', 'CL', 'K', 'MG', 'CA', 'ZN', 'NA+', 'CL-', 'K+',
    }

    resnames = set()
    try:
        with open(pdb_path, 'r') as f:
            for line in f:
                if line.startswith('ATOM') or line.startswith('HETATM'):
                    resname = line[17:20].strip()
                    if resname and resname.upper() not in standard_residues:
                        resnames.add(resname)
        return resnames
    except Exception as e:
        return {f"ERROR: {e}"}

def parse_filename(filename):
    """Parse {Mutant}_{ligand}_rep{N} pattern."""
    # Remove extension
    base = filename.replace('.dcd', '').replace('_solvated.pdb', '')

    parts = base.split('_')
    if len(parts) >= 3:
        mutant = parts[0]  # R702W or G908R
        ligand = parts[1]  # febuxostat or natural
        rep_part = parts[2]  # rep1, rep2, rep3

        if rep_part.startswith('rep'):
            rep = int(rep_part.replace('rep', ''))
            return mutant, ligand, rep

    return None, None, None

def discover_files():
    """Find all DCD/PDB pairs across all directories."""
    discoveries = []

    for traj_dir in TRAJECTORY_DIRS:
        if not traj_dir.exists():
            print(f"WARNING: Directory not found: {traj_dir}")
            continue

        # Find all DCD files
        for dcd_file in traj_dir.glob("*.dcd"):
            mutant, ligand, rep = parse_filename(dcd_file.name)

            if mutant is None:
                print(f"WARNING: Could not parse filename: {dcd_file.name}")
                continue

            # Find matching solvated PDB
            pdb_name = dcd_file.stem + "_solvated.pdb"
            pdb_file = traj_dir / pdb_name

            if not pdb_file.exists():
                print(f"WARNING: Missing topology for {dcd_file.name}")
                continue

            # Get file sizes
            dcd_size_gb = dcd_file.stat().st_size / (1024**3)
            pdb_size_mb = pdb_file.stat().st_size / (1024**2)

            # Get atom counts
            dcd_atoms, dcd_frames = get_dcd_atom_count(dcd_file)
            pdb_atoms = get_pdb_atom_count(pdb_file)

            # Check match
            if dcd_atoms and pdb_atoms:
                atom_match = "YES" if dcd_atoms == pdb_atoms else f"NO ({dcd_atoms} vs {pdb_atoms})"
            else:
                atom_match = "ERROR"

            # Get ligand resnames
            ligand_resnames = get_unique_resnames(pdb_file)

            discoveries.append({
                'mutant': mutant,
                'ligand': ligand,
                'rep': rep,
                'dcd_path': str(dcd_file),
                'pdb_path': str(pdb_file),
                'dcd_size_gb': dcd_size_gb,
                'pdb_size_mb': pdb_size_mb,
                'dcd_atoms': dcd_atoms,
                'pdb_atoms': pdb_atoms,
                'atom_match': atom_match,
                'frames': dcd_frames,
                'ligand_resnames': ligand_resnames,
            })

    # Sort by mutant, ligand, rep
    discoveries.sort(key=lambda x: (x['mutant'], x['ligand'], x['rep']))

    return discoveries

def main():
    print("=" * 80)
    print("PHASE A1 ANALYSIS - CHECKPOINT 1: FILE DISCOVERY")
    print("=" * 80)
    print()

    # Create analysis directory
    analysis_dir = BASE_DIR / "analysis"
    analysis_dir.mkdir(exist_ok=True)

    # Discover files
    print("Scanning directories...")
    for d in TRAJECTORY_DIRS:
        print(f"  - {d}")
    print()

    discoveries = discover_files()

    print(f"Found {len(discoveries)} trajectory/topology pairs")
    print()

    # Print discovery table
    print("=" * 120)
    print(f"{'ID':<4} {'Mutant':<8} {'Ligand':<12} {'Rep':<4} {'DCD (GB)':<10} {'PDB (MB)':<10} {'Atoms Match':<15} {'Frames':<8} {'Ligand Resnames'}")
    print("=" * 120)

    all_ligand_resnames = {}

    for i, d in enumerate(discoveries, 1):
        resnames_str = ', '.join(d['ligand_resnames']) if d['ligand_resnames'] else 'NONE FOUND'
        print(f"{i:<4} {d['mutant']:<8} {d['ligand']:<12} {d['rep']:<4} {d['dcd_size_gb']:<10.2f} {d['pdb_size_mb']:<10.1f} {d['atom_match']:<15} {d['frames'] if d['frames'] else 'N/A':<8} {resnames_str}")

        # Collect all ligand resnames
        key = f"{d['mutant']}_{d['ligand']}"
        if key not in all_ligand_resnames:
            all_ligand_resnames[key] = set()
        all_ligand_resnames[key].update(d['ligand_resnames'])

    print("=" * 120)
    print()

    # Summary of ligand resnames
    print("LIGAND RESNAMES DETECTED:")
    print("-" * 40)
    for key, resnames in sorted(all_ligand_resnames.items()):
        print(f"  {key}: {', '.join(resnames) if resnames else 'NONE'}")
    print()

    # Verification summary
    print("VERIFICATION SUMMARY:")
    print("-" * 40)
    total = len(discoveries)
    matched = sum(1 for d in discoveries if d['atom_match'] == 'YES')
    print(f"  Total pairs found: {total}/12 expected")
    print(f"  Atom count matches: {matched}/{total}")

    if total != 12:
        print(f"  WARNING: Expected 12 simulations, found {total}!")
    if matched != total:
        print(f"  WARNING: Some atom counts don't match!")

    print()

    # Save discovery log
    log_path = analysis_dir / "file_discovery_log.txt"
    with open(log_path, 'w') as f:
        f.write("PHASE A1 FILE DISCOVERY LOG\n")
        f.write("=" * 80 + "\n\n")

        f.write(f"{'ID':<4} {'Mutant':<8} {'Ligand':<12} {'Rep':<4} {'DCD (GB)':<10} {'PDB (MB)':<10} {'Atoms Match':<15} {'Frames':<8}\n")
        f.write("-" * 80 + "\n")

        for i, d in enumerate(discoveries, 1):
            f.write(f"{i:<4} {d['mutant']:<8} {d['ligand']:<12} {d['rep']:<4} {d['dcd_size_gb']:<10.2f} {d['pdb_size_mb']:<10.1f} {d['atom_match']:<15} {d['frames'] if d['frames'] else 'N/A':<8}\n")

        f.write("\n\nLIGAND RESNAMES:\n")
        for key, resnames in sorted(all_ligand_resnames.items()):
            f.write(f"  {key}: {', '.join(resnames) if resnames else 'NONE'}\n")

        f.write("\n\nFULL PATHS:\n")
        for i, d in enumerate(discoveries, 1):
            f.write(f"\n{i}. {d['mutant']}_{d['ligand']}_rep{d['rep']}:\n")
            f.write(f"   DCD: {d['dcd_path']}\n")
            f.write(f"   PDB: {d['pdb_path']}\n")

    print(f"Discovery log saved to: {log_path}")

    # Save ligand resnames
    resnames_path = analysis_dir / "ligand_resnames_detected.txt"
    with open(resnames_path, 'w') as f:
        f.write("LIGAND RESNAMES DETECTED IN PHASE A1 TOPOLOGIES\n")
        f.write("=" * 50 + "\n\n")
        for key, resnames in sorted(all_ligand_resnames.items()):
            f.write(f"{key}: {', '.join(resnames) if resnames else 'NONE'}\n")

    print(f"Ligand resnames saved to: {resnames_path}")
    print()
    print("=" * 80)
    print("CHECKPOINT 1 COMPLETE - WAITING FOR APPROVAL")
    print("=" * 80)

    return discoveries

if __name__ == "__main__":
    discoveries = main()
