#!/usr/bin/env python3
"""
PHASE A1 ANALYSIS - Protein-Ligand Contact Analysis
Identifies residues in contact with ligand and calculates contact frequencies.
"""

import os
import sys
from pathlib import Path
import platform
import warnings
warnings.filterwarnings('ignore')

if 'microsoft' in platform.uname().release.lower() or os.path.exists('/mnt/c'):
    BASE_DIR = Path("/mnt/c/Users/vasud/nod2-screening-data/PHASE_A1_mutant_MD")
else:
    BASE_DIR = Path(r"C:\Users\vasud\nod2-screening-data\PHASE_A1_mutant_MD")

try:
    import MDAnalysis as mda
    from MDAnalysis.analysis import distances
    import numpy as np
    from collections import defaultdict
except ImportError:
    print("ERROR: MDAnalysis not installed!")
    sys.exit(1)

# ============================================================
# PARAMETERS (LOCKED)
# ============================================================
LIGAND_SELECTION = "resname UNK"
PROTEIN_SELECTION = "protein"
CONTACT_CUTOFF = 4.0  # Angstroms

# Key pocket residues to track
KEY_RESIDUES = [1008, 1011, 1037]  # GLU1008, ASP1011, ARG1037

TRAJECTORY_LOCATIONS = {
    'trajectories': BASE_DIR / "trajectories",
    '5080': BASE_DIR / "vast_downloads" / "5080",
    '5090': BASE_DIR / "vast_downloads" / "5090",
    'new_5090': BASE_DIR / "vast_downloads" / "new_5090",
}

SIMULATIONS = [
    ("G908R", "febuxostat", 1, "trajectories"),
    ("G908R", "febuxostat", 2, "5080"),
    ("G908R", "febuxostat", 3, "5080"),
    ("G908R", "natural", 1, "5090"),
    ("G908R", "natural", 2, "new_5090"),
    ("G908R", "natural", 3, "new_5090"),
    ("R702W", "febuxostat", 1, "trajectories"),
    ("R702W", "febuxostat", 2, "trajectories"),
    ("R702W", "febuxostat", 3, "trajectories"),
    ("R702W", "natural", 1, "5080"),
    ("R702W", "natural", 2, "5080"),
    ("R702W", "natural", 3, "new_5090"),
]

def get_file_paths(mutant, ligand, rep, folder_key):
    folder = TRAJECTORY_LOCATIONS[folder_key]
    base_name = f"{mutant}_{ligand}_rep{rep}"
    return folder / f"{base_name}.dcd", folder / f"{base_name}_solvated.pdb"

def analyze_contacts(pdb_path, dcd_path, stride=1):
    """Analyze protein-ligand contacts."""
    print(f"  Loading: {pdb_path.name}")
    u = mda.Universe(str(pdb_path), str(dcd_path))

    ligand = u.select_atoms(LIGAND_SELECTION)
    protein = u.select_atoms(PROTEIN_SELECTION)

    print(f"  Ligand atoms: {len(ligand)} | Protein atoms: {len(protein)}")

    if len(ligand) == 0 or len(protein) == 0:
        return None, None

    # Get unique residues in protein
    protein_residues = {}
    for atom in protein:
        resid = atom.resid
        if resid not in protein_residues:
            protein_residues[resid] = {'resname': atom.resname, 'atoms': []}
        protein_residues[resid]['atoms'].append(atom.index)

    # Track contacts per residue
    contact_counts = defaultdict(int)
    total_frames = 0
    frame_count = 0

    # Read frames sequentially (DCD files have issues with random access)
    for ts in u.trajectory:
        frame_count += 1
        if (frame_count - 1) % stride != 0:
            continue

        total_frames += 1

        # Get ligand positions
        lig_pos = ligand.positions

        # Check each residue
        for resid, info in protein_residues.items():
            res_atoms = u.atoms[info['atoms']]
            res_pos = res_atoms.positions

            # Calculate minimum distance between ligand and this residue
            # Use box=u.dimensions for PBC minimum image convention
            dist_matrix = distances.distance_array(lig_pos, res_pos, box=u.dimensions)
            min_dist = dist_matrix.min()

            if min_dist < CONTACT_CUTOFF:
                contact_counts[resid] += 1

        if total_frames % 500 == 0:
            print(f"    Processed {total_frames} frames...")

    # Calculate frequencies
    contact_frequencies = {}
    for resid, count in contact_counts.items():
        freq = (count / total_frames) * 100
        resname = protein_residues[resid]['resname']
        contact_frequencies[resid] = {'resname': resname, 'frequency': freq, 'frames': count}

    return contact_frequencies, total_frames

def main():
    print("=" * 80)
    print("PHASE A1 ANALYSIS - PROTEIN-LIGAND CONTACTS")
    print("=" * 80)
    print(f"Contact cutoff: {CONTACT_CUTOFF} Å")
    print(f"Key residues: {KEY_RESIDUES}")
    print()

    analysis_dir = BASE_DIR / "analysis"
    analysis_dir.mkdir(exist_ok=True)

    all_results = []
    key_residue_data = []

    for i, (mutant, ligand, rep, folder_key) in enumerate(SIMULATIONS, 1):
        print(f"[{i}/12] {mutant}_{ligand}_rep{rep}")
        dcd_path, pdb_path = get_file_paths(mutant, ligand, rep, folder_key)

        if not dcd_path.exists() or not pdb_path.exists():
            print(f"  ERROR: Files not found")
            # Still write 0% for key residues when files not found
            for key_res in KEY_RESIDUES:
                key_residue_data.append({
                    'mutant': mutant, 'ligand': ligand, 'rep': rep,
                    'resid': key_res, 'resname': 'FILE_NOT_FOUND', 'frequency': 0.0
                })
            print()
            continue

        try:
            contacts, total_frames = analyze_contacts(pdb_path, dcd_path)

            if contacts:
                # Get persistent contacts (>50%)
                persistent = {k: v for k, v in contacts.items() if v['frequency'] > 50}
                print(f"  Total contacts: {len(contacts)} | Persistent (>50%): {len(persistent)}")

                # Store all contacts
                for resid, data in contacts.items():
                    all_results.append({
                        'mutant': mutant, 'ligand': ligand, 'rep': rep,
                        'resid': resid, 'resname': data['resname'],
                        'frequency': data['frequency'], 'frames': data['frames'],
                        'total_frames': total_frames
                    })

                # Print top 10 contacts
                sorted_contacts = sorted(contacts.items(), key=lambda x: x[1]['frequency'], reverse=True)[:10]
                contact_strs = [f"{c[1]['resname']}{c[0]}({c[1]['frequency']:.0f}%)" for c in sorted_contacts]
                print(f"  Top contacts: {', '.join(contact_strs)}")
            else:
                print(f"  No contacts found")

            # ALWAYS track key residues (even if no contacts found)
            for key_res in KEY_RESIDUES:
                if contacts and key_res in contacts:
                    freq = contacts[key_res]['frequency']
                    resname = contacts[key_res]['resname']
                else:
                    freq = 0.0
                    resname = "N/A"

                key_residue_data.append({
                    'mutant': mutant, 'ligand': ligand, 'rep': rep,
                    'resid': key_res, 'resname': resname, 'frequency': freq
                })

        except Exception as e:
            print(f"  ERROR: {e}")
            # Still write 0% for key residues on error
            for key_res in KEY_RESIDUES:
                key_residue_data.append({
                    'mutant': mutant, 'ligand': ligand, 'rep': rep,
                    'resid': key_res, 'resname': 'ERROR', 'frequency': 0.0
                })

        print()

    # Save full contact frequencies
    csv_path = analysis_dir / "contact_frequencies.csv"
    with open(csv_path, 'w') as f:
        f.write("mutant,ligand,rep,resid,resname,frequency_percent,contact_frames,total_frames\n")
        for r in all_results:
            f.write(f"{r['mutant']},{r['ligand']},{r['rep']},{r['resid']},{r['resname']},{r['frequency']:.2f},{r['frames']},{r['total_frames']}\n")
    print(f"Full contacts saved to: {csv_path}")

    # Print key residue summary
    print("\n" + "=" * 80)
    print("KEY RESIDUE CONTACT FREQUENCIES (%)")
    print("=" * 80)

    # Aggregate by mutant/ligand
    print(f"\n{'Mutant':<8} {'Ligand':<12} {'Rep':<4} {'GLU1008':<10} {'ASP1011':<10} {'ARG1037':<10}")
    print("-" * 60)

    for i, (mutant, ligand, rep, folder_key) in enumerate(SIMULATIONS):
        row_data = [d for d in key_residue_data if d['mutant'] == mutant and d['ligand'] == ligand and d['rep'] == rep]
        if row_data:
            vals = {d['resid']: d['frequency'] for d in row_data}
            print(f"{mutant:<8} {ligand:<12} {rep:<4} {vals.get(1008, 0):<10.1f} {vals.get(1011, 0):<10.1f} {vals.get(1037, 0):<10.1f}")

    # Save key residue data
    key_csv_path = analysis_dir / "key_residue_contacts.csv"
    with open(key_csv_path, 'w') as f:
        f.write("mutant,ligand,rep,resid,resname,frequency_percent\n")
        for r in key_residue_data:
            f.write(f"{r['mutant']},{r['ligand']},{r['rep']},{r['resid']},{r['resname']},{r['frequency']:.2f}\n")
    print(f"\nKey residue data saved to: {key_csv_path}")

    print("\nContact analysis complete!")

if __name__ == "__main__":
    main()
