#!/usr/bin/env python3
"""
Validate AlphaFold NOD2 structure against experimental PDB 5IRM.

Calculates backbone CA RMSD after structural alignment.
Output: receptor/alphafold_validation.json

Note: AlphaFold structure is renumbered 1-1040, while PDB 5IRM uses original
UniProt numbering (193-1020). We align by superimposing all CA atoms using
least-squares fitting, not by residue number matching.

Author: Tanmay Vasudeva
"""

import json
import os
import urllib.request
from pathlib import Path

try:
    from Bio.PDB import PDBParser, Superimposer
    from Bio.PDB.Polypeptide import is_aa
except ImportError:
    print("Installing biopython...")
    os.system("pip install biopython")
    from Bio.PDB import PDBParser, Superimposer
    from Bio.PDB.Polypeptide import is_aa

import numpy as np

# Paths
BASE_DIR = Path(__file__).parent
RECEPTOR_DIR = BASE_DIR / "receptor"
RECEPTOR_DIR.mkdir(exist_ok=True)

ALPHAFOLD_PDB = BASE_DIR / "cross_validation" / "human_NOD2_LRR.pdb"
PDB_5IRM_PATH = RECEPTOR_DIR / "5IRM.pdb"
OUTPUT_JSON = RECEPTOR_DIR / "alphafold_validation.json"


def download_pdb(pdb_id, output_path):
    """Download PDB file from RCSB."""
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    print(f"Downloading {pdb_id} from RCSB...")
    urllib.request.urlretrieve(url, output_path)
    print(f"Saved to {output_path}")


def get_ca_atoms(structure, chain_id=None):
    """Extract CA atoms from structure."""
    ca_atoms = []
    residue_info = []

    for model in structure:
        for chain in model:
            if chain_id and chain.id != chain_id:
                continue
            for residue in chain:
                if is_aa(residue) and 'CA' in residue:
                    ca_atoms.append(residue['CA'])
                    residue_info.append({
                        'chain': chain.id,
                        'resnum': residue.id[1],
                        'resname': residue.resname
                    })

    return ca_atoms, residue_info


def get_sequence_from_atoms(ca_atoms, residue_info):
    """Get one-letter sequence from residue info."""
    three_to_one = {
        'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F',
        'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LYS': 'K', 'LEU': 'L',
        'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R',
        'SER': 'S', 'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y'
    }
    seq = ''
    for r in residue_info:
        seq += three_to_one.get(r['resname'], 'X')
    return seq


def find_best_alignment_window(af_ca, exp_ca, window_size=200):
    """
    Find the best overlapping window between structures.
    Since structures may cover different regions, find the window
    with lowest RMSD after superposition.
    """
    best_rmsd = float('inf')
    best_af_start = 0
    best_exp_start = 0

    # Try different starting positions
    for af_start in range(0, len(af_ca) - window_size + 1, 50):
        for exp_start in range(0, len(exp_ca) - window_size + 1, 50):
            af_subset = af_ca[af_start:af_start + window_size]
            exp_subset = exp_ca[exp_start:exp_start + window_size]

            if len(af_subset) != len(exp_subset):
                continue

            try:
                sup = Superimposer()
                sup.set_atoms(exp_subset, af_subset)
                if sup.rms < best_rmsd:
                    best_rmsd = sup.rms
                    best_af_start = af_start
                    best_exp_start = exp_start
            except:
                continue

    return best_af_start, best_exp_start, best_rmsd


def main():
    print("=" * 60)
    print("ALPHAFOLD STRUCTURE VALIDATION")
    print("Comparing AlphaFold NOD2-LRR to experimental PDB 5IRM")
    print("=" * 60)

    # Download PDB 5IRM if not exists
    if not PDB_5IRM_PATH.exists():
        download_pdb("5IRM", PDB_5IRM_PATH)

    # Check AlphaFold structure exists
    if not ALPHAFOLD_PDB.exists():
        print(f"ERROR: AlphaFold structure not found at {ALPHAFOLD_PDB}")
        return

    # Parse structures
    parser = PDBParser(QUIET=True)

    print(f"\nLoading AlphaFold: {ALPHAFOLD_PDB}")
    alphafold = parser.get_structure("alphafold", ALPHAFOLD_PDB)

    print(f"Loading PDB 5IRM: {PDB_5IRM_PATH}")
    experimental = parser.get_structure("5IRM", PDB_5IRM_PATH)

    # Get CA atoms - only chain A for 5IRM
    af_ca, af_res = get_ca_atoms(alphafold)
    exp_ca, exp_res = get_ca_atoms(experimental, chain_id='A')

    print(f"\nAlphaFold CA atoms: {len(af_ca)}")
    print(f"PDB 5IRM chain A CA atoms: {len(exp_ca)}")

    # The structures have different numbering and potentially cover different regions
    # We'll use all atoms and let superimposer find best fit

    # Use the smaller structure's length
    n_atoms = min(len(af_ca), len(exp_ca))

    print(f"\nUsing {n_atoms} residues for alignment (smaller structure)")

    # Take contiguous blocks from each structure
    # AlphaFold starts at residue 1, 5IRM starts at ~193
    # The LRR domain in both should have similar fold

    # Try to find best overlapping window
    print("\nFinding best structural alignment...")
    window_size = min(300, n_atoms)
    af_start, exp_start, window_rmsd = find_best_alignment_window(af_ca, exp_ca, window_size)

    print(f"Best alignment: AF[{af_start}:{af_start+window_size}] vs 5IRM[{exp_start}:{exp_start+window_size}]")
    print(f"Window RMSD ({window_size} residues): {window_rmsd:.2f} A")

    # Now do full alignment using the discovered offset
    # Determine the offset and align maximum residues
    offset = exp_start - af_start

    if offset >= 0:
        # 5IRM starts later
        af_aligned = af_ca[:len(exp_ca)-offset] if offset > 0 else af_ca[:len(exp_ca)]
        exp_aligned = exp_ca[offset:offset+len(af_aligned)] if offset > 0 else exp_ca[:len(af_ca)]
    else:
        # AlphaFold starts later
        exp_aligned = exp_ca[:len(af_ca)+offset]
        af_aligned = af_ca[-offset:-offset+len(exp_aligned)]

    # Ensure equal length
    min_len = min(len(af_aligned), len(exp_aligned))
    af_aligned = af_aligned[:min_len]
    exp_aligned = exp_aligned[:min_len]

    print(f"\nAligning {len(af_aligned)} residues...")

    # Final superposition
    super_imposer = Superimposer()
    super_imposer.set_atoms(exp_aligned, af_aligned)
    rmsd = super_imposer.rms

    print(f"\n{'='*60}")
    print(f"RESULT: Backbone CA RMSD = {rmsd:.2f} A")
    print(f"(over {len(af_aligned)} aligned residues)")
    print(f"{'='*60}")

    # Interpretation
    if rmsd < 1.0:
        quality = "EXCELLENT - Nearly identical structures"
    elif rmsd < 2.0:
        quality = "GOOD - Minor structural differences"
    elif rmsd < 3.0:
        quality = "MODERATE - Noticeable differences"
    else:
        quality = "POOR - Significant structural deviation"

    print(f"Interpretation: {quality}")

    # Save results
    results = {
        "alphafold_structure": str(ALPHAFOLD_PDB),
        "experimental_structure": "PDB 5IRM chain A",
        "alphafold_total_residues": len(af_ca),
        "experimental_total_residues": len(exp_ca),
        "aligned_residues": len(af_aligned),
        "backbone_ca_rmsd_angstrom": round(rmsd, 2),
        "quality": quality,
        "notes": f"Best alignment window ({window_size} res) RMSD: {window_rmsd:.2f} A"
    }

    with open(OUTPUT_JSON, 'w') as f:
        json.dump(results, f, indent=2)

    print(f"\nResults saved to: {OUTPUT_JSON}")

    return rmsd


if __name__ == "__main__":
    rmsd = main()
