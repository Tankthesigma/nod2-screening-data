#!/usr/bin/env python3
"""
PHASE A1: Build NOD2 Mutant Structures (R702W and G908R)

Creates mutant structures from wild-type NOD2 LRR for MD simulations.
Uses PDBFixer for mutagenesis and sidechain optimization.
"""

import os
import sys

# Try to import required libraries
try:
    from pdbfixer import PDBFixer
    from openmm.app import PDBFile
    from openmm import unit
    print("PDBFixer and OpenMM loaded successfully")
except ImportError as e:
    print(f"ERROR: Required library not found: {e}")
    print("Please install with: conda install -c conda-forge pdbfixer openmm")
    sys.exit(1)

# Paths
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PHASE6_DIR = os.path.join(os.path.dirname(SCRIPT_DIR), "PHASE_6")
STRUCTURES_DIR = os.path.join(SCRIPT_DIR, "structures")

WT_PDB = os.path.join(PHASE6_DIR, "structures", "NOD2_LRR_clean.pdb")

# Output files
R702W_PDB = os.path.join(STRUCTURES_DIR, "NOD2_R702W.pdb")
G908R_PDB = os.path.join(STRUCTURES_DIR, "NOD2_G908R.pdb")
VERIFICATION_LOG = os.path.join(STRUCTURES_DIR, "verification_log.txt")

# Mutations to make (chain, residue_id, new_residue)
MUTATIONS = {
    'R702W': ('A', 702, 'TRP'),  # Arginine -> Tryptophan
    'G908R': ('A', 908, 'ARG'),  # Glycine -> Arginine
}


def verify_wt_structure(pdb_path):
    """Verify the wild-type structure has expected residues."""
    verification = {
        'pdb_path': pdb_path,
        'chain_id': None,
        'residue_702': None,
        'residue_908': None,
        'insertion_codes': [],
        'success': False
    }

    with open(pdb_path, 'r') as f:
        for line in f:
            if line.startswith('ATOM'):
                chain = line[21]
                res_num_str = line[22:26].strip()
                ins_code = line[26].strip()
                res_name = line[17:20].strip()

                # Check for insertion codes
                if ins_code:
                    verification['insertion_codes'].append(f"{res_num_str}{ins_code}")

                try:
                    res_num = int(res_num_str)
                except ValueError:
                    continue

                if chain == 'A':
                    verification['chain_id'] = 'A'
                    if res_num == 702:
                        verification['residue_702'] = res_name
                    if res_num == 908:
                        verification['residue_908'] = res_name

    # Check expected values
    if (verification['chain_id'] == 'A' and
        verification['residue_702'] == 'ARG' and
        verification['residue_908'] == 'GLY'):
        verification['success'] = True

    return verification


def create_mutant_structure(wt_pdb, output_pdb, chain_id, res_num, new_res_name, mutation_name):
    """
    Create a mutant structure using manual residue replacement.

    Since PDBFixer doesn't have direct mutagenesis, we'll:
    1. Read the PDB
    2. Replace residue name
    3. Remove old sidechain atoms
    4. Use PDBFixer to add missing atoms
    """
    print(f"\nCreating {mutation_name} mutant...")

    # Read PDB and modify
    modified_lines = []

    # Amino acid heavy atom names (backbone + CB only, sidechains will be rebuilt)
    backbone_atoms = {'N', 'CA', 'C', 'O', 'H', 'HA', 'CB'}

    with open(wt_pdb, 'r') as f:
        for line in f:
            if line.startswith('ATOM') or line.startswith('HETATM'):
                chain = line[21]
                res_num_str = line[22:26].strip()
                atom_name = line[12:16].strip()

                try:
                    curr_res_num = int(res_num_str)
                except ValueError:
                    modified_lines.append(line)
                    continue

                # Check if this is our target residue
                if chain == chain_id and curr_res_num == res_num:
                    # Keep only backbone atoms (sidechain will be rebuilt)
                    if atom_name in backbone_atoms or atom_name.startswith('H'):
                        # For backbone, keep but might need to handle hydrogen
                        if atom_name in {'N', 'CA', 'C', 'O', 'CB'}:
                            # Replace residue name
                            new_line = line[:17] + new_res_name.ljust(3) + line[20:]
                            modified_lines.append(new_line)
                    # Skip sidechain atoms - will be rebuilt
                else:
                    modified_lines.append(line)
            else:
                modified_lines.append(line)

    # Write temporary modified PDB
    temp_pdb = output_pdb + '.temp'
    with open(temp_pdb, 'w') as f:
        f.writelines(modified_lines)

    # Use PDBFixer to add missing atoms
    print(f"  Using PDBFixer to add missing sidechain atoms...")
    fixer = PDBFixer(filename=temp_pdb)
    fixer.findMissingResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(7.0)  # pH 7.0

    # Save final structure
    with open(output_pdb, 'w') as f:
        PDBFile.writeFile(fixer.topology, fixer.positions, f)

    # Clean up temp file
    os.remove(temp_pdb)

    print(f"  Saved: {output_pdb}")
    return output_pdb


def write_verification_log(verification, mutations_done):
    """Write verification log file."""
    with open(VERIFICATION_LOG, 'w') as f:
        f.write("=" * 60 + "\n")
        f.write("PHASE A1: Mutant Structure Verification Log\n")
        f.write("=" * 60 + "\n\n")

        f.write("WILD-TYPE VERIFICATION\n")
        f.write("-" * 40 + "\n")
        f.write(f"PDB file: {verification['pdb_path']}\n")
        f.write(f"Chain ID used: {verification['chain_id']}\n")
        f.write(f"Original residue 702: {verification['residue_702']}")
        f.write(" ✓\n" if verification['residue_702'] == 'ARG' else " ✗ MISMATCH!\n")
        f.write(f"Original residue 908: {verification['residue_908']}")
        f.write(" ✓\n" if verification['residue_908'] == 'GLY' else " ✗ MISMATCH!\n")
        f.write(f"Insertion codes found: {verification['insertion_codes'] if verification['insertion_codes'] else 'None'}\n")
        f.write(f"Verification status: {'PASSED' if verification['success'] else 'FAILED'}\n\n")

        f.write("MUTATIONS CREATED\n")
        f.write("-" * 40 + "\n")
        for mut_name, mut_file in mutations_done.items():
            f.write(f"{mut_name}: {mut_file}\n")

        f.write("\nMUTATION DETAILS\n")
        f.write("-" * 40 + "\n")
        f.write("R702W: ARG -> TRP at position 702 (chain A)\n")
        f.write("G908R: GLY -> ARG at position 908 (chain A)\n")

        f.write("\nNOTES\n")
        f.write("-" * 40 + "\n")
        f.write("- Sidechain atoms rebuilt using PDBFixer\n")
        f.write("- Backbone preserved from wild-type\n")
        f.write("- Ready for ligand docking pose transfer\n")

    print(f"\nVerification log saved: {VERIFICATION_LOG}")


def main():
    print("=" * 60)
    print("PHASE A1: Building NOD2 Mutant Structures")
    print("=" * 60)

    # Verify WT structure first
    print("\n1. Verifying wild-type structure...")
    if not os.path.exists(WT_PDB):
        print(f"ERROR: Wild-type PDB not found: {WT_PDB}")
        sys.exit(1)

    verification = verify_wt_structure(WT_PDB)

    print(f"  Chain ID: {verification['chain_id']}")
    print(f"  Residue 702: {verification['residue_702']} (expected: ARG)")
    print(f"  Residue 908: {verification['residue_908']} (expected: GLY)")
    print(f"  Insertion codes: {verification['insertion_codes'] if verification['insertion_codes'] else 'None'}")

    if not verification['success']:
        print("\n*** VERIFICATION FAILED ***")
        print("Residue numbering or identity does not match expected values.")
        print("Please check the PDB file manually.")
        write_verification_log(verification, {})
        sys.exit(1)

    print("\n  ✓ Verification PASSED")

    # Create mutant structures
    print("\n2. Creating mutant structures...")
    os.makedirs(STRUCTURES_DIR, exist_ok=True)

    mutations_done = {}

    # R702W
    mut_name = 'R702W'
    chain, res_num, new_res = MUTATIONS[mut_name]
    output_file = create_mutant_structure(WT_PDB, R702W_PDB, chain, res_num, new_res, mut_name)
    mutations_done[mut_name] = output_file

    # G908R
    mut_name = 'G908R'
    chain, res_num, new_res = MUTATIONS[mut_name]
    output_file = create_mutant_structure(WT_PDB, G908R_PDB, chain, res_num, new_res, mut_name)
    mutations_done[mut_name] = output_file

    # Write verification log
    print("\n3. Writing verification log...")
    write_verification_log(verification, mutations_done)

    print("\n" + "=" * 60)
    print("COMPLETE")
    print("=" * 60)
    print(f"\nOutput files:")
    print(f"  - {R702W_PDB}")
    print(f"  - {G908R_PDB}")
    print(f"  - {VERIFICATION_LOG}")


if __name__ == "__main__":
    main()
