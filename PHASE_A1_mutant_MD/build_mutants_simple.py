#!/usr/bin/env python3
"""
PHASE A1: Build NOD2 Mutant Structures (R702W and G908R) - Simple Version

Creates mutant structures by modifying residue names.
Missing sidechain atoms will be added during MD setup (OpenMM/Modeller).
"""

import os
import sys

# Paths
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PHASE6_DIR = os.path.join(os.path.dirname(SCRIPT_DIR), "PHASE_6")
STRUCTURES_DIR = os.path.join(SCRIPT_DIR, "structures")

WT_PDB = os.path.join(PHASE6_DIR, "structures", "NOD2_LRR_clean.pdb")

# Output files
R702W_PDB = os.path.join(STRUCTURES_DIR, "NOD2_R702W.pdb")
G908R_PDB = os.path.join(STRUCTURES_DIR, "NOD2_G908R.pdb")
VERIFICATION_LOG = os.path.join(STRUCTURES_DIR, "verification_log.txt")

# Amino acid 3-letter codes
AA_CODES = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
    'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
    'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
    'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'
}


def verify_wt_structure(pdb_path):
    """Verify the wild-type structure has expected residues."""
    verification = {
        'pdb_path': pdb_path,
        'chain_id': None,
        'residue_702': None,
        'residue_908': None,
        'residue_1008': None,  # Pocket residue
        'residue_1011': None,  # Pocket residue
        'residue_1037': None,  # Pocket residue
        'insertion_codes': [],
        'total_residues': 0,
        'success': False
    }

    seen_residues = set()

    with open(pdb_path, 'r') as f:
        for line in f:
            if line.startswith('ATOM'):
                chain = line[21]
                res_num_str = line[22:26].strip()
                ins_code = line[26].strip()
                res_name = line[17:20].strip()

                # Check for insertion codes
                if ins_code:
                    ic_entry = f"{res_num_str}{ins_code}"
                    if ic_entry not in verification['insertion_codes']:
                        verification['insertion_codes'].append(ic_entry)

                try:
                    res_num = int(res_num_str)
                except ValueError:
                    continue

                if chain == 'A':
                    verification['chain_id'] = 'A'
                    seen_residues.add(res_num)

                    if res_num == 702:
                        verification['residue_702'] = res_name
                    if res_num == 908:
                        verification['residue_908'] = res_name
                    if res_num == 1008:
                        verification['residue_1008'] = res_name
                    if res_num == 1011:
                        verification['residue_1011'] = res_name
                    if res_num == 1037:
                        verification['residue_1037'] = res_name

    verification['total_residues'] = len(seen_residues)

    # Check expected values
    if (verification['chain_id'] == 'A' and
        verification['residue_702'] == 'ARG' and
        verification['residue_908'] == 'GLY'):
        verification['success'] = True

    return verification


def create_mutant_pdb(wt_pdb, output_pdb, chain_id, res_num, old_res, new_res, mutation_name):
    """
    Create a mutant PDB by replacing residue name and removing sidechain atoms.

    Parameters:
    - wt_pdb: Path to wild-type PDB
    - output_pdb: Path for output mutant PDB
    - chain_id: Chain containing the residue
    - res_num: Residue number to mutate
    - old_res: Original residue 3-letter code (for verification)
    - new_res: New residue 3-letter code
    - mutation_name: Name for logging (e.g., "R702W")
    """
    print(f"\nCreating {mutation_name} mutant...")
    print(f"  Mutation: {old_res}{res_num} -> {new_res}{res_num}")

    # Atoms to keep for the mutated residue (backbone only)
    # Sidechain atoms will be rebuilt during MD setup
    backbone_atoms = {'N', 'CA', 'C', 'O'}

    modified_lines = []
    atoms_removed = 0
    residue_found = False

    with open(wt_pdb, 'r') as f:
        for line in f:
            if line.startswith('ATOM'):
                chain = line[21]
                res_num_str = line[22:26].strip()
                atom_name = line[12:16].strip()
                res_name = line[17:20].strip()

                try:
                    curr_res_num = int(res_num_str)
                except ValueError:
                    modified_lines.append(line)
                    continue

                # Check if this is our target residue
                if chain == chain_id and curr_res_num == res_num:
                    residue_found = True

                    # Verify it's the expected residue
                    if res_name != old_res:
                        print(f"  WARNING: Expected {old_res} at position {res_num}, found {res_name}")

                    # Keep only backbone atoms
                    if atom_name in backbone_atoms:
                        # Replace residue name in the line
                        new_line = line[:17] + new_res.ljust(3) + line[20:]
                        modified_lines.append(new_line)
                    else:
                        atoms_removed += 1
                else:
                    modified_lines.append(line)
            elif line.startswith('TER'):
                modified_lines.append(line)
            elif line.startswith('END'):
                modified_lines.append(line)
            elif line.startswith('CRYST') or line.startswith('REMARK'):
                modified_lines.append(line)
            else:
                modified_lines.append(line)

    if not residue_found:
        print(f"  ERROR: Residue {res_num} not found on chain {chain_id}!")
        return None

    # Write output PDB
    with open(output_pdb, 'w') as f:
        # Add header comment about mutation
        f.write(f"REMARK   1 MUTANT STRUCTURE: {mutation_name}\n")
        f.write(f"REMARK   2 MUTATION: {old_res}{res_num}{new_res[0]} on chain {chain_id}\n")
        f.write(f"REMARK   3 SIDECHAIN ATOMS REMOVED - WILL BE REBUILT DURING MD SETUP\n")
        f.writelines(modified_lines)

    print(f"  Sidechain atoms removed: {atoms_removed}")
    print(f"  Saved: {output_pdb}")

    return output_pdb


def write_verification_log(verification, mutations_done):
    """Write detailed verification log."""
    with open(VERIFICATION_LOG, 'w') as f:
        f.write("=" * 60 + "\n")
        f.write("PHASE A1: Mutant Structure Verification Log\n")
        f.write("=" * 60 + "\n\n")

        f.write("WILD-TYPE STRUCTURE VERIFICATION\n")
        f.write("-" * 40 + "\n")
        f.write(f"PDB file: {verification['pdb_path']}\n")
        f.write(f"Chain ID used: {verification['chain_id']}\n")
        f.write(f"Total residues: {verification['total_residues']}\n\n")

        f.write("Mutation Site Verification:\n")
        status_702 = "[OK]" if verification['residue_702'] == 'ARG' else "[X] MISMATCH"
        status_908 = "[OK]" if verification['residue_908'] == 'GLY' else "[X] MISMATCH"
        f.write(f"  Original residue 702: {verification['residue_702']} (expected ARG) {status_702}\n")
        f.write(f"  Original residue 908: {verification['residue_908']} (expected GLY) {status_908}\n\n")

        f.write("Binding Pocket Residues (for reference):\n")
        f.write(f"  Residue 1008: {verification['residue_1008']} (GLU in WT)\n")
        f.write(f"  Residue 1011: {verification['residue_1011']} (ASP in WT)\n")
        f.write(f"  Residue 1037: {verification['residue_1037']} (ARG in WT)\n\n")

        ic_str = ', '.join(verification['insertion_codes']) if verification['insertion_codes'] else 'None'
        f.write(f"Insertion codes: {ic_str}\n\n")

        f.write(f"Overall verification: {'PASSED' if verification['success'] else 'FAILED'}\n\n")

        f.write("=" * 60 + "\n")
        f.write("MUTATIONS CREATED\n")
        f.write("=" * 60 + "\n\n")

        for mut_name, mut_file in mutations_done.items():
            if mut_file:
                f.write(f"{mut_name}:\n")
                f.write(f"  Output file: {mut_file}\n")
                if mut_name == 'R702W':
                    f.write("  Change: ARG (R) -> TRP (W) at position 702\n")
                    f.write("  Notes: R702W is a Crohn's-associated variant (~10% of patients)\n")
                    f.write("  Distance from pocket: 79.4 Å (distal, should not affect binding)\n")
                elif mut_name == 'G908R':
                    f.write("  Change: GLY (G) -> ARG (R) at position 908\n")
                    f.write("  Notes: G908R is a Crohn's-associated variant (~4.6% of patients)\n")
                    f.write("  Distance from pocket: 34.0 Å (distal, should not affect binding)\n")
                f.write("\n")

        f.write("=" * 60 + "\n")
        f.write("TECHNICAL NOTES\n")
        f.write("=" * 60 + "\n\n")
        f.write("1. Backbone atoms (N, CA, C, O) preserved from wild-type\n")
        f.write("2. Sidechain atoms removed (will be rebuilt during MD setup)\n")
        f.write("3. OpenMM/PDBFixer will add missing atoms during system preparation\n")
        f.write("4. Ligand pose from Phase 6 WT simulations will be transferred\n")

    print(f"\nVerification log saved: {VERIFICATION_LOG}")


def main():
    print("=" * 60)
    print("PHASE A1: Building NOD2 Mutant Structures")
    print("=" * 60)

    # Create output directory
    os.makedirs(STRUCTURES_DIR, exist_ok=True)

    # Verify WT structure first
    print("\n1. Verifying wild-type structure...")
    if not os.path.exists(WT_PDB):
        print(f"ERROR: Wild-type PDB not found: {WT_PDB}")
        sys.exit(1)

    verification = verify_wt_structure(WT_PDB)

    print(f"\n  Structure: {os.path.basename(WT_PDB)}")
    print(f"  Chain ID: {verification['chain_id']}")
    print(f"  Total residues: {verification['total_residues']}")
    print(f"\n  Residue 702: {verification['residue_702']} (expected: ARG)")
    print(f"  Residue 908: {verification['residue_908']} (expected: GLY)")
    print(f"\n  Pocket residues:")
    print(f"    1008: {verification['residue_1008']}")
    print(f"    1011: {verification['residue_1011']}")
    print(f"    1037: {verification['residue_1037']}")

    if verification['insertion_codes']:
        print(f"\n  Insertion codes found: {verification['insertion_codes']}")
    else:
        print(f"\n  Insertion codes: None")

    if not verification['success']:
        print("\n*** VERIFICATION FAILED ***")
        print("Residue numbering or identity does not match expected values.")
        write_verification_log(verification, {})
        sys.exit(1)

    print("\n  [OK] Verification PASSED")

    # Create mutant structures
    print("\n" + "=" * 60)
    print("2. Creating mutant structures...")
    print("=" * 60)

    mutations_done = {}

    # R702W: ARG -> TRP
    output_file = create_mutant_pdb(
        WT_PDB, R702W_PDB,
        chain_id='A', res_num=702,
        old_res='ARG', new_res='TRP',
        mutation_name='R702W'
    )
    mutations_done['R702W'] = output_file

    # G908R: GLY -> ARG
    output_file = create_mutant_pdb(
        WT_PDB, G908R_PDB,
        chain_id='A', res_num=908,
        old_res='GLY', new_res='ARG',
        mutation_name='G908R'
    )
    mutations_done['G908R'] = output_file

    # Write verification log
    print("\n" + "=" * 60)
    print("3. Writing verification log...")
    print("=" * 60)
    write_verification_log(verification, mutations_done)

    print("\n" + "=" * 60)
    print("CHECKPOINT 1 COMPLETE")
    print("=" * 60)
    print(f"\nOutput files:")
    for name, path in mutations_done.items():
        if path:
            print(f"  [OK] {name}: {path}")
    print(f"  [OK] Log: {VERIFICATION_LOG}")

    print("\nNEXT STEP: Transfer ligand poses from Phase 6 WT structures")


if __name__ == "__main__":
    main()
