#!/usr/bin/env python3
"""
Fix corrupted receptor PDB - run ONCE before simulations.
Removes duplicates, fixes missing atoms, adds hydrogens.
"""

from pdbfixer import PDBFixer
from openmm.app import PDBFile
import os

STRUCTURES_DIR = "structures"

# List of complex PDBs to fix
COMPLEX_PDBS = [
    "complex_febuxostat.pdb",
    "complex_ursodiol.pdb",
    "complex_budesonide.pdb",
    "complex_natural_cid10592.pdb",
    "complex_apo.pdb",
    "complex_decoy.pdb",
]

def fix_pdb(input_file, output_file):
    """Fix a single PDB file."""
    print(f"\nFixing: {input_file}")

    fixer = PDBFixer(filename=input_file)

    # Remove heterogens (ligands, waters) - we'll add ligand back separately
    fixer.removeHeterogens(keepWater=False)
    print("  Removed heterogens")

    # Find and fix missing residues
    fixer.findMissingResidues()
    print(f"  Missing residues: {len(fixer.missingResidues)}")

    # Find and replace nonstandard residues
    fixer.findNonstandardResidues()
    print(f"  Nonstandard residues: {len(fixer.nonstandardResidues)}")
    fixer.replaceNonstandardResidues()

    # Find and add missing atoms
    fixer.findMissingAtoms()
    print(f"  Missing atoms: {sum(len(atoms) for atoms in fixer.missingAtoms.values())}")
    fixer.addMissingAtoms()

    # Add hydrogens at pH 7.0
    fixer.addMissingHydrogens(7.0)

    # Save fixed PDB
    with open(output_file, 'w') as f:
        PDBFile.writeFile(fixer.topology, fixer.positions, f)

    print(f"  Saved: {output_file} ({fixer.topology.getNumAtoms()} atoms)")

def main():
    print("=" * 60)
    print("FIXING RECEPTOR PDB FILES")
    print("=" * 60)

    os.makedirs(STRUCTURES_DIR, exist_ok=True)

    for pdb in COMPLEX_PDBS:
        input_path = f"{STRUCTURES_DIR}/{pdb}"
        output_path = f"{STRUCTURES_DIR}/{pdb.replace('.pdb', '_fixed.pdb')}"

        if os.path.exists(input_path):
            fix_pdb(input_path, output_path)
        else:
            print(f"\nSkipping (not found): {input_path}")

    print("\n" + "=" * 60)
    print("DONE! Use *_fixed.pdb files for simulations.")
    print("=" * 60)

if __name__ == "__main__":
    main()
