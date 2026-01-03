#!/usr/bin/env python3
"""
Download fresh NOD2 structure from AlphaFold and prepare for MD.
"""

import urllib.request
import os

# Human NOD2 AlphaFold structure
URL = "https://alphafold.ebi.ac.uk/files/AF-Q9HC29-F1-model_v4.pdb"
STRUCTURES_DIR = "structures"
RAW_OUTPUT = f"{STRUCTURES_DIR}/NOD2_alphafold.pdb"
CLEAN_OUTPUT = f"{STRUCTURES_DIR}/NOD2_LRR_clean.pdb"

def main():
    os.makedirs(STRUCTURES_DIR, exist_ok=True)

    # Download from AlphaFold
    print("=" * 60)
    print("DOWNLOADING NOD2 FROM ALPHAFOLD")
    print("=" * 60)
    print(f"URL: {URL}")
    print(f"Output: {RAW_OUTPUT}")

    urllib.request.urlretrieve(URL, RAW_OUTPUT)
    print("Download complete!")

    # Fix and clean structure
    print("\nPreparing structure with PDBFixer...")
    from pdbfixer import PDBFixer
    from openmm.app import PDBFile

    fixer = PDBFixer(filename=RAW_OUTPUT)

    # Find and fix any issues
    fixer.findMissingResidues()
    print(f"  Missing residues: {len(fixer.missingResidues)}")

    fixer.findMissingAtoms()
    print(f"  Missing atoms: {sum(len(atoms) for atoms in fixer.missingAtoms.values())}")
    fixer.addMissingAtoms()

    fixer.addMissingHydrogens(7.0)
    print(f"  Added hydrogens at pH 7.0")

    # Save clean structure
    with open(CLEAN_OUTPUT, 'w') as f:
        PDBFile.writeFile(fixer.topology, fixer.positions, f)

    print(f"\nSaved: {CLEAN_OUTPUT}")
    print(f"Atoms: {fixer.topology.getNumAtoms()}")

    print("\n" + "=" * 60)
    print("CLEAN RECEPTOR READY!")
    print("=" * 60)

if __name__ == "__main__":
    main()
