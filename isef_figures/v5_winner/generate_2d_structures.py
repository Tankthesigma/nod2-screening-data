#!/usr/bin/env python3
"""
Generate 2D ligand structure images for ISEF poster.
"""

from rdkit import Chem
from rdkit.Chem import Draw, AllChem
from rdkit.Chem.Draw import rdMolDraw2D
from pathlib import Path

# Output directory
OUTPUT_DIR = Path(r"C:\Users\vasud\nod2-screening-data\isef_figures\v5_winner\ppt_photos")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# Ligand data
LIGANDS = {
    "febuxostat": {
        "smiles": "Cc1nc(-c2ccc(OCC(C)C)c(C#N)c2)sc1C(=O)O",
        "name": "Febuxostat",
        "mw": "316.37",
        "formula": "C16H16N2O3S"
    },
    "bufadienolide": {
        "smiles": "C[C@]12CCC3C([C@]1(CC[C@@H]2C4=COC(=O)C=C4)O)CC[C@]5([C@@]3(CC[C@@H](C5)O)C)O",
        "name": "Bufadienolide (CID 10120)",
        "mw": "402.53",
        "formula": "C24H34O5"
    }
}

def generate_2d_image(smiles, name, output_path, size=(500, 500)):
    """Generate high-quality 2D structure image."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        print(f"ERROR: Could not parse SMILES for {name}")
        return False

    # Generate 2D coordinates
    AllChem.Compute2DCoords(mol)

    # Create drawer with white background
    drawer = rdMolDraw2D.MolDraw2DCairo(size[0], size[1])

    # Set drawing options
    opts = drawer.drawOptions()
    opts.backgroundColour = (1, 1, 1)  # White background
    opts.bondLineWidth = 2.5

    # Draw molecule
    drawer.DrawMolecule(mol)
    drawer.FinishDrawing()

    # Save to file
    with open(output_path, 'wb') as f:
        f.write(drawer.GetDrawingText())

    print(f"Saved: {output_path}")
    return True

def main():
    print("=" * 60)
    print("GENERATING 2D LIGAND STRUCTURES")
    print("=" * 60)

    for key, data in LIGANDS.items():
        output_path = OUTPUT_DIR / f"ligand_2d_{key}.png"
        print(f"\n{data['name']}:")
        print(f"  Formula: {data['formula']}")
        print(f"  MW: {data['mw']} g/mol")
        print(f"  SMILES: {data['smiles'][:50]}...")

        success = generate_2d_image(
            data['smiles'],
            data['name'],
            output_path,
            size=(500, 500)
        )

        if success:
            print(f"  Status: SUCCESS")
        else:
            print(f"  Status: FAILED")

    print("\n" + "=" * 60)
    print("DONE")
    print("=" * 60)

if __name__ == "__main__":
    main()
