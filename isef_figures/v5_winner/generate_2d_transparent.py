#!/usr/bin/env python3
"""Generate 2D ligand structures with TRANSPARENT backgrounds."""

from rdkit import Chem
from rdkit.Chem import Draw, AllChem
from rdkit.Chem.Draw import rdMolDraw2D
from pathlib import Path

OUTPUT_DIR = Path(r"C:\Users\vasud\nod2-screening-data\isef_figures\v5_winner\ppt_photos\transparent_bg")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

LIGANDS = {
    "febuxostat": {
        "smiles": "Cc1nc(-c2ccc(OCC(C)C)c(C#N)c2)sc1C(=O)O",
        "name": "Febuxostat",
    },
    "bufadienolide": {
        "smiles": "C[C@]12CCC3C([C@]1(CC[C@@H]2C4=COC(=O)C=C4)O)CC[C@]5([C@@]3(CC[C@@H](C5)O)C)O",
        "name": "Bufadienolide",
    }
}

def generate_transparent_2d(smiles, name, output_path, size=(600, 600)):
    """Generate 2D structure with transparent background."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        print(f"ERROR: Could not parse SMILES for {name}")
        return False
    
    AllChem.Compute2DCoords(mol)
    
    # Use Cairo drawer for PNG with alpha channel
    drawer = rdMolDraw2D.MolDraw2DCairo(size[0], size[1])
    
    # Set transparent background (alpha = 0)
    opts = drawer.drawOptions()
    opts.clearBackground = False  # Don't fill background
    opts.bondLineWidth = 3.0
    
    drawer.DrawMolecule(mol)
    drawer.FinishDrawing()
    
    with open(output_path, 'wb') as f:
        f.write(drawer.GetDrawingText())
    
    print(f"Saved: {output_path}")
    return True

def main():
    print("=" * 60)
    print("GENERATING TRANSPARENT 2D LIGAND STRUCTURES")
    print("=" * 60)
    
    for key, data in LIGANDS.items():
        output_path = OUTPUT_DIR / f"{key}_transparent.png"
        print(f"\n{data['name']}:")
        generate_transparent_2d(data['smiles'], data['name'], output_path)
    
    print("\n" + "=" * 60)
    print("DONE - Files in transparent_bg/")
    print("=" * 60)

if __name__ == "__main__":
    main()
