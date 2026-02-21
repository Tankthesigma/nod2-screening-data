#!/usr/bin/env python3
"""
Generate MDP (Muramyl Dipeptide) 2D structure image for pathway diagram.
MDP is the bacterial signal that triggers NOD2 activation.
"""

try:
    from rdkit import Chem
    from rdkit.Chem import Draw, AllChem
    from pathlib import Path

    # MDP SMILES (Muramyl Dipeptide)
    # N-acetylmuramyl-L-alanyl-D-isoglutamine
    mdp_smiles = "CC(=O)N[C@@H]1[C@H](O)[C@@H](CO)O[C@@H](O)[C@@H]1OC(C)C(=O)N[C@@H](C)C(=O)N[C@@H](CCC(N)=O)C(=O)O"

    mol = Chem.MolFromSmiles(mdp_smiles)
    if mol is None:
        # Simplified version
        mdp_smiles = "CC(=O)NC1C(O)C(CO)OC(O)C1OC(C)C(=O)NC(C)C(=O)NC(CCC(N)=O)C(=O)O"
        mol = Chem.MolFromSmiles(mdp_smiles)

    if mol:
        # Generate 2D coordinates
        AllChem.Compute2DCoords(mol)

        # Draw with orange carbons
        output_path = Path(r"C:\Users\vasud\nod2-screening-data\isef_figures\pathway_diagram\MDP_3D_ISEF.png")

        img = Draw.MolToImage(mol, size=(800, 800), kekulize=True)
        img.save(str(output_path))
        print(f"Created: {output_path}")
    else:
        print("Could not parse MDP SMILES")

except ImportError:
    print("RDKit not available. Install with: conda install -c conda-forge rdkit")
    print("Or use ChemDraw/MarvinSketch to create MDP structure manually.")
