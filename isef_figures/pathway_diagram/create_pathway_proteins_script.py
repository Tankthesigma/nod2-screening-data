#!/usr/bin/env python3
"""
Create PyMOL scripts for pathway protein renders (MDP, RIPK2, TAK1, NF-kB).
These are needed for the NOD2 signaling pathway diagram on the ISEF poster.

Output: pathway_proteins.pml - Run in PyMOL to generate all renders
"""

script_content = '''# ============================================================================
# ISEF 2026 - NOD2 Signaling Pathway Protein Renders
# Tanmay Vasudeva - Texas Virtual Academy at Hallsville
# ============================================================================
#
# This script generates renders of downstream signaling proteins:
# - RIPK2 (PDB: 4C8B) - full and dim versions
# - NF-kB (PDB: 1NFI) - full and dim versions
#
# Note: MDP molecule should be created separately with RDKit or ChemDraw
# Note: TAK1 (5JH6) is a large complex, can be added if needed
#
# Run in PyMOL: @pathway_proteins.pml
# ============================================================================

# ============================================================================
# QUALITY SETTINGS
# ============================================================================
set antialias, 2
set ray_trace_mode, 1
set ray_shadows, 1
set ray_trace_fog, 0
set depth_cue, 0
set cartoon_fancy_helices, 1
set cartoon_smooth_loops, 1
set ambient, 0.4
set spec_reflect, 0.5
set spec_power, 200
set ray_opaque_background, 0

bg_color white

# ============================================================================
# RIPK2 - Full (Active pathway - bright)
# ============================================================================
reinitialize

# Fetch RIPK2 structure from PDB
fetch 4C8B, ripk2, async=0

# Clean up - remove waters and ligands
remove solvent
remove organic

# Show as cartoon
hide everything
show cartoon, ripk2

# Color forest green (matches NOD2 LRR)
color forest, ripk2

# Orient and zoom
orient ripk2
zoom ripk2, 5

# Render
ray 1200, 1200
png C:/Users/vasud/nod2-screening-data/isef_figures/pathway_diagram/RIPK2_full_ISEF.png, dpi=300

# ============================================================================
# RIPK2 - Dim (Inactive/R702W pathway - faded)
# ============================================================================
# Same structure, add transparency

set cartoon_transparency, 0.6, ripk2
color gray70, ripk2

# Render dim version
ray 1200, 1200
png C:/Users/vasud/nod2-screening-data/isef_figures/pathway_diagram/RIPK2_dim_ISEF.png, dpi=300

# ============================================================================
# NF-kB - Full (Active pathway - bright)
# ============================================================================
reinitialize

# Fetch NF-kB structure (p50-p65 heterodimer with DNA)
fetch 1NFI, nfkb, async=0

# Remove waters
remove solvent

# Show as cartoon
hide everything
show cartoon, nfkb

# Color purple
color purple, nfkb

# Remove DNA if desired (optional - keep for context)
# remove nfkb and chain C

# Orient and zoom
orient nfkb
zoom nfkb, 5

# Render
ray 1200, 1200
png C:/Users/vasud/nod2-screening-data/isef_figures/pathway_diagram/NFkB_full_ISEF.png, dpi=300

# ============================================================================
# NF-kB - Dim (Inactive/R702W pathway - faded)
# ============================================================================

set cartoon_transparency, 0.6, nfkb
color gray70, nfkb

# Render dim version
ray 1200, 1200
png C:/Users/vasud/nod2-screening-data/isef_figures/pathway_diagram/NFkB_dim_ISEF.png, dpi=300

# ============================================================================
# IKK/TAK1 - Optional (can add if needed)
# ============================================================================
# reinitialize
# fetch 5JH6, tak1, async=0
# ... similar setup ...

# ============================================================================
# DONE
# ============================================================================
print "Pathway protein renders generated!"
print "Output directory: C:/Users/vasud/nod2-screening-data/isef_figures/pathway_diagram/"
print ""
print "Generated images:"
print "  - RIPK2_full_ISEF.png"
print "  - RIPK2_dim_ISEF.png"
print "  - NFkB_full_ISEF.png"
print "  - NFkB_dim_ISEF.png"
'''

# Write the script
with open(r'C:\\Users\\vasud\\nod2-screening-data\\isef_figures\\pathway_diagram\\pathway_proteins.pml', 'w', encoding='utf-8') as f:
    f.write(script_content)

print("Created: pathway_proteins.pml")

# Also create a Python script to generate MDP 2D structure using RDKit
mdp_script = '''#!/usr/bin/env python3
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
        output_path = Path(r"C:\\Users\\vasud\\nod2-screening-data\\isef_figures\\pathway_diagram\\MDP_3D_ISEF.png")

        img = Draw.MolToImage(mol, size=(800, 800), kekulize=True)
        img.save(str(output_path))
        print(f"Created: {output_path}")
    else:
        print("Could not parse MDP SMILES")

except ImportError:
    print("RDKit not available. Install with: conda install -c conda-forge rdkit")
    print("Or use ChemDraw/MarvinSketch to create MDP structure manually.")
'''

with open(r'C:\\Users\\vasud\\nod2-screening-data\\isef_figures\\pathway_diagram\\generate_mdp_structure.py', 'w', encoding='utf-8') as f:
    f.write(mdp_script)

print("Created: generate_mdp_structure.py")
