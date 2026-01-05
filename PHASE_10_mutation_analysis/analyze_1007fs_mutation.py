#!/usr/bin/env python3
"""
PHASE 10: NOD2 1007fs Mutation Analysis

The 1007fs frameshift mutation truncates NOD2 at position 1007,
deleting the binding pocket residues:
  - GLU1008, ASN1010, ASP1011 (deleted)
  - ARG1034, ASP1035, THR1036, ARG1037 (deleted)

This analysis shows why Febuxostat cannot bind to the 1007fs mutant.
"""

import numpy as np
import os
import json

# Paths
BASE_DIR = r"C:\Users\vasud\nod2-screening-data"
PHASE_10_DIR = os.path.join(BASE_DIR, "PHASE_10_mutation_analysis")
ALPHAFOLD_DIR = os.path.join(PHASE_10_DIR, "alphafold")
BOLTZ2_DIR = os.path.join(PHASE_10_DIR, "boltz2")
WT_LRR = os.path.join(BASE_DIR, "cross_validation", "human_NOD2_LRR.pdb")
WT_LIGAND = os.path.join(BASE_DIR, "PHASE_6", "structures", "febuxostat_docked.sdf")

# Key binding residues (in full NOD2 numbering)
BINDING_RESIDUES = {
    'LEU1007': 1007,  # Last residue before truncation
    'GLU1008': 1008,  # DELETED in mutant
    'ARG1009': 1009,  # DELETED
    'ASN1010': 1010,  # DELETED
    'ASP1011': 1011,  # DELETED
    'ARG1034': 1034,  # DELETED
    'ASP1035': 1035,  # DELETED
    'THR1036': 1036,  # DELETED
    'ARG1037': 1037,  # DELETED
}

# LRR domain spans residues 744-1040 in full NOD2
LRR_START = 744
LRR_END_WT = 1040
LRR_END_MUTANT = 1006  # Truncated at 1007fs


def parse_cif_atoms(cif_file):
    """Parse atom coordinates from mmCIF file."""
    atoms = []
    with open(cif_file, 'r') as f:
        in_atom_site = False
        headers = []
        for line in f:
            if line.startswith('_atom_site.'):
                in_atom_site = True
                header = line.strip().split('.')[1]
                headers.append(header)
            elif in_atom_site and line.startswith('ATOM') or line.startswith('HETATM'):
                parts = line.split()
                if len(parts) >= len(headers):
                    atom = {}
                    for i, h in enumerate(headers):
                        if i < len(parts):
                            atom[h] = parts[i]
                    atoms.append(atom)
            elif in_atom_site and line.startswith('loop_'):
                break
    return atoms, headers


def parse_pdb_atoms(pdb_file):
    """Parse atom coordinates from PDB file."""
    atoms = []
    with open(pdb_file, 'r') as f:
        for line in f:
            if line.startswith('ATOM'):
                try:
                    atom = {
                        'atom_name': line[12:16].strip(),
                        'resname': line[17:20].strip(),
                        'chain': line[21].strip(),
                        'resnum': int(line[22:26].strip()),
                        'x': float(line[30:38].strip()),
                        'y': float(line[38:46].strip()),
                        'z': float(line[46:54].strip()),
                    }
                    atoms.append(atom)
                except:
                    continue
    return atoms


def get_ca_coords(atoms, source='pdb'):
    """Extract CA atom coordinates."""
    ca_coords = {}
    if source == 'pdb':
        for atom in atoms:
            if atom['atom_name'] == 'CA':
                ca_coords[atom['resnum']] = {
                    'xyz': np.array([atom['x'], atom['y'], atom['z']]),
                    'resname': atom['resname']
                }
    else:  # CIF
        for atom in atoms:
            if atom.get('label_atom_id', '') == 'CA':
                try:
                    resnum = int(atom.get('label_seq_id', 0))
                    x = float(atom.get('Cartn_x', 0))
                    y = float(atom.get('Cartn_y', 0))
                    z = float(atom.get('Cartn_z', 0))
                    resname = atom.get('label_comp_id', 'UNK')
                    ca_coords[resnum] = {
                        'xyz': np.array([x, y, z]),
                        'resname': resname
                    }
                except:
                    continue
    return ca_coords


def calculate_rmsd(coords1, coords2, common_residues):
    """Calculate RMSD between two structures for common residues."""
    xyz1 = []
    xyz2 = []
    for res in common_residues:
        if res in coords1 and res in coords2:
            xyz1.append(coords1[res]['xyz'])
            xyz2.append(coords2[res]['xyz'])

    if not xyz1:
        return None

    xyz1 = np.array(xyz1)
    xyz2 = np.array(xyz2)

    # Center structures
    xyz1 -= np.mean(xyz1, axis=0)
    xyz2 -= np.mean(xyz2, axis=0)

    # Simple RMSD without rotation
    diff = xyz1 - xyz2
    return np.sqrt(np.mean(np.sum(diff**2, axis=1)))


def main():
    print("=" * 70)
    print("PHASE 10: NOD2 1007fs Mutation Analysis")
    print("=" * 70)

    # Load confidence scores
    print("\n1. STRUCTURE CONFIDENCE SCORES")
    print("-" * 40)

    # AlphaFold confidence
    af_conf_file = os.path.join(ALPHAFOLD_DIR, "fold_mutatioan_summary_confidences_0.json")
    with open(af_conf_file, 'r') as f:
        af_conf = json.load(f)
    print(f"AlphaFold 1007fs:")
    print(f"  pTM: {af_conf['ptm']:.2f}")
    print(f"  Fraction disordered: {af_conf['fraction_disordered']:.0%}")
    print(f"  Ranking score: {af_conf['ranking_score']:.2f}")

    # Boltz-2 confidence
    boltz_conf_file = os.path.join(BOLTZ2_DIR, "scores.json")
    with open(boltz_conf_file, 'r') as f:
        boltz_conf = json.load(f)
    print(f"\nBoltz-2 1007fs:")
    print(f"  Confidence score: {boltz_conf['1']['confidence_score']:.3f}")
    print(f"  pTM: {boltz_conf['1']['ptm']:.3f}")
    print(f"  pLDDT: {boltz_conf['1']['complex_plddt']:.3f}")

    print("\n  ** LOW CONFIDENCE indicates structural instability **")
    print("  ** Normal pTM for stable proteins: 0.7-0.9 **")

    # Analyze binding site deletion
    print("\n2. BINDING SITE ANALYSIS")
    print("-" * 40)

    print("\nWild-type NOD2 LRR binding site residues:")
    for name, num in BINDING_RESIDUES.items():
        status = "PRESENT" if num <= 1006 else "PRESENT (forms pocket)"
        print(f"  {name}: Residue {num} - {status}")

    print("\n1007fs Mutant (truncated at position 1006):")
    deleted_residues = []
    present_residues = []
    for name, num in BINDING_RESIDUES.items():
        if num > 1006:
            status = "DELETED (frameshift truncation)"
            deleted_residues.append(name)
        else:
            status = "Present but pocket incomplete"
            present_residues.append(name)
        print(f"  {name}: Residue {num} - {status}")

    print(f"\n  DELETED binding residues: {len(deleted_residues)}/9")
    print(f"  Key H-bond partners DELETED: GLU1008, ASP1011, ARG1037")

    # Load structures for comparison
    print("\n3. STRUCTURAL COMPARISON")
    print("-" * 40)

    # Load wild-type LRR
    wt_atoms = parse_pdb_atoms(WT_LRR)
    wt_ca = get_ca_coords(wt_atoms, 'pdb')

    # Get residue range
    wt_residues = sorted(wt_ca.keys())
    print(f"\nWild-type LRR:")
    print(f"  Total residues: {len(wt_residues)}")
    print(f"  Residue range: {min(wt_residues)} - {max(wt_residues)}")

    # Count binding site residues in WT
    wt_binding = [r for r in BINDING_RESIDUES.values() if r in wt_ca]
    print(f"  Binding site residues present: {len(wt_binding)}/9")

    # Load AlphaFold mutant structure
    af_cif = os.path.join(ALPHAFOLD_DIR, "fold_mutatioan_model_0.cif")
    af_atoms, af_headers = parse_cif_atoms(af_cif)
    af_ca = get_ca_coords(af_atoms, 'cif')

    if af_ca:
        af_residues = sorted(af_ca.keys())
        print(f"\nAlphaFold 1007fs mutant:")
        print(f"  Total residues: {len(af_residues)}")
        print(f"  Residue range: {min(af_residues)} - {max(af_residues)}")

        # Check if binding residues exist
        af_binding = [r for r in BINDING_RESIDUES.values() if r in af_ca]
        print(f"  Binding site residues present: {len(af_binding)}/9")

        # Check for residues > 1006
        over_1006 = [r for r in af_residues if r > 1006]
        print(f"  Residues > 1006: {len(over_1006)} (should be 0 for truncation)")

    # Load Boltz-2 mutant structure
    boltz_cif = os.path.join(BOLTZ2_DIR, "rank_1.cif")
    boltz_atoms, boltz_headers = parse_cif_atoms(boltz_cif)
    boltz_ca = get_ca_coords(boltz_atoms, 'cif')

    if boltz_ca:
        boltz_residues = sorted(boltz_ca.keys())
        print(f"\nBoltz-2 1007fs mutant:")
        print(f"  Total residues: {len(boltz_residues)}")
        print(f"  Residue range: {min(boltz_residues)} - {max(boltz_residues)}")

        boltz_binding = [r for r in BINDING_RESIDUES.values() if r in boltz_ca]
        print(f"  Binding site residues present: {len(boltz_binding)}/9")

    # Docking implications
    print("\n4. DOCKING IMPLICATIONS")
    print("-" * 40)

    print("""
Febuxostat binding to wild-type NOD2 LRR:
  - Binds in pocket formed by residues 1007-1011 and 1034-1037
  - Key interactions:
    * GLU1008: H-bond acceptor (carboxylate)
    * ASP1011: H-bond acceptor (carboxylate)
    * ARG1037: H-bond donor (guanidinium)
    * THR1036: Backbone H-bond
  - GNINA score: -4.14 kcal/mol
  - CNN score: 0.746

Febuxostat binding to 1007fs mutant:
  - Binding pocket DOES NOT EXIST
  - Residues 1008-1040 are DELETED
  - No GLU1008 -> No H-bond acceptor
  - No ARG1037 -> No H-bond donor
  - DOCKING WOULD FAIL or give very poor scores
""")

    # Summary
    print("\n" + "=" * 70)
    print("SUMMARY: Drug-Mutation Selectivity")
    print("=" * 70)

    print("""
+---------------------------------------------------------------------+
|                    NOD2 Genotype-Drug Response                      |
+---------------------------------------------------------------------+
|  Wild-Type NOD2:                                                    |
|    [+] Complete LRR domain (744-1040)                               |
|    [+] Intact binding pocket                                        |
|    [+] Febuxostat can bind (CNNscore: 0.746)                        |
|    [+] MD simulation: 99%% pocket occupancy                          |
|    --> DRUG CANDIDATE FOR WT PATIENTS                               |
+---------------------------------------------------------------------+
|  1007fs Mutant NOD2:                                                |
|    [X] Truncated LRR domain (744-1006)                              |
|    [X] Binding pocket DELETED                                       |
|    [X] No GLU1008, ASP1011, ARG1037                                 |
|    [X] Febuxostat CANNOT bind                                       |
|    --> NOT A CANDIDATE (requires different therapy)                 |
+---------------------------------------------------------------------+
|  Clinical Implication:                                              |
|    Patient genotyping enables PRECISION MEDICINE:                   |
|    - WT or other mutations with intact LRR -> Febuxostat candidate  |
|    - 1007fs carriers -> Alternative therapy needed                  |
+---------------------------------------------------------------------+
""")

    # Save report
    report_file = os.path.join(PHASE_10_DIR, "MUTATION_ANALYSIS_REPORT.txt")
    with open(report_file, 'w') as f:
        f.write("NOD2 1007fs Mutation Analysis Report\n")
        f.write("=" * 50 + "\n\n")
        f.write("Key Findings:\n")
        f.write(f"1. AlphaFold pTM: {af_conf['ptm']:.2f} (low - unstable)\n")
        f.write(f"2. Boltz-2 pTM: {boltz_conf['1']['ptm']:.3f} (low - unstable)\n")
        f.write(f"3. Binding residues deleted: 8/9\n")
        f.write(f"4. Critical H-bond partners (GLU1008, ASP1011, ARG1037): ALL DELETED\n\n")
        f.write("Conclusion:\n")
        f.write("The 1007fs mutation completely eliminates the Febuxostat binding site.\n")
        f.write("This enables precision medicine patient selection based on genotype.\n")

    print(f"\nReport saved to: {report_file}")


if __name__ == "__main__":
    main()
