#!/usr/bin/env python
"""Deep dive into WT vs Mut system differences"""
import numpy as np
import os

BASE = "C:/Users/vasud/nod2-screening-data/fep_complete/fep_pmx"

print("="*70)
print("DEEP DIVE: WT vs MUT SYSTEM DIFFERENCES")
print("="*70)

# 1. ATOM COUNTS
print("\n1. ATOM COUNTS")
print("-"*70)
for sys in ['wt_complex', 'mut_complex', 'solvent']:
    pos_path = f"{BASE}/{sys}/positions.npy"
    topo_path = f"{BASE}/{sys}/topology.pdb"

    if os.path.exists(pos_path):
        pos = np.load(pos_path)
        print(f"{sys}: {pos.shape[0]} atoms (from positions.npy)")

    if os.path.exists(topo_path):
        with open(topo_path) as f:
            atom_count = sum(1 for line in f if line.startswith('ATOM') or line.startswith('HETATM'))
        print(f"{sys}: {atom_count} atoms (from topology.pdb)")

# 2. RESIDUE 702
print("\n2. RESIDUE 702")
print("-"*70)
for sys in ['wt_complex', 'mut_complex']:
    topo_path = f"{BASE}/{sys}/topology.pdb"
    if os.path.exists(topo_path):
        print(f"\n{sys}:")
        with open(topo_path) as f:
            res702_atoms = []
            for line in f:
                if (line.startswith('ATOM') or line.startswith('HETATM')):
                    try:
                        res_num = int(line[22:26].strip())
                        if res_num == 702:
                            res_name = line[17:20].strip()
                            atom_name = line[12:16].strip()
                            res702_atoms.append((res_name, atom_name))
                    except:
                        pass
            if res702_atoms:
                res_name = res702_atoms[0][0]
                print(f"  Residue 702 = {res_name} ({len(res702_atoms)} atoms)")
                print(f"  Atoms: {[a[1] for a in res702_atoms[:10]]}...")
            else:
                print("  Residue 702 NOT FOUND")

# 3. BOX DIMENSIONS
print("\n3. BOX DIMENSIONS")
print("-"*70)
for sys in ['wt_complex', 'mut_complex', 'solvent']:
    topo_path = f"{BASE}/{sys}/topology.pdb"
    if os.path.exists(topo_path):
        with open(topo_path) as f:
            for line in f:
                if line.startswith('CRYST1'):
                    a = float(line[6:15])
                    b = float(line[15:24])
                    c = float(line[24:33])
                    print(f"{sys}: {a:.2f} x {b:.2f} x {c:.2f} Angstrom")
                    break

# 4. LIGAND POSE
print("\n4. LIGAND POSE")
print("-"*70)
for sys in ['wt_complex', 'mut_complex', 'solvent']:
    topo_path = f"{BASE}/{sys}/topology.pdb"
    pos_path = f"{BASE}/{sys}/positions.npy"

    if os.path.exists(topo_path) and os.path.exists(pos_path):
        positions = np.load(pos_path)

        # Find ligand atoms (UNL or LIG residue)
        ligand_indices = []
        res702_indices = []

        with open(topo_path) as f:
            for i, line in enumerate(f):
                if line.startswith('ATOM') or line.startswith('HETATM'):
                    res_name = line[17:20].strip()
                    if res_name in ['UNL', 'LIG', 'MOL']:
                        ligand_indices.append(i)
                    try:
                        res_num = int(line[22:26].strip())
                        if res_num == 702:
                            res702_indices.append(i)
                    except:
                        pass

        if ligand_indices:
            # This won't work directly - need to map PDB line to atom index
            # Let me try a different approach
            pass

print("\n  (Need MDAnalysis for detailed pose analysis)")

# 5. PROTONATION / CHARGE
print("\n5. SYSTEM CHARGE")
print("-"*70)
# Check alchemical system XML for total charge
for sys in ['wt_complex', 'mut_complex', 'solvent']:
    xml_path = f"{BASE}/{sys}/alchemical_system.xml"
    if os.path.exists(xml_path):
        with open(xml_path) as f:
            content = f.read()
            # Count NonbondedForce particles with charge
            # This is complex - skip for now
        print(f"{sys}: Check XML manually or use OpenMM")

# 6. SOLVENT SYSTEM
print("\n6. SOLVENT SYSTEM")
print("-"*70)
solv_pos = f"{BASE}/solvent/positions.npy"
wt_pos = f"{BASE}/wt_complex/positions.npy"
if os.path.exists(solv_pos) and os.path.exists(wt_pos):
    solv = np.load(solv_pos)
    wt = np.load(wt_pos)
    print(f"Solvent system: {solv.shape[0]} atoms")
    print(f"WT complex: {wt.shape[0]} atoms")
    print(f"Difference: {wt.shape[0] - solv.shape[0]} atoms (should be protein)")

# 7. RAW dG INTERPRETATION
print("\n7. RAW ΔG_decouple INTERPRETATION")
print("-"*70)
print("ΔG_decouple values:")
print("  Solvent:     51.07 kcal/mol (cost to remove ligand from water)")
print("  WT_complex:  55.07 kcal/mol (cost to remove ligand from WT protein)")
print("  Mut_complex: 54.07 kcal/mol (cost to remove ligand from Mut protein)")
print()
print("Interpretation:")
print("  WT costs 4.00 kcal/mol MORE than solvent → ligand prefers WT")
print("  Mut costs 3.00 kcal/mol MORE than solvent → ligand prefers Mut (less)")
print("  WT - Mut = 1.00 kcal/mol → WT binds 1 kcal/mol tighter")
print()
print("Physical sense for ARG→TRP:")
print("  ARG702: positive charge, can form salt bridges/H-bonds")
print("  TRP702: neutral aromatic, pi-stacking but no charge")
print("  If febuxostat has negative groups (carboxylate), losing ARG")
print("  would weaken binding → CONSISTENT with +1 kcal/mol ΔΔG")
