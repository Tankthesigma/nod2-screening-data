#!/usr/bin/env python
"""Direct mutation setup using OpenMMTools (without perses PointMutationExecutor)."""
import sys
import os

sys.stdout.reconfigure(line_buffering=True)
os.environ["PYTHONUNBUFFERED"] = "1"

print("="*60)
print("DIRECT MUTATION SETUP TEST")
print("="*60)
print()

# First, check if we can avoid OpenEye entirely
print("[1] Checking for OpenEye import...")
try:
    # Try to prevent OpenEye from loading
    sys.modules['openeye'] = None
    sys.modules['openeye.oechem'] = None
    print("  [INFO] Blocked OpenEye modules")
except:
    pass

print()
print("[2] Importing OpenMM and OpenFF...")
try:
    from openmm import app, unit
    from openmm import System, LangevinIntegrator
    from openff.toolkit import Molecule, Topology, ForceField
    print("  [PASS] Basic imports OK")
except Exception as e:
    print(f"  [FAIL] Import error: {e}")
    sys.exit(1)

print()
print("[3] Loading protein PDB...")
try:
    pdb = app.PDBFile("fep_study/inputs/WT_protein_for_complex.pdb")
    print(f"  [PASS] Loaded PDB: {pdb.topology.getNumAtoms()} atoms")
except Exception as e:
    print(f"  [FAIL] PDB load error: {e}")
    sys.exit(1)

print()
print("[4] Finding ARG 702...")
try:
    target_res = None
    for res in pdb.topology.residues():
        if res.id == '702' and res.name in ['ARG', 'ARN']:
            target_res = res
            break

    if target_res is None:
        print("  [FAIL] Could not find ARG 702")
        sys.exit(1)

    print(f"  [PASS] Found {target_res.name} 702 in chain {target_res.chain.id}")
    print(f"  Atoms in residue: {sum(1 for _ in target_res.atoms())}")
except Exception as e:
    print(f"  [FAIL] Residue search error: {e}")
    sys.exit(1)

print()
print("[5] Attempting direct OpenMMTools import...")
try:
    from openmmtools import alchemy
    from openmmtools.alchemy import AlchemicalState, AlchemicalRegion
    print("  [PASS] OpenMMTools alchemy imported")
except Exception as e:
    print(f"  [FAIL] OpenMMTools import error: {e}")
    import traceback
    traceback.print_exc()
    sys.exit(1)

print()
print("[6] Loading force field for protein...")
try:
    forcefield = app.ForceField('amber14/protein.ff14SB.xml', 'amber14/tip3p.xml')
    print("  [PASS] Loaded amber14/protein.ff14SB.xml")
except Exception as e:
    print(f"  [FAIL] Forcefield error: {e}")
    sys.exit(1)

print()
print("[7] Creating basic system (without solvation)...")
try:
    # Create system without water (just protein)
    system = forcefield.createSystem(
        pdb.topology,
        nonbondedMethod=app.NoCutoff,  # No cutoff for unsolvated
        constraints=app.HBonds
    )
    print(f"  [PASS] System created: {system.getNumParticles()} particles")
except Exception as e:
    print(f"  [FAIL] System creation error: {e}")
    import traceback
    traceback.print_exc()
    sys.exit(1)

print()
print("="*60)
print("[SUCCESS] Direct system creation works!")
print()
print("The issue is that perses requires OpenEye for some of its")
print("functionality. To proceed with FEP, you need either:")
print()
print("1. Get a free OpenEye academic license:")
print("   https://www.eyesopen.com/academic-licensing")
print()
print("2. Use a different FEP tool that doesn't require OpenEye")
print()
print("3. Set up the hybrid topology manually using OpenMMTools")
print("="*60)
