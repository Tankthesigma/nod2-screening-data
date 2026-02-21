#!/usr/bin/env python
"""Test perses with RDKit patch."""
import sys
import os

sys.stdout.reconfigure(line_buffering=True)
os.environ["PYTHONUNBUFFERED"] = "1"

print("="*60)
print("PATCHED PERSES TEST")
print("="*60)
print()

# CRITICAL: Apply patch BEFORE importing perses
print("[1] Applying perses patch...")
import patch_perses
print()

# Now import perses
print("[2] Importing perses...")
try:
    from perses.app.relative_point_mutation_setup import PointMutationExecutor
    import openmm.unit as unit
    print("  [PASS] Imported PointMutationExecutor")
except Exception as e:
    print(f"  [FAIL] Import error: {e}")
    import traceback
    traceback.print_exc()
    sys.exit(1)

print()
print("[3] Testing protein-only mutation setup...")
print("  Protein: fep_study/inputs/WT_protein_for_complex.pdb")
print("  Mutation: ARG 702 -> TRP (chain A)")
print()

kwargs = {
    "protein_filename": "fep_study/inputs/WT_protein_for_complex.pdb",
    "mutation_chain_id": "A",
    "mutation_residue_id": "702",
    "proposed_residue": "TRP",
    "conduct_endstate_validation": False,
    "ionic_strength": 0.15 * unit.molar,
}

print("  Parameters:")
for k, v in kwargs.items():
    print(f"    {k}: {v}")
print()

print("  [INFO] Calling PointMutationExecutor...")
print("  [INFO] This may take 2-5 minutes...")
print()

try:
    executor = PointMutationExecutor(**kwargs)
    print("  [PASS] PointMutationExecutor created successfully!")
    print(f"  Hybrid system: {executor.hybrid_system.getNumParticles()} particles")
except Exception as e:
    print(f"  [FAIL] Exception: {type(e).__name__}")
    print(f"  Message: {e}")
    print()
    print("  Full traceback:")
    import traceback
    traceback.print_exc()
    sys.exit(1)

print()
print("="*60)
print("[SUCCESS] Patched perses test PASSED!")
print("="*60)
