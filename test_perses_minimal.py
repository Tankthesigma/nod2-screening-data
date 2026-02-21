#!/usr/bin/env python
"""Minimal test of PointMutationExecutor."""
import sys
import os
import traceback as tb

sys.stdout.reconfigure(line_buffering=True)
sys.stderr.reconfigure(line_buffering=True)
os.environ["PYTHONUNBUFFERED"] = "1"

print("="*60)
print("MINIMAL PERSES TEST")
print("="*60)
print()

# Step 1: Import perses
print("[1] Importing perses...")
try:
    from perses.app.relative_point_mutation_setup import PointMutationExecutor
    import openmm.unit as unit
    print("  [PASS] Imported PointMutationExecutor")
except ImportError as e:
    print(f"  [FAIL] Import error: {e}")
    sys.exit(1)

# Step 2: Test with protein-only (no ligand) first
print()
print("[2] Testing protein-only mutation setup...")
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

print("  [INFO] Calling PointMutationExecutor (protein-only)...")
print("  [INFO] This may take 2-5 minutes...")
print()

import signal

def handler(signum, frame):
    print(f"\n  [SIGNAL] Received signal {signum}")
    sys.exit(128 + signum)

for sig in [signal.SIGTERM, signal.SIGINT, signal.SIGABRT]:
    try:
        signal.signal(sig, handler)
    except:
        pass

try:
    executor = PointMutationExecutor(**kwargs)
    print("  [PASS] PointMutationExecutor created successfully!")
    print(f"  Hybrid system: {executor.hybrid_system.getNumParticles()} particles")
except Exception as e:
    print(f"  [FAIL] Exception: {type(e).__name__}")
    print(f"  Message: {e}")
    print()
    print("  Full traceback:")
    tb.print_exc()
    sys.exit(1)

print()
print("="*60)
print("[SUCCESS] Protein-only mutation test PASSED")
print("="*60)
