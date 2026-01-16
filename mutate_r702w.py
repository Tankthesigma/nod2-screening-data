#!/usr/bin/env python
"""
Mutate ARG 702 to TRP using PDBFixer.

This script performs the R702W mutation on the protein structure.
"""
import os
import sys

print("="*60)
print("R702W MUTATION SCRIPT")
print("="*60)

try:
    from pdbfixer import PDBFixer
    from openmm.app import PDBFile

    print("[INFO] Using PDBFixer for mutation...")

    fixer = PDBFixer(filename="fep_pmx/wt/protein.pdb")

    # Find target chain and residue
    target_chain_id = None
    target_found = False

    for chain in fixer.topology.chains():
        for residue in chain.residues():
            try:
                res_id = int(str(residue.id).strip())
            except (ValueError, AttributeError):
                continue

            if res_id == 702 and residue.name in ['ARG', 'ARN']:
                print(f"  Found {residue.name} 702 in chain '{chain.id}'")
                target_chain_id = chain.id
                target_found = True
                break
        if target_found:
            break

    if not target_found:
        print("[FAIL] Could not find ARG 702!")
        sys.exit(1)

    # PDBFixer mutation format: list of "OLDRES-RESID-NEWRES" strings, chain_id separate
    mutation_str = "ARG-702-TRP"
    print(f"  Applying mutation: {mutation_str} on chain {target_chain_id}")
    fixer.applyMutations([mutation_str], target_chain_id)

    # Find and add missing atoms for the new residue
    fixer.findMissingResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(7.0)

    # Save mutant
    with open("fep_pmx/mutant/protein.pdb", "w") as f:
        PDBFile.writeFile(fixer.topology, fixer.positions, f)
    print("[PASS] Saved mutant structure to fep_pmx/mutant/protein.pdb")

    # Verify mutation worked
    print()
    print("[INFO] Verifying mutation...")
    fixer2 = PDBFixer(filename="fep_pmx/mutant/protein.pdb")
    verified = False
    for chain in fixer2.topology.chains():
        for residue in chain.residues():
            try:
                res_id = int(str(residue.id).strip())
            except (ValueError, AttributeError):
                continue
            if res_id == 702:
                print(f"  Residue 702 is now: {residue.name}")
                if residue.name == 'TRP':
                    print("[PASS] R702W mutation successful!")
                    verified = True
                else:
                    print(f"[FAIL] Expected TRP, got {residue.name}")
                break
        if verified:
            break

    if not verified:
        print("[FAIL] Could not verify mutation!")
        sys.exit(1)

except ImportError as e:
    print(f"[FAIL] PDBFixer not available: {e}")
    print("  Install with: pip install pdbfixer")
    sys.exit(1)

except Exception as e:
    print(f"[FAIL] Error: {e}")
    import traceback
    traceback.print_exc()
    sys.exit(1)

print()
print("="*60)
print("Mutation complete! Now re-run setup_fep_pmx.py to continue.")
print("="*60)
