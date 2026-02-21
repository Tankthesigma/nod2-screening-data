#!/usr/bin/env python
"""Test ligand parameterization with OpenFF."""
import sys
sys.stdout.reconfigure(line_buffering=True)

print("Testing ligand parameterization...")
print()

# Test 1: Load SDF with RDKit
print("[1] Loading SDF with RDKit...")
try:
    from rdkit import Chem
    mol = Chem.SDMolSupplier("fep_study/inputs/febuxostat.sdf")[0]
    if mol is None:
        print("[FAIL] Could not load SDF")
        sys.exit(1)
    print(f"  [PASS] Loaded molecule: {mol.GetNumAtoms()} atoms, {mol.GetNumHeavyAtoms()} heavy")
    smiles = Chem.MolToSmiles(mol)
    print(f"  SMILES: {smiles[:80]}...")
except Exception as e:
    print(f"[FAIL] RDKit error: {e}")
    sys.exit(1)

# Test 2: Create OpenFF Molecule
print()
print("[2] Creating OpenFF Molecule from RDKit mol...")
try:
    from openff.toolkit.topology import Molecule
    off_mol = Molecule.from_rdkit(mol, allow_undefined_stereo=True)
    print(f"  [PASS] OpenFF Molecule created: {off_mol.n_atoms} atoms")
except Exception as e:
    print(f"[FAIL] OpenFF Molecule error: {e}")
    import traceback
    traceback.print_exc()
    sys.exit(1)

# Test 3: Generate conformer (if needed)
print()
print("[3] Checking conformers...")
try:
    n_conf = off_mol.n_conformers
    print(f"  Conformers: {n_conf}")
    if n_conf == 0:
        print("  [INFO] Generating conformer...")
        off_mol.generate_conformers(n_conformers=1)
        print(f"  [PASS] Generated {off_mol.n_conformers} conformer(s)")
except Exception as e:
    print(f"[FAIL] Conformer error: {e}")
    import traceback
    traceback.print_exc()
    sys.exit(1)

# Test 4: Parameterize with OpenFF
print()
print("[4] Parameterizing with OpenFF 2.1.0 (Sage)...")
try:
    from openff.toolkit.typing.engines.smirnoff import ForceField
    ff = ForceField("openff-2.1.0.offxml")
    print("  [PASS] Loaded OpenFF 2.1.0 force field")

    from openff.toolkit.topology import Topology
    topology = Topology.from_molecules([off_mol])
    print("  [PASS] Created topology")

    # This is the key step that often fails
    interchange = ff.create_interchange(topology)
    print("  [PASS] Created interchange (parameterization complete)")
    print(f"  Topology atoms: {interchange.topology.n_atoms}")
except Exception as e:
    print(f"[FAIL] Parameterization error: {e}")
    import traceback
    traceback.print_exc()
    sys.exit(1)

print()
print("="*60)
print("[SUCCESS] Ligand parameterization test PASSED")
print("="*60)
