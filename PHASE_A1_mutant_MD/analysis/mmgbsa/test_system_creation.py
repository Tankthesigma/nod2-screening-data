#!/usr/bin/env python
"""Test MM-GBSA system creation"""

import os
import sys
import warnings
warnings.filterwarnings('ignore')

# Check OpenMM
try:
    from openmm import *
    from openmm.app import *
    from openmm.unit import *
    print(f"OpenMM version: {version.version}")
except ImportError as e:
    print(f"OpenMM not available: {e}")
    sys.exit(1)

import MDAnalysis as mda

base = r'C:\Users\vasud\nod2-screening-data'
pdb_path = os.path.join(base, r'PHASE_6\trajectories\febuxostat_rep1_solvated.pdb')

print(f"\n{'='*60}")
print("TEST: System Creation for MM-GBSA")
print(f"{'='*60}")

# Step 1: Load with MDAnalysis and strip
print("\n1. Loading and stripping PDB...")
u = mda.Universe(pdb_path)
print(f"   Total atoms: {len(u.atoms)}")

# Select protein + ligand (UNK)
protein = u.select_atoms('protein')
ligand = u.select_atoms('resname UNK')
complex_atoms = u.select_atoms('protein or resname UNK')

print(f"   Protein atoms: {len(protein)}")
print(f"   Ligand atoms: {len(ligand)}")
print(f"   Complex atoms: {len(complex_atoms)}")

# Save stripped PDB
stripped_pdb = os.path.join(base, r'PHASE_A1_mutant_MD\analysis\mmgbsa\test_stripped.pdb')
os.makedirs(os.path.dirname(stripped_pdb), exist_ok=True)
complex_atoms.write(stripped_pdb)
print(f"   Saved stripped PDB: {stripped_pdb}")

# Step 2: Try to load with OpenMM
print("\n2. Loading stripped PDB with OpenMM...")
try:
    pdb = PDBFile(stripped_pdb)
    print(f"   Loaded: {pdb.topology.getNumAtoms()} atoms")

    # List residue names
    resnames = set()
    for residue in pdb.topology.residues():
        resnames.add(residue.name)
    print(f"   Residue types: {sorted(resnames)}")
except Exception as e:
    print(f"   ERROR: {e}")
    sys.exit(1)

# Step 3: Try to create system with GBSA
print("\n3. Creating GBSA system...")
try:
    # Use AMBER force field with implicit solvent
    forcefield = ForceField('amber14-all.xml', 'implicit/obc2.xml')

    system = forcefield.createSystem(
        pdb.topology,
        nonbondedMethod=NoCutoff,  # No cutoff for implicit solvent
        constraints=HBonds
    )

    print(f"   System created successfully!")
    print(f"   Number of forces: {system.getNumForces()}")

    # List forces
    for i in range(system.getNumForces()):
        force = system.getForce(i)
        print(f"   - {force.__class__.__name__}")

except Exception as e:
    print(f"   ERROR creating system: {e}")
    print("\n   This likely means the ligand (UNK) is not recognized by AMBER force field.")
    print("   We need to parameterize the ligand separately.")
    sys.exit(1)

# Step 4: Test energy calculation
print("\n4. Testing energy calculation...")
try:
    integrator = VerletIntegrator(0.001*picoseconds)
    platform = Platform.getPlatformByName('CPU')
    context = Context(system, integrator, platform)
    context.setPositions(pdb.positions)

    state = context.getState(getEnergy=True)
    energy_kJ = state.getPotentialEnergy().value_in_unit(kilojoules_per_mole)
    energy_kcal = energy_kJ / 4.184

    print(f"   Potential Energy: {energy_kcal:.2f} kcal/mol")

    if abs(energy_kcal) < 1e6:
        print("   Energy is reasonable - system is valid!")
    else:
        print("   WARNING: Energy is very large - possible parameterization issue")

except Exception as e:
    print(f"   ERROR: {e}")
    sys.exit(1)

print(f"\n{'='*60}")
print("TEST COMPLETE - System creation works!")
print(f"{'='*60}")
