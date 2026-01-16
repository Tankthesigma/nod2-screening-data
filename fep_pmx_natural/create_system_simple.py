#!/usr/bin/env python
"""Simple system creation script - just creates the basic OpenMM system."""
import numpy as np
from pathlib import Path
from openmm import *
from openmm.app import *
from openmm import unit
from openmmforcefields.generators import GAFFTemplateGenerator
from openff.toolkit import Molecule

BASE = Path("/mnt/c/Users/vasud/nod2-screening-data/fep_pmx_natural")
PDB = "/mnt/c/Users/vasud/nod2-screening-data/PHASE_6/trajectories/natural_top_rep1_solvated.pdb"
SDF = "/mnt/c/Users/vasud/nod2-screening-data/PHASE_6/structures/natural_top_docked.sdf"

print("="*60)
print("SYSTEM CREATION")
print("="*60)

# Load ligand
print("\n1. Loading ligand...")
mol = Molecule.from_file(SDF)[0]
print(f"   Ligand: {mol.n_atoms} atoms, {mol.hill_formula}")

# Load PDB
print("\n2. Loading solvated PDB...")
pdb = PDBFile(PDB)
print(f"   System: {pdb.topology.getNumAtoms()} atoms")

# Create GAFF generator
print("\n3. Creating GAFF generator...")
gaff = GAFFTemplateGenerator(molecules=mol, forcefield="gaff-2.11")

# Create forcefield
print("\n4. Setting up force field...")
ff = ForceField("amber14-all.xml", "amber14/tip3pfb.xml")
ff.registerTemplateGenerator(gaff.generator)

# Create system
print("\n5. Creating system (may take ~60s)...")
import time
t0 = time.time()
system = ff.createSystem(
    pdb.topology,
    nonbondedMethod=PME,
    nonbondedCutoff=1.0*unit.nanometers,
    constraints=HBonds,
    hydrogenMass=1.5*unit.amu
)
print(f"   System created in {time.time()-t0:.1f}s")
print(f"   Particles: {system.getNumParticles()}")

# Find ligand indices
print("\n6. Finding ligand atoms...")
lig_idx = []
for res in pdb.topology.residues():
    if res.name == 'UNK':
        for a in res.atoms():
            lig_idx.append(a.index)
print(f"   Ligand indices: {lig_idx[0]} to {lig_idx[-1]} ({len(lig_idx)} atoms)")

# Save basic system
print("\n7. Saving system...")
for sys_name in ['wt_complex', 'mut_complex']:
    out_dir = BASE / sys_name
    out_dir.mkdir(parents=True, exist_ok=True)

    # System XML
    with open(out_dir / "system.xml", 'w') as f:
        f.write(XmlSerializer.serialize(system))
    print(f"   Saved: {out_dir / 'system.xml'}")

    # Topology
    with open(out_dir / "topology.pdb", 'w') as f:
        PDBFile.writeFile(pdb.topology, pdb.positions, f)
    print(f"   Saved: {out_dir / 'topology.pdb'}")

    # Positions
    pos = np.array([[p.x, p.y, p.z] for p in pdb.positions.value_in_unit(unit.nanometers)])
    np.save(out_dir / "positions.npy", pos)
    print(f"   Saved: {out_dir / 'positions.npy'}")

    # Ligand indices
    np.save(out_dir / "ligand_indices.npy", np.array(lig_idx))
    print(f"   Saved: {out_dir / 'ligand_indices.npy'}")

print("\n" + "="*60)
print("DONE - Basic system created")
print("="*60)
