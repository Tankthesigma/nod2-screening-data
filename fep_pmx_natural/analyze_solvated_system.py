#!/usr/bin/env python
"""Analyze the solvated system from Phase 6 MD."""

from openmm import *
from openmm.app import *

pdb_file = 'C:/Users/vasud/nod2-screening-data/PHASE_6/trajectories/natural_top_rep1_solvated.pdb'

# Load topology from PDB
print("Loading solvated PDB...")
pdb = PDBFile(pdb_file)
print(f'Topology: {pdb.topology.getNumAtoms()} atoms, {pdb.topology.getNumResidues()} residues')

# Count chains
chains = list(pdb.topology.chains())
print(f'Chains: {len(chains)}')
for c in chains:
    residues = list(c.residues())
    print(f'  Chain {c.id}: {len(residues)} residues, first={residues[0].name}, last={residues[-1].name}')

# Find ligand
print("\nSearching for ligand (UNK)...")
for res in pdb.topology.residues():
    if res.name == 'UNK':
        atoms = list(res.atoms())
        print(f"  Found UNK at chain {res.chain.id}, residue {res.index}")
        print(f"  Atoms: {len(atoms)}")
        for a in atoms[:5]:
            print(f"    {a.name} ({a.element.symbol})")
        print(f"    ... and {len(atoms)-5} more")
        break

# Count waters
n_water = sum(1 for r in pdb.topology.residues() if r.name == 'HOH')
print(f"\nWater molecules: {n_water}")

# Count ions
n_cl = sum(1 for r in pdb.topology.residues() if r.name == 'CL')
n_na = sum(1 for r in pdb.topology.residues() if r.name == 'NA')
print(f"Chloride ions: {n_cl}")
print(f"Sodium ions: {n_na}")

print("\n[DONE]")
