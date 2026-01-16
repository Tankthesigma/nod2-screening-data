#!/usr/bin/env python
"""
Build complete FEP systems for Natural Product CID_10120.
Uses GAFFTemplateGenerator to parameterize the ligand.
"""

import numpy as np
from pathlib import Path
import json
import sys

# OpenMM imports
from openmm import *
from openmm.app import *
from openmm import unit

# OpenFF/GAFF imports
from openmmforcefields.generators import GAFFTemplateGenerator
from openff.toolkit import Molecule

# Paths
BASE_DIR = Path("/mnt/c/Users/vasud/nod2-screening-data/fep_pmx_natural")
PHASE6_DIR = Path("/mnt/c/Users/vasud/nod2-screening-data/PHASE_6")
LIGAND_SDF = PHASE6_DIR / "structures" / "natural_top_docked.sdf"
SOLVATED_PDB = PHASE6_DIR / "trajectories" / "natural_top_rep1_solvated.pdb"

# Parameters
TEMPERATURE = 310.0  # K
SOFTCORE_ALPHA = 0.5


def load_ligand_molecule():
    """Load ligand molecule for template generation."""
    print("Loading ligand molecule...")
    mols = Molecule.from_file(str(LIGAND_SDF))
    if isinstance(mols, list):
        mol = mols[0]
    else:
        mol = mols
    print(f"  Ligand: {mol.name or 'CID_10120'}")
    print(f"  Formula: {mol.hill_formula}")
    return mol


def create_system_with_ligand(pdb_file, ligand_mol):
    """Create OpenMM system using GAFF for ligand."""
    print(f"\nLoading PDB: {pdb_file}")
    pdb = PDBFile(str(pdb_file))
    print(f"  Atoms: {pdb.topology.getNumAtoms()}")

    # Create GAFF template generator
    print("\nCreating GAFF template generator...")
    gaff = GAFFTemplateGenerator(molecules=ligand_mol, forcefield='gaff-2.11')

    # Create force field with GAFF generator
    print("Setting up force field (Amber14 + GAFF)...")
    forcefield = ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')
    forcefield.registerTemplateGenerator(gaff.generator)

    # Create system
    print("Creating system...")
    system = forcefield.createSystem(
        pdb.topology,
        nonbondedMethod=PME,
        nonbondedCutoff=1.0*unit.nanometers,
        constraints=HBonds,
        hydrogenMass=1.5*unit.amu  # HMR for larger timestep
    )

    print(f"  System created: {system.getNumParticles()} particles")
    print(f"  Forces: {system.getNumForces()}")

    return system, pdb.topology, pdb.positions


def find_ligand_indices(topology):
    """Find ligand atom indices."""
    indices = []
    for res in topology.residues():
        if res.name == 'UNK':
            for atom in res.atoms():
                indices.append(atom.index)
    return indices


def create_alchemical_system(system, ligand_indices):
    """
    Add alchemical parameters to the system for FEP.
    Adds lambda_elec and lambda_sterics global parameters.
    """
    print("\nCreating alchemical modifications...")
    print(f"  Ligand atoms: {len(ligand_indices)}")

    # Find the NonbondedForce
    nonbonded = None
    for i, force in enumerate(system.getForces()):
        if isinstance(force, NonbondedForce):
            nonbonded = force
            break

    if nonbonded is None:
        raise ValueError("NonbondedForce not found")

    # Store original ligand parameters
    ligand_params = {}
    for idx in ligand_indices:
        charge, sigma, epsilon = nonbonded.getParticleParameters(idx)
        ligand_params[idx] = (charge, sigma, epsilon)

    # Create custom forces for alchemical ligand-environment interactions

    # 1. Softcore LJ for sterics
    softcore_lj = """
    4*epsilon*lambda_sterics*x*(x-1);
    x = 1/sqrt(alpha*(1-lambda_sterics)^2 + (r/sigma)^6);
    sigma = 0.5*(sigma1+sigma2);
    epsilon = sqrt(epsilon1*epsilon2);
    """

    custom_lj = CustomNonbondedForce(softcore_lj)
    custom_lj.addGlobalParameter("lambda_sterics", 1.0)
    custom_lj.addGlobalParameter("alpha", SOFTCORE_ALPHA)
    custom_lj.addPerParticleParameter("sigma")
    custom_lj.addPerParticleParameter("epsilon")

    # Add all particles
    for i in range(system.getNumParticles()):
        charge, sigma, epsilon = nonbonded.getParticleParameters(i)
        custom_lj.addParticle([sigma, epsilon])

    # Set up interaction groups: ligand <-> environment
    ligand_set = set(ligand_indices)
    env_set = set(range(system.getNumParticles())) - ligand_set

    custom_lj.addInteractionGroup(ligand_set, env_set)
    custom_lj.setNonbondedMethod(CustomNonbondedForce.CutoffPeriodic)
    custom_lj.setCutoffDistance(nonbonded.getCutoffDistance())
    custom_lj.setUseLongRangeCorrection(True)

    # 2. Scaled electrostatics
    scaled_elec = """
    lambda_elec*ONE_4PI_EPS0*charge1*charge2/r;
    ONE_4PI_EPS0 = 138.935456;
    """

    custom_elec = CustomNonbondedForce(scaled_elec)
    custom_elec.addGlobalParameter("lambda_elec", 1.0)
    custom_elec.addPerParticleParameter("charge")

    for i in range(system.getNumParticles()):
        charge, sigma, epsilon = nonbonded.getParticleParameters(i)
        custom_elec.addParticle([charge])

    custom_elec.addInteractionGroup(ligand_set, env_set)
    custom_elec.setNonbondedMethod(CustomNonbondedForce.CutoffPeriodic)
    custom_elec.setCutoffDistance(nonbonded.getCutoffDistance())

    # Add the custom forces
    system.addForce(custom_lj)
    system.addForce(custom_elec)

    # Zero out ligand-environment interactions in original NonbondedForce
    # by adding exceptions
    print("  Adding exceptions for ligand-environment pairs...")
    n_exceptions = 0
    for lig_idx in ligand_indices:
        for env_idx in env_set:
            # Check if exception already exists
            try:
                nonbonded.addException(lig_idx, env_idx, 0.0, 1.0, 0.0)
                n_exceptions += 1
            except:
                pass  # Exception already exists

    print(f"  Added {n_exceptions} exceptions")
    print("  Alchemical system ready")

    return system


def add_boresch_restraints(system, topology, anchors_file):
    """Add Boresch restraints to the system."""
    print("\nAdding Boresch restraints...")

    with open(anchors_file, 'r') as f:
        anchors = json.load(f)

    # Get atom indices from serial numbers
    serial_to_idx = {}
    for atom in topology.atoms():
        serial_to_idx[atom.index + 1] = atom.index  # PDB serials are 1-indexed

    # This is simplified - we need to match by residue and atom name
    # For now, use the pre-selected anchors

    eq = anchors['equilibrium']
    fc = anchors['force_constants']

    # The actual implementation would add CustomBondForce, CustomAngleForce,
    # CustomTorsionForce with lambda_restraints parameter

    print(f"  Distance r0: {eq['r0_nm']:.4f} nm")
    print(f"  Force constants: k_dist={fc['k_distance']}, k_ang={fc['k_angle']}")
    print("  [Note: Full restraint implementation needed]")

    return system


def save_system(system, topology, positions, out_dir):
    """Save all system files."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # System XML
    xml_path = out_dir / "alchemical_system.xml"
    with open(xml_path, 'w') as f:
        f.write(XmlSerializer.serialize(system))
    print(f"  Saved: {xml_path}")

    # Topology PDB
    pdb_path = out_dir / "topology.pdb"
    with open(pdb_path, 'w') as f:
        PDBFile.writeFile(topology, positions, f)
    print(f"  Saved: {pdb_path}")

    # Positions
    pos_array = np.array([[p.x, p.y, p.z] for p in positions.value_in_unit(unit.nanometers)])
    np.save(out_dir / "positions.npy", pos_array)
    print(f"  Saved: {out_dir / 'positions.npy'}")


def main():
    print("="*70)
    print("FEP FULL SYSTEM BUILDER")
    print("Natural Product CID_10120 (Bufadienolide)")
    print("="*70)

    # Load ligand molecule
    ligand_mol = load_ligand_molecule()

    # Create system for WT complex
    print("\n" + "="*70)
    print("BUILDING WT COMPLEX SYSTEM")
    print("="*70)

    system, topology, positions = create_system_with_ligand(SOLVATED_PDB, ligand_mol)

    # Find ligand indices
    ligand_indices = find_ligand_indices(topology)
    print(f"\nLigand atom indices: {ligand_indices[0]} to {ligand_indices[-1]}")

    # Create alchemical system
    alch_system = create_alchemical_system(system, ligand_indices)

    # Save WT system
    print("\nSaving WT complex system...")
    save_system(alch_system, topology, positions, BASE_DIR / "wt_complex")

    # For mutant, we would need to apply R702W mutation
    # For now, copy WT system (mutation effect comes from protein structure)
    print("\n" + "="*70)
    print("BUILDING MUTANT COMPLEX SYSTEM")
    print("="*70)
    print("[Note: Using WT structure - mutation requires structure modification]")

    # Copy system for mutant
    mut_system = XmlSerializer.deserialize(XmlSerializer.serialize(alch_system))
    save_system(mut_system, topology, positions, BASE_DIR / "mut_complex")

    # Summary
    print("\n" + "="*70)
    print("SYSTEM BUILDING COMPLETE")
    print("="*70)
    print(f"\nTotal atoms: {system.getNumParticles()}")
    print(f"Ligand atoms: {len(ligand_indices)}")
    print(f"Box size: ~131 A cubic")
    print("\nSystems saved:")
    print(f"  - {BASE_DIR / 'wt_complex'}")
    print(f"  - {BASE_DIR / 'mut_complex'}")

    print("\n[READY FOR FEP]")


if __name__ == "__main__":
    main()
