#!/usr/bin/env python
"""
Build FEP systems for Natural Product CID_10120 using OpenFF toolkit.
Run this in WSL with the fep conda environment.
"""

import numpy as np
from pathlib import Path
import json

# OpenMM imports
from openmm import *
from openmm.app import *
from openmm import unit

# OpenFF imports
from openff.toolkit import Molecule, ForceField as OFFForceField
from openff.interchange import Interchange

# Paths (Windows paths mapped to WSL)
BASE_DIR = Path("/mnt/c/Users/vasud/nod2-screening-data/fep_pmx_natural")
PHASE6_DIR = Path("/mnt/c/Users/vasud/nod2-screening-data/PHASE_6")
LIGAND_SDF = PHASE6_DIR / "structures" / "natural_top_docked.sdf"
SOLVATED_PDB = PHASE6_DIR / "trajectories" / "natural_top_rep1_solvated.pdb"

# FEP Parameters
TEMPERATURE = 310.0 * unit.kelvin
TIMESTEP = 2.0 * unit.femtoseconds
SOFTCORE_ALPHA = 0.5

def load_ligand():
    """Load and parameterize ligand with OpenFF."""
    print("Loading ligand from SDF...")
    mols = Molecule.from_file(str(LIGAND_SDF))
    # If multiple molecules, take the first one
    if isinstance(mols, list):
        mol = mols[0]
        print(f"  Loaded {len(mols)} molecules, using first one")
    else:
        mol = mols
    print(f"  Ligand: {mol.name or 'CID_10120'}")
    print(f"  Atoms: {mol.n_atoms}")
    print(f"  Formula: {mol.hill_formula}")
    return mol

def create_ligand_system(mol):
    """Create OpenMM system for ligand using OpenFF Sage."""
    print("\nParameterizing ligand with OpenFF 2.1.0 (Sage)...")

    # Load OpenFF force field
    off_ff = OFFForceField("openff-2.1.0.offxml")

    # Create interchange and convert to OpenMM
    interchange = Interchange.from_smirnoff(off_ff, mol.to_topology())

    # Get OpenMM system
    system = interchange.to_openmm()
    topology = interchange.to_openmm_topology()
    positions = interchange.positions.to_openmm()

    print(f"  System created: {system.getNumParticles()} particles")

    return system, topology, positions

def load_solvated_system():
    """Load the solvated system from Phase 6 MD."""
    print("\nLoading solvated system...")
    pdb = PDBFile(str(SOLVATED_PDB))
    print(f"  Atoms: {pdb.topology.getNumAtoms()}")
    return pdb

def identify_ligand_atoms(topology):
    """Find ligand atom indices in the solvated system."""
    ligand_indices = []
    for res in topology.residues():
        if res.name == 'UNK':  # Ligand is labeled UNK
            for atom in res.atoms():
                ligand_indices.append(atom.index)
    print(f"  Found {len(ligand_indices)} ligand atoms")
    return ligand_indices

def create_alchemical_system(base_system, ligand_indices):
    """
    Create alchemical system with softcore potentials for FEP.

    This modifies the system to allow decoupling of the ligand.
    """
    print("\nCreating alchemical system...")

    # Clone the system
    system = XmlSerializer.deserialize(XmlSerializer.serialize(base_system))

    # Find nonbonded force
    nonbonded = None
    for force in system.getForces():
        if isinstance(force, NonbondedForce):
            nonbonded = force
            break

    if nonbonded is None:
        raise ValueError("No NonbondedForce found in system")

    # Create custom nonbonded force for alchemical interactions
    # Softcore LJ: 4*epsilon*lambda_sterics * x^12 * (x^12 - 1)
    # where x = 1/sqrt(alpha*(1-lambda_sterics) + (r/sigma)^6)

    energy_expression = """
    4*epsilon*lambda_sterics*x*(x-1.0);
    x = 1.0/sqrt(alpha*(1.0-lambda_sterics) + y);
    y = (r/sigma)^6;
    sigma = 0.5*(sigma1+sigma2);
    epsilon = sqrt(epsilon1*epsilon2);
    """

    alchemical_force = CustomNonbondedForce(energy_expression)
    alchemical_force.addGlobalParameter("lambda_sterics", 1.0)
    alchemical_force.addGlobalParameter("alpha", SOFTCORE_ALPHA)

    alchemical_force.addPerParticleParameter("sigma")
    alchemical_force.addPerParticleParameter("epsilon")

    # Add particles
    for i in range(system.getNumParticles()):
        charge, sigma, epsilon = nonbonded.getParticleParameters(i)
        alchemical_force.addParticle([sigma, epsilon])

    # Set interaction groups
    ligand_set = set(ligand_indices)
    solvent_set = set(range(system.getNumParticles())) - ligand_set

    alchemical_force.addInteractionGroup(ligand_set, solvent_set)

    # Set cutoff
    alchemical_force.setNonbondedMethod(CustomNonbondedForce.CutoffPeriodic)
    alchemical_force.setCutoffDistance(nonbonded.getCutoffDistance())

    system.addForce(alchemical_force)

    # Also need to handle electrostatics alchemically
    # Add custom electrostatics force
    elec_expression = """
    lambda_elec*ONE_4PI_EPS0*chargeprod/r;
    chargeprod = charge1*charge2;
    ONE_4PI_EPS0 = 138.935456;
    """

    alchemical_elec = CustomNonbondedForce(elec_expression)
    alchemical_elec.addGlobalParameter("lambda_elec", 1.0)
    alchemical_elec.addPerParticleParameter("charge")

    for i in range(system.getNumParticles()):
        charge, sigma, epsilon = nonbonded.getParticleParameters(i)
        alchemical_elec.addParticle([charge])

    alchemical_elec.addInteractionGroup(ligand_set, solvent_set)
    alchemical_elec.setNonbondedMethod(CustomNonbondedForce.CutoffPeriodic)
    alchemical_elec.setCutoffDistance(nonbonded.getCutoffDistance())

    system.addForce(alchemical_elec)

    # Zero out ligand-solvent interactions in original nonbonded
    for i in ligand_indices:
        for j in range(system.getNumParticles()):
            if j not in ligand_set:
                nonbonded.addException(i, j, 0.0, 1.0, 0.0, replace=True)

    print("  Alchemical system created")
    print(f"  Global parameters: lambda_sterics, lambda_elec, alpha")

    return system

def save_system(system, topology, positions, sys_name):
    """Save system files for FEP."""
    out_dir = BASE_DIR / sys_name
    out_dir.mkdir(parents=True, exist_ok=True)

    # Save system XML
    with open(out_dir / "alchemical_system.xml", 'w') as f:
        f.write(XmlSerializer.serialize(system))
    print(f"  Saved: {out_dir / 'alchemical_system.xml'}")

    # Save topology PDB
    with open(out_dir / "topology.pdb", 'w') as f:
        PDBFile.writeFile(topology, positions, f)
    print(f"  Saved: {out_dir / 'topology.pdb'}")

    # Save positions
    pos_array = np.array([[p.x, p.y, p.z] for p in positions.value_in_unit(unit.nanometers)])
    np.save(out_dir / "positions.npy", pos_array)
    print(f"  Saved: {out_dir / 'positions.npy'}")

def main():
    print("="*70)
    print("FEP SYSTEM BUILDER - Natural Product CID_10120")
    print("="*70)

    # Load ligand and create parameters
    mol = load_ligand()

    # Load solvated system
    solvated = load_solvated_system()

    # Find ligand atoms
    ligand_indices = identify_ligand_atoms(solvated.topology)

    # Create force field and system
    print("\nCreating system with Amber14 + OpenFF...")
    forcefield = ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')

    # Try to create system (this will fail for UNK residue)
    # We need a different approach - use the OpenFF system for ligand

    print("\n[INFO] Full system building requires combining:")
    print("  1. Amber14 parameters for protein + water")
    print("  2. OpenFF parameters for ligand")
    print("  3. This is complex - using simplified approach")

    # For now, let's just verify the ligand can be parameterized
    lig_system, lig_top, lig_pos = create_ligand_system(mol)

    # Save ligand system for reference
    lig_dir = BASE_DIR / "ligand_params"
    lig_dir.mkdir(parents=True, exist_ok=True)

    with open(lig_dir / "ligand_system.xml", 'w') as f:
        f.write(XmlSerializer.serialize(lig_system))
    print(f"\n  Saved ligand system: {lig_dir / 'ligand_system.xml'}")

    # Save ligand positions
    if hasattr(lig_pos, 'value_in_unit'):
        lig_pos_array = np.array(lig_pos.value_in_unit(unit.nanometers))
    else:
        lig_pos_array = np.array(lig_pos)
    np.save(lig_dir / "ligand_positions.npy", lig_pos_array)
    print(f"  Saved ligand positions: {lig_dir / 'ligand_positions.npy'}")

    print("\n" + "="*70)
    print("LIGAND PARAMETERIZATION COMPLETE")
    print("="*70)
    print("\nNext step: Combine with solvated protein system")
    print("This requires openmmforcefields.generators.GAFFTemplateGenerator")

if __name__ == "__main__":
    main()
