#!/usr/bin/env python
"""Add alchemical modifications to the base system - Windows version."""
import numpy as np
from pathlib import Path
from openmm import *
from openmm.app import *
from openmm import unit

BASE = Path("C:/Users/vasud/nod2-screening-data/fep_pmx_natural")
SOFTCORE_ALPHA = 0.5

def add_alchemical_forces(system, ligand_indices):
    """Add alchemical forces to the system."""
    print("Adding alchemical modifications...")

    # Find NonbondedForce
    nonbonded = None
    for f in system.getForces():
        if isinstance(f, NonbondedForce):
            nonbonded = f
            break

    if nonbonded is None:
        raise ValueError("NonbondedForce not found")

    n_particles = system.getNumParticles()
    ligand_set = set(ligand_indices)
    env_set = set(range(n_particles)) - ligand_set

    print(f"  Ligand atoms: {len(ligand_set)}")
    print(f"  Environment atoms: {len(env_set)}")

    # 1. Softcore LJ
    print("  Creating softcore LJ force...")
    sc_lj = """
    4*epsilon*lambda_sterics*x*(x-1);
    x = 1/sqrt(alpha*(1-lambda_sterics)^2 + (r/sigma)^6);
    sigma = 0.5*(sigma1+sigma2);
    epsilon = sqrt(epsilon1*epsilon2);
    """
    custom_lj = CustomNonbondedForce(sc_lj)
    custom_lj.addGlobalParameter("lambda_sterics", 1.0)
    custom_lj.addGlobalParameter("alpha", SOFTCORE_ALPHA)
    custom_lj.addPerParticleParameter("sigma")
    custom_lj.addPerParticleParameter("epsilon")

    for i in range(n_particles):
        q, sigma, eps = nonbonded.getParticleParameters(i)
        custom_lj.addParticle([sigma, eps])

    custom_lj.addInteractionGroup(ligand_set, env_set)
    custom_lj.setNonbondedMethod(CustomNonbondedForce.CutoffPeriodic)
    custom_lj.setCutoffDistance(nonbonded.getCutoffDistance())
    custom_lj.setUseLongRangeCorrection(True)

    # Copy exclusions from NonbondedForce to CustomNonbondedForce
    print("  Copying exclusions to softcore LJ...")
    for i in range(nonbonded.getNumExceptions()):
        p1, p2, qq, sig, eps = nonbonded.getExceptionParameters(i)
        custom_lj.addExclusion(p1, p2)

    system.addForce(custom_lj)

    # 2. Scaled electrostatics
    print("  Creating scaled electrostatics force...")
    sc_elec = """
    lambda_elec*ONE_4PI_EPS0*charge1*charge2/r;
    ONE_4PI_EPS0 = 138.935456;
    """
    custom_elec = CustomNonbondedForce(sc_elec)
    custom_elec.addGlobalParameter("lambda_elec", 1.0)
    custom_elec.addPerParticleParameter("charge")

    for i in range(n_particles):
        q, sigma, eps = nonbonded.getParticleParameters(i)
        custom_elec.addParticle([q])

    custom_elec.addInteractionGroup(ligand_set, env_set)
    custom_elec.setNonbondedMethod(CustomNonbondedForce.CutoffPeriodic)
    custom_elec.setCutoffDistance(nonbonded.getCutoffDistance())

    # Copy exclusions from NonbondedForce to custom electrostatics
    print("  Copying exclusions to scaled electrostatics...")
    for i in range(nonbonded.getNumExceptions()):
        p1, p2, qq, sig, eps = nonbonded.getExceptionParameters(i)
        custom_elec.addExclusion(p1, p2)

    system.addForce(custom_elec)

    # 3. Zero out ligand-environment interactions in original force
    # AND add them as exclusions to custom forces
    print("  Zeroing ligand-environment interactions...")

    # Get existing exceptions
    existing_exceptions = set()
    for i in range(nonbonded.getNumExceptions()):
        p1, p2, qq, sig, eps = nonbonded.getExceptionParameters(i)
        existing_exceptions.add((min(p1,p2), max(p1,p2)))

    # Add exceptions for ligand-environment pairs to NonbondedForce
    # AND add as exclusions to custom forces
    count = 0
    for lig_idx in ligand_indices:
        for env_idx in env_set:
            pair = (min(lig_idx, env_idx), max(lig_idx, env_idx))
            if pair not in existing_exceptions:
                nonbonded.addException(lig_idx, env_idx, 0.0, 1.0, 0.0)
                custom_lj.addExclusion(lig_idx, env_idx)
                custom_elec.addExclusion(lig_idx, env_idx)
                count += 1

    print(f"  Added {count} exceptions/exclusions")
    print("  Alchemical modifications complete")

    return system


def main():
    print("="*60)
    print("ADDING ALCHEMICAL MODIFICATIONS")
    print("="*60)

    for sys_name in ['wt_complex', 'mut_complex']:
        print(f"\nProcessing {sys_name}...")
        sys_dir = BASE / sys_name

        # Load system
        print("  Loading system.xml...")
        with open(sys_dir / "system.xml", 'r') as f:
            system = XmlSerializer.deserialize(f.read())
        print(f"  Particles: {system.getNumParticles()}")

        # Load ligand indices
        lig_idx = np.load(sys_dir / "ligand_indices.npy").tolist()
        print(f"  Ligand indices: {len(lig_idx)}")

        # Add alchemical forces
        alch_system = add_alchemical_forces(system, lig_idx)

        # Save alchemical system
        print("  Saving alchemical_system.xml...")
        with open(sys_dir / "alchemical_system.xml", 'w') as f:
            f.write(XmlSerializer.serialize(alch_system))
        print(f"  Saved: {sys_dir / 'alchemical_system.xml'}")

    print("\n" + "="*60)
    print("ALCHEMICAL SYSTEMS READY")
    print("="*60)


if __name__ == "__main__":
    main()
