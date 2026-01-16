#!/usr/bin/env python
"""Build alchemical system using openmmtools for proper handling."""
import numpy as np
from pathlib import Path
from openmm import *
from openmm.app import *
from openmm import unit

try:
    from openmmtools.alchemy import AbsoluteAlchemicalFactory, AlchemicalRegion, AlchemicalState
    HAS_OPENMMTOOLS = True
except ImportError:
    HAS_OPENMMTOOLS = False
    print("openmmtools not available, using manual approach")

BASE = Path("C:/Users/vasud/nod2-screening-data/fep_pmx_natural")


def add_alchemical_manual(system, ligand_indices):
    """Manual alchemical setup without exclusion conflicts."""
    print("Adding alchemical modifications (manual)...")

    # Find NonbondedForce
    nonbonded = None
    nb_index = None
    for i, f in enumerate(system.getForces()):
        if isinstance(f, NonbondedForce):
            nonbonded = f
            nb_index = i
            break

    if nonbonded is None:
        raise ValueError("NonbondedForce not found")

    n_particles = system.getNumParticles()
    ligand_set = set(ligand_indices)

    print(f"  Ligand atoms: {len(ligand_set)}")
    print(f"  Total atoms: {n_particles}")

    # Store original ligand parameters
    original_params = {}
    for i in ligand_indices:
        q, sigma, eps = nonbonded.getParticleParameters(i)
        original_params[i] = (q, sigma, eps)

    # Create alchemical LJ force (ligand-env only, uses context parameters)
    print("  Creating alchemical LJ force...")
    # Standard softcore LJ
    alch_lj = CustomNonbondedForce(
        "4*epsilon*lambda_sterics*x*(x-1);"
        "x = 1/softcore_term;"
        "softcore_term = sqrt(softcore_alpha*(1-lambda_sterics)^2 + (r/sigma)^6);"
        "sigma = 0.5*(sigma1+sigma2);"
        "epsilon = sqrt(epsilon1*epsilon2)"
    )
    alch_lj.addGlobalParameter("lambda_sterics", 1.0)
    alch_lj.addGlobalParameter("softcore_alpha", 0.5)
    alch_lj.addPerParticleParameter("sigma")
    alch_lj.addPerParticleParameter("epsilon")

    for i in range(n_particles):
        q, sigma, eps = nonbonded.getParticleParameters(i)
        alch_lj.addParticle([sigma, eps])

    # Only ligand-environment interactions
    env_set = set(range(n_particles)) - ligand_set
    alch_lj.addInteractionGroup(ligand_set, env_set)

    # Copy exclusions from NonbondedForce
    print("  Copying exclusions...")
    for i in range(nonbonded.getNumExceptions()):
        p1, p2, qq, sig, eps = nonbonded.getExceptionParameters(i)
        alch_lj.addExclusion(p1, p2)

    alch_lj.setNonbondedMethod(CustomNonbondedForce.CutoffPeriodic)
    alch_lj.setCutoffDistance(nonbonded.getCutoffDistance())
    alch_lj.setUseLongRangeCorrection(True)
    alch_lj.setForceGroup(1)

    # Create alchemical electrostatics force
    print("  Creating alchemical electrostatics force...")
    alch_elec = CustomNonbondedForce(
        "lambda_elec*ONE_4PI_EPS0*charge1*charge2/r;"
        "ONE_4PI_EPS0 = 138.935456"
    )
    alch_elec.addGlobalParameter("lambda_elec", 1.0)
    alch_elec.addPerParticleParameter("charge")

    for i in range(n_particles):
        q, sigma, eps = nonbonded.getParticleParameters(i)
        alch_elec.addParticle([q])

    alch_elec.addInteractionGroup(ligand_set, env_set)

    # Copy same exclusions
    for i in range(nonbonded.getNumExceptions()):
        p1, p2, qq, sig, eps = nonbonded.getExceptionParameters(i)
        alch_elec.addExclusion(p1, p2)

    alch_elec.setNonbondedMethod(CustomNonbondedForce.CutoffPeriodic)
    alch_elec.setCutoffDistance(nonbonded.getCutoffDistance())
    alch_elec.setForceGroup(2)

    # Now zero out ligand-env in original NonbondedForce by setting ligand charges/eps to 0
    # BUT we need to preserve ligand-ligand interactions
    # Solution: Set ligand parameters to 0 and add explicit ligand-ligand 1-4 pairs back
    print("  Zeroing ligand parameters in NonbondedForce...")

    # First, store all ligand-ligand exception parameters
    ligand_exceptions = {}
    for i in range(nonbonded.getNumExceptions()):
        p1, p2, qq, sig, eps = nonbonded.getExceptionParameters(i)
        if p1 in ligand_set and p2 in ligand_set:
            ligand_exceptions[(min(p1,p2), max(p1,p2))] = (qq, sig, eps)

    # Zero out ligand particle parameters
    for i in ligand_indices:
        nonbonded.setParticleParameters(i, 0.0, 1.0*unit.nanometers, 0.0)

    # Create a separate force for ligand-ligand interactions
    print("  Creating ligand-ligand NonbondedForce...")
    ligand_nb = CustomNonbondedForce(
        "4*epsilon*((sigma/r)^12 - (sigma/r)^6) + ONE_4PI_EPS0*charge1*charge2/r;"
        "sigma = 0.5*(sigma1+sigma2);"
        "epsilon = sqrt(epsilon1*epsilon2);"
        "ONE_4PI_EPS0 = 138.935456"
    )
    ligand_nb.addPerParticleParameter("sigma")
    ligand_nb.addPerParticleParameter("epsilon")
    ligand_nb.addPerParticleParameter("charge")

    for i in range(n_particles):
        if i in ligand_set:
            q, sigma, eps = original_params[i]
            ligand_nb.addParticle([sigma, eps, q])
        else:
            ligand_nb.addParticle([1.0*unit.nanometers, 0.0, 0.0])

    ligand_nb.addInteractionGroup(ligand_set, ligand_set)

    # Copy exclusions
    for i in range(nonbonded.getNumExceptions()):
        p1, p2, qq, sig, eps = nonbonded.getExceptionParameters(i)
        ligand_nb.addExclusion(p1, p2)

    ligand_nb.setNonbondedMethod(CustomNonbondedForce.CutoffPeriodic)
    ligand_nb.setCutoffDistance(nonbonded.getCutoffDistance())
    ligand_nb.setForceGroup(3)

    # Add all forces
    system.addForce(alch_lj)
    system.addForce(alch_elec)
    system.addForce(ligand_nb)

    print("  Alchemical modifications complete")
    return system


def add_alchemical_openmmtools(system, ligand_indices):
    """Use openmmtools for alchemical setup."""
    print("Adding alchemical modifications using openmmtools...")

    # Define alchemical region
    alchemical_region = AlchemicalRegion(
        alchemical_atoms=ligand_indices,
        annihilate_electrostatics=True,
        annihilate_sterics=True,
        softcore_alpha=0.5,
        softcore_a=1,
        softcore_b=1,
        softcore_c=6
    )

    # Create factory
    factory = AbsoluteAlchemicalFactory()

    # Create alchemical system
    alchemical_system = factory.create_alchemical_system(system, alchemical_region)

    print("  Alchemical system created with openmmtools")
    return alchemical_system


def main():
    print("="*60)
    print("BUILDING ALCHEMICAL SYSTEM (PROPER)")
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
        if HAS_OPENMMTOOLS:
            alch_system = add_alchemical_openmmtools(system, lig_idx)
        else:
            alch_system = add_alchemical_manual(system, lig_idx)

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
