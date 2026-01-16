#!/usr/bin/env python
"""
Setup solvent leg FEP system for CID_10120 (Bufadienolide).

This script creates a solvated ligand system for the solvent leg of FEP
calculations, using openmmtools for alchemical transformations to ensure
compatibility with the complex leg simulations.

Key design decisions:
1. Extract ligand coordinates from the complex (consistent with complex setup)
2. Solvate using Modeller with explicit ligand topology
3. Build combined system preserving ligand force field parameters
4. Use openmmtools AbsoluteAlchemicalFactory (same as complex)

Files created:
- solvent/system.xml
- solvent/alchemical_system.xml
- solvent/topology.pdb
- solvent/positions.npy
- solvent/ligand_indices.npy

CRITICAL MATCHING REQUIREMENTS:
- softcore_alpha = 0.5 (matched by openmmtools defaults)
- softcore_a = 1, softcore_b = 1, softcore_c = 6
- PME electrostatics via ParticleOffsets
- CutoffPeriodic for sterics CustomNonbondedForce
- useLongRangeCorrection = True
- Temperature 310 K
"""
import numpy as np
from pathlib import Path
from openmm import *
from openmm.app import *
from openmm import unit

# Try to import openmmtools
try:
    from openmmtools.alchemy import AbsoluteAlchemicalFactory, AlchemicalRegion
    HAS_OPENMMTOOLS = True
    print("Using openmmtools for alchemical transformations")
except ImportError:
    HAS_OPENMMTOOLS = False
    print("WARNING: openmmtools not available - using manual approach")

# Configuration - must match complex setup exactly
TEMPERATURE = 310.0 * unit.kelvin
PRESSURE = 1.0 * unit.atmospheres
BOX_PADDING = 1.2 * unit.nanometers
IONIC_STRENGTH = 0.15 * unit.molar
CUTOFF = 1.0 * unit.nanometers
SOFTCORE_ALPHA = 0.5

BASE = Path(__file__).parent.resolve()


def extract_ligand_from_complex():
    """Extract ligand coordinates and topology from the complex system."""
    print("Extracting ligand from complex...")

    complex_dir = BASE / "wt_complex"

    # Load complex topology
    pdb = PDBFile(str(complex_dir / "topology.pdb"))

    # Load complex positions
    positions = np.load(complex_dir / "positions.npy")
    all_positions = [Vec3(p[0], p[1], p[2]) * unit.nanometer for p in positions]

    # Load ligand indices
    ligand_indices = np.load(complex_dir / "ligand_indices.npy").tolist()

    print(f"  Complex atoms: {pdb.topology.getNumAtoms()}")
    print(f"  Ligand atom indices: {ligand_indices[0]} to {ligand_indices[-1]}")
    print(f"  Ligand atoms: {len(ligand_indices)}")

    # Extract ligand atoms from topology
    # Find the ligand residue (UNK in chain B)
    ligand_chain = None
    ligand_residue = None
    for chain in pdb.topology.chains():
        for residue in chain.residues():
            if residue.name == 'UNK':
                ligand_chain = chain
                ligand_residue = residue
                break
        if ligand_residue:
            break

    if ligand_residue is None:
        raise ValueError("Could not find UNK residue in complex topology")

    print(f"  Found ligand: {ligand_residue.name} in chain {ligand_chain.id}")

    # Build ligand-only topology
    lig_topology = Topology()
    new_chain = lig_topology.addChain(id='A')
    new_residue = lig_topology.addResidue('UNK', new_chain)

    # Map old atom indices to new
    old_to_new = {}
    new_atoms = []

    for i, old_idx in enumerate(ligand_indices):
        old_atom = list(pdb.topology.atoms())[old_idx]
        new_atom = lig_topology.addAtom(old_atom.name, old_atom.element, new_residue)
        old_to_new[old_idx] = i
        new_atoms.append(new_atom)

    # Extract ligand positions
    lig_positions = [all_positions[i] for i in ligand_indices]

    print(f"  Extracted {len(lig_positions)} ligand positions")

    return lig_topology, lig_positions, ligand_indices


def load_ligand_system():
    """Load the ligand OpenMM system with force field parameters."""
    print("Loading ligand force field parameters...")

    lig_dir = BASE / "ligand_params"

    with open(lig_dir / "ligand_system.xml", 'r') as f:
        lig_system = XmlSerializer.deserialize(f.read())

    n_atoms = lig_system.getNumParticles()
    print(f"  Ligand system has {n_atoms} particles")

    return lig_system


def compute_ligand_net_charge(lig_system):
    """Compute net charge of ligand from NonbondedForce parameters."""
    for force in lig_system.getForces():
        if isinstance(force, NonbondedForce):
            total_charge = 0.0
            for i in range(force.getNumParticles()):
                q, _, _ = force.getParticleParameters(i)
                total_charge += q.value_in_unit(unit.elementary_charge)
            return round(total_charge)  # Round to nearest integer
    return 0


def solvate_ligand(lig_topology, lig_positions, lig_system):
    """Solvate the ligand in TIP3P water.

    Since Modeller can't handle UNK residues without templates, we:
    1. Create a water box separately
    2. Center the ligand in the box
    3. Delete waters that overlap with ligand

    Note: We use tip3p.xml for water GEOMETRY placement only.
    The actual TIP3P-FB PARAMETERS are applied in build_combined_system().
    """
    print("Solvating ligand...")

    # Compute ligand net charge for info
    lig_charge = compute_ligand_net_charge(lig_system)
    print(f"  Ligand net charge: {lig_charge}")

    # Get ligand coordinates
    coords = np.array([[p.x, p.y, p.z] for p in lig_positions])
    min_coords = coords.min(axis=0)
    max_coords = coords.max(axis=0)
    ligand_extent = max_coords - min_coords
    ligand_center = (min_coords + max_coords) / 2

    # Box size = ligand extent + 2*padding
    padding_nm = BOX_PADDING.value_in_unit(unit.nanometer)
    box_size = max(ligand_extent) + 2 * padding_nm
    box_size = max(box_size, 3.0)  # Minimum 3 nm box

    print(f"  Ligand extent: {ligand_extent}")
    print(f"  Box size: {box_size:.2f} nm")

    # Create empty topology and add water box
    forcefield = ForceField('tip3p.xml')
    empty_topology = Topology()
    empty_positions = []

    modeller = Modeller(empty_topology, empty_positions)
    modeller.addSolvent(
        forcefield,
        model='tip3p',
        boxSize=Vec3(box_size, box_size, box_size) * unit.nanometer
    )

    # Get water positions (strip units)
    water_positions = np.array([[p.x, p.y, p.z] for p in modeller.positions])
    # If positions have units, extract values
    if hasattr(modeller.positions[0].x, 'value_in_unit'):
        water_positions = np.array([
            [p.x.value_in_unit(unit.nanometer), p.y.value_in_unit(unit.nanometer), p.z.value_in_unit(unit.nanometer)]
            for p in modeller.positions
        ])
    water_topology = modeller.topology
    box_center = np.array([box_size / 2, box_size / 2, box_size / 2])

    # Shift ligand to box center
    ligand_shift = box_center - ligand_center
    shifted_lig_coords = coords + ligand_shift

    print(f"  Shifted ligand to box center")

    # Find waters that overlap with ligand (within 2.5 A of any ligand atom)
    min_dist = 0.25  # nm
    waters_to_keep = []
    water_atoms = list(water_topology.atoms())

    # Group water atoms into molecules (O, H, H)
    n_waters = len(water_atoms) // 3
    for w in range(n_waters):
        o_idx = w * 3
        o_pos = water_positions[o_idx]

        # Check distance to all ligand atoms
        dists = np.sqrt(np.sum((shifted_lig_coords - o_pos)**2, axis=1))
        if dists.min() > min_dist:
            waters_to_keep.append(w)

    n_removed = n_waters - len(waters_to_keep)
    print(f"  Removed {n_removed} overlapping waters")

    # Build combined topology: ligand + remaining waters
    combined_topology = Topology()
    combined_topology.setPeriodicBoxVectors(water_topology.getPeriodicBoxVectors())

    # Add ligand
    lig_chain = combined_topology.addChain(id='A')
    lig_res = combined_topology.addResidue('UNK', lig_chain)
    for atom in lig_topology.atoms():
        combined_topology.addAtom(atom.name, atom.element, lig_res)

    # Add waters
    water_chain = combined_topology.addChain(id='B')
    for w in waters_to_keep:
        water_res = combined_topology.addResidue('HOH', water_chain)
        for i in range(3):
            orig_atom = water_atoms[w * 3 + i]
            combined_topology.addAtom(orig_atom.name, orig_atom.element, water_res)

    # Build combined positions: shifted ligand + remaining waters
    combined_positions = []

    # Add ligand positions
    for pos in shifted_lig_coords:
        combined_positions.append(Vec3(pos[0], pos[1], pos[2]) * unit.nanometer)

    # Add water positions
    for w in waters_to_keep:
        for i in range(3):
            pos = water_positions[w * 3 + i]
            combined_positions.append(Vec3(pos[0], pos[1], pos[2]) * unit.nanometer)

    n_total = combined_topology.getNumAtoms()
    n_ligand = lig_topology.getNumAtoms()

    print(f"  Total atoms after solvation: {n_total}")
    print(f"  Ligand: {n_ligand}")
    print(f"  Water molecules: {len(waters_to_keep)}")

    box = combined_topology.getPeriodicBoxVectors()
    box_x = box[0][0].value_in_unit(unit.nanometer)
    box_y = box[1][1].value_in_unit(unit.nanometer)
    box_z = box[2][2].value_in_unit(unit.nanometer)
    print(f"  Box: {box_x:.2f} x {box_y:.2f} x {box_z:.2f} nm")
    print(f"  NOTE: No ions added - pure water solvation")

    return combined_topology, combined_positions


def build_combined_system(solvated_topology, solvated_positions, lig_system):
    """Build a combined system using ligand parameters + water parameters."""
    print("Building combined system...")

    n_lig_atoms = lig_system.getNumParticles()
    n_total = solvated_topology.getNumAtoms()
    n_water = n_total - n_lig_atoms

    # Create the combined system
    system = System()

    # Set box vectors
    box_vectors = solvated_topology.getPeriodicBoxVectors()
    system.setDefaultPeriodicBoxVectors(*box_vectors)

    # Get ligand NonbondedForce for parameters
    lig_nb = None
    for force in lig_system.getForces():
        if isinstance(force, NonbondedForce):
            lig_nb = force
            break

    # Create system for water only to get parameters
    water_ff = ForceField('tip3p.xml')

    # Create a minimal water topology for parameter extraction
    water_topology = Topology()
    water_chain = water_topology.addChain()

    # Count waters and ions from solvated topology
    atom_list = list(solvated_topology.atoms())
    residue_list = list(solvated_topology.residues())

    # Add particles with correct masses
    print("  Adding particles...")
    for i, atom in enumerate(atom_list):
        if i < n_lig_atoms:
            mass = lig_system.getParticleMass(i)
        else:
            mass = atom.element.mass
        system.addParticle(mass)

    # Create NonbondedForce - MUST match complex leg exactly
    nonbonded = NonbondedForce()
    nonbonded.setNonbondedMethod(NonbondedForce.PME)
    nonbonded.setCutoffDistance(CUTOFF)
    nonbonded.setEwaldErrorTolerance(0.0005)
    nonbonded.setUseSwitchingFunction(False)  # Must match complex (useSwitchingFunction="0")
    nonbonded.setUseDispersionCorrection(True)

    # Add ligand particles
    print("  Adding ligand nonbonded parameters...")
    for i in range(n_lig_atoms):
        q, sig, eps = lig_nb.getParticleParameters(i)
        nonbonded.addParticle(q, sig, eps)

    # Add water/ion particles
    print("  Adding water/ion nonbonded parameters...")
    # TIP3P-FB water parameters (MUST match complex!)
    # From wt_complex/system.xml:
    # Oxygen: q=-0.848448690103, sigma=0.317796456355, epsilon=0.652143528104
    # Hydrogen: q=+0.4242243450515, sigma=1.0 (dummy), epsilon=0
    # Na+: q=1, sigma=0.2439280690268249, epsilon=0.3658460312
    # Cl-: q=-1, sigma=0.4477656957373345, epsilon=0.14891274399999999

    tip3p_fb_o = (-0.848448690103*unit.elementary_charge, 0.317796456355*unit.nanometer, 0.652143528104*unit.kilojoule_per_mole)
    tip3p_fb_h = (0.4242243450515*unit.elementary_charge, 1.0*unit.nanometer, 0.0*unit.kilojoule_per_mole)
    na_params = (1.0*unit.elementary_charge, 0.2439280690268249*unit.nanometer, 0.3658460312*unit.kilojoule_per_mole)
    cl_params = (-1.0*unit.elementary_charge, 0.4477656957373345*unit.nanometer, 0.14891274399999999*unit.kilojoule_per_mole)

    for i in range(n_lig_atoms, n_total):
        atom = atom_list[i]
        res_name = atom.residue.name
        atom_name = atom.name

        if res_name == 'HOH':
            if atom_name == 'O':
                nonbonded.addParticle(*tip3p_fb_o)
            else:  # H1, H2
                nonbonded.addParticle(*tip3p_fb_h)
        elif res_name in ('NA', 'Na+'):
            nonbonded.addParticle(*na_params)
        elif res_name in ('CL', 'Cl-'):
            nonbonded.addParticle(*cl_params)
        else:
            # This should only happen if there's an unexpected residue type
            raise ValueError(f"Unknown residue: {res_name} (atom {atom_name}, index {i})")

    # Add ligand exceptions (1-4 scaling)
    print("  Adding ligand exceptions...")
    for i in range(lig_nb.getNumExceptions()):
        p1, p2, qq, sig, eps = lig_nb.getExceptionParameters(i)
        nonbonded.addException(p1, p2, qq, sig, eps)

    # Add water exclusions (within each water molecule)
    print("  Adding water exclusions...")
    for res in residue_list:
        if res.name == 'HOH':
            atoms = list(res.atoms())
            atom_indices = [a.index for a in atoms]
            # Exclude all pairs within water molecule
            for j in range(len(atom_indices)):
                for k in range(j+1, len(atom_indices)):
                    nonbonded.addException(atom_indices[j], atom_indices[k], 0, 1*unit.nanometer, 0)

    system.addForce(nonbonded)

    # Add ligand bonded forces
    print("  Adding ligand bonded forces...")
    for force in lig_system.getForces():
        if isinstance(force, HarmonicBondForce):
            new_force = HarmonicBondForce()
            for i in range(force.getNumBonds()):
                p1, p2, length, k = force.getBondParameters(i)
                new_force.addBond(p1, p2, length, k)
            system.addForce(new_force)
            print(f"    Bonds: {new_force.getNumBonds()}")

        elif isinstance(force, HarmonicAngleForce):
            new_force = HarmonicAngleForce()
            for i in range(force.getNumAngles()):
                p1, p2, p3, angle, k = force.getAngleParameters(i)
                new_force.addAngle(p1, p2, p3, angle, k)
            system.addForce(new_force)
            print(f"    Angles: {new_force.getNumAngles()}")

        elif isinstance(force, PeriodicTorsionForce):
            new_force = PeriodicTorsionForce()
            for i in range(force.getNumTorsions()):
                p1, p2, p3, p4, periodicity, phase, k = force.getTorsionParameters(i)
                new_force.addTorsion(p1, p2, p3, p4, periodicity, phase, k)
            system.addForce(new_force)
            print(f"    Torsions: {new_force.getNumTorsions()}")

    # Add water constraints (TIP3P-FB geometry - MUST match complex!)
    print("  Adding water constraints...")
    # TIP3P-FB geometry from wt_complex/system.xml constraints
    # Note: TIP3P-FB has DIFFERENT geometry from standard TIP3P!
    o_h_dist = 0.101181082494 * unit.nanometer  # O-H distance (TIP3P-FB)
    h_h_dist = 0.16386837572186558 * unit.nanometer  # H-H distance (TIP3P-FB)

    for res in residue_list:
        if res.name == 'HOH':
            atoms = {a.name: a.index for a in res.atoms()}
            if 'O' in atoms and 'H1' in atoms and 'H2' in atoms:
                system.addConstraint(atoms['O'], atoms['H1'], o_h_dist)
                system.addConstraint(atoms['O'], atoms['H2'], o_h_dist)
                system.addConstraint(atoms['H1'], atoms['H2'], h_h_dist)

    # Add ligand constraints
    for i in range(lig_system.getNumConstraints()):
        p1, p2, dist = lig_system.getConstraintParameters(i)
        system.addConstraint(p1, p2, dist)

    print(f"    Total constraints: {system.getNumConstraints()}")

    # Add barostat for NPT equilibration (same as febuxostat solvent)
    barostat = MonteCarloBarostat(PRESSURE, TEMPERATURE, 25)
    system.addForce(barostat)
    print("  Added MonteCarloBarostat (frequency=25)")

    # Add CMMotionRemover to match complex
    cmm = CMMotionRemover(1)
    system.addForce(cmm)
    print("  Added CMMotionRemover")

    print(f"  System complete: {system.getNumParticles()} particles, {system.getNumForces()} forces")

    return system


def add_alchemical_openmmtools(system, ligand_indices):
    """Use openmmtools to create alchemical system (preferred method)."""
    print("Creating alchemical system with openmmtools...")

    alchemical_region = AlchemicalRegion(
        alchemical_atoms=ligand_indices,
        annihilate_electrostatics=True,
        annihilate_sterics=True,
        softcore_alpha=SOFTCORE_ALPHA,
        softcore_a=1,
        softcore_b=1,
        softcore_c=6
    )

    factory = AbsoluteAlchemicalFactory()
    alchemical_system = factory.create_alchemical_system(system, alchemical_region)

    print(f"  Alchemical system created")
    print(f"  Forces: {alchemical_system.getNumForces()}")

    # List forces
    for i, f in enumerate(alchemical_system.getForces()):
        print(f"    {i}: {f.__class__.__name__}")

    return alchemical_system


def add_alchemical_manual(system, ligand_indices):
    """Manual alchemical setup - NOT RECOMMENDED!

    This is disabled because the manual approach has known issues:
    1. Ligand-ligand electrostatics use cutoff (no PME), causing mismatch with complex leg
    2. Potential double-counting with ParticleOffsets

    openmmtools is REQUIRED for correct alchemical setup.
    """
    raise RuntimeError(
        "openmmtools is REQUIRED for alchemical system setup.\n"
        "The manual fallback has known issues that cause incorrect free energies.\n"
        "Please install openmmtools: pip install openmmtools\n"
        "Or: conda install -c conda-forge openmmtools"
    )


def verify_lambda_schedule():
    """Verify lambda schedule exists and matches complex."""
    print("\nVerifying lambda schedule...")

    solvent_schedule = BASE / "solvent" / "lambda_schedule.npy"
    complex_schedule = BASE / "wt_complex" / "lambda_schedule.npy"

    if not solvent_schedule.exists():
        raise FileNotFoundError(f"Missing: {solvent_schedule}")

    sol_sched = np.load(solvent_schedule)
    cplx_sched = np.load(complex_schedule)

    print(f"  Solvent shape: {sol_sched.shape}")
    print(f"  Complex shape: {cplx_sched.shape}")

    if sol_sched.shape == cplx_sched.shape:
        if np.allclose(sol_sched, cplx_sched):
            print("  [OK] Schedules match exactly")
        else:
            if np.allclose(sol_sched[:, :2], cplx_sched[:, :2]):
                print("  [OK] Electrostatics and sterics match (restraints may differ)")
            else:
                print("  WARNING: Schedules differ!")
    else:
        print("  WARNING: Schedule shapes differ!")

    return sol_sched


def main():
    print("=" * 70)
    print("SETUP SOLVENT LEG FEP FOR CID_10120 (DIHYDROCORTISOL)")
    print("=" * 70)
    print(f"openmmtools available: {HAS_OPENMMTOOLS}")
    print()

    # Fail early if openmmtools is not available
    if not HAS_OPENMMTOOLS:
        raise RuntimeError(
            "openmmtools is REQUIRED for correct alchemical setup.\n"
            "Please install: pip install openmmtools\n"
            "Or: conda install -c conda-forge openmmtools"
        )

    solvent_dir = BASE / "solvent"
    solvent_dir.mkdir(exist_ok=True)

    # 1. Extract ligand from complex (ensures consistency)
    lig_topology, lig_positions, complex_lig_indices = extract_ligand_from_complex()

    # 2. Load ligand force field parameters
    lig_system = load_ligand_system()

    # Validate atom counts match
    n_lig_topology = lig_topology.getNumAtoms()
    n_lig_system = lig_system.getNumParticles()
    if n_lig_topology != n_lig_system:
        raise ValueError(
            f"Ligand atom count mismatch!\n"
            f"  From complex topology: {n_lig_topology}\n"
            f"  From ligand_system.xml: {n_lig_system}\n"
            f"This indicates the ligand_params may not match the complex."
        )
    print(f"  [OK] Ligand atom counts match: {n_lig_topology}")

    # 3. Solvate ligand
    solvated_topology, solvated_positions = solvate_ligand(lig_topology, lig_positions, lig_system)

    # 4. Build combined system
    n_lig = lig_topology.getNumAtoms()
    system = build_combined_system(solvated_topology, solvated_positions, lig_system)

    # Ligand indices in the new system (always 0 to n_lig-1)
    ligand_indices = list(range(n_lig))

    # 5. Save base system
    print("\nSaving base system...")
    with open(solvent_dir / "system.xml", 'w') as f:
        f.write(XmlSerializer.serialize(system))
    print(f"  Saved: {solvent_dir / 'system.xml'}")

    # 6. Create alchemical system
    if HAS_OPENMMTOOLS:
        alch_system = add_alchemical_openmmtools(system, ligand_indices)
    else:
        alch_system = add_alchemical_manual(system, ligand_indices)

    # 7. Save alchemical system
    print("\nSaving alchemical system...")
    with open(solvent_dir / "alchemical_system.xml", 'w') as f:
        f.write(XmlSerializer.serialize(alch_system))
    print(f"  Saved: {solvent_dir / 'alchemical_system.xml'}")

    # 8. Save topology
    print("\nSaving topology...")
    # Convert positions to plain Vec3 list (units were applied in solvate_ligand)
    # PDBFile expects plain Vec3 with separate unit specification or no units
    clean_positions = []
    for p in solvated_positions:
        # Each p is Quantity[Vec3] from Vec3(...) * unit.nanometer
        # Extract the values in nanometers
        if hasattr(p, 'value_in_unit'):
            pos_nm = p.value_in_unit(unit.nanometer)
            clean_positions.append(Vec3(pos_nm[0], pos_nm[1], pos_nm[2]))
        else:
            clean_positions.append(p)

    with open(solvent_dir / "topology.pdb", 'w') as f:
        # Pass positions in nanometers
        PDBFile.writeFile(solvated_topology, clean_positions * unit.nanometer, f)
    print(f"  Saved: {solvent_dir / 'topology.pdb'}")

    # 9. Save positions (as nanometers)
    print("\nSaving positions...")
    # Use already cleaned positions
    pos_array = np.array([[p.x, p.y, p.z] for p in clean_positions])
    np.save(solvent_dir / "positions.npy", pos_array)
    print(f"  Saved: {solvent_dir / 'positions.npy'}")
    print(f"  Shape: {pos_array.shape}")

    # 10. Save ligand indices
    print("\nSaving ligand indices...")
    np.save(solvent_dir / "ligand_indices.npy", np.array(ligand_indices))
    print(f"  Saved: {solvent_dir / 'ligand_indices.npy'}")

    # 11. Verify lambda schedule
    verify_lambda_schedule()

    # Summary
    print("\n" + "=" * 70)
    print("SOLVENT SYSTEM SETUP COMPLETE")
    print("=" * 70)
    print(f"\nFiles in {solvent_dir}:")
    print(f"  - system.xml")
    print(f"  - alchemical_system.xml")
    print(f"  - topology.pdb")
    print(f"  - positions.npy")
    print(f"  - ligand_indices.npy")
    print(f"\nSystem statistics:")
    print(f"  Total atoms: {system.getNumParticles()}")
    print(f"  Ligand atoms: {len(ligand_indices)}")
    print(f"  Solvent atoms: {system.getNumParticles() - len(ligand_indices)}")
    print(f"\nAlchemical parameters:")
    print(f"  softcore_alpha = {SOFTCORE_ALPHA}")
    print(f"  NO Boresch restraints (solvent leg)")

    # Check run_window.py files
    print("\nChecking window directories...")
    missing = []
    for i in range(20):
        if not (solvent_dir / f"window_{i:02d}" / "run_window.py").exists():
            missing.append(i)

    if missing:
        print(f"  WARNING: Missing run_window.py for windows: {missing[:5]}...")
    else:
        print("  [OK] All 20 windows have run_window.py")

    print("\nReady for FEP simulations!")


if __name__ == "__main__":
    main()
