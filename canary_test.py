#!/usr/bin/env python
"""
CANARY TEST - Find minimum softcore_alpha that works for window_19 (sterics=0.0)

Tests: 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 1.00
Uses solvent system (smallest, fastest)

For each alpha:
1. Regenerate alchemical system with that softcore_alpha
2. Run window_19 with short production
3. Report pass/fail
"""
import os
import sys
import numpy as np
import argparse
import shutil

from openmm import XmlSerializer, LangevinMiddleIntegrator, Platform, MonteCarloBarostat, Vec3
from openmm.app import PDBFile, Simulation, ForceField
from openmm import unit

# Energy thresholds for sanity checks
MAX_ENERGY_KJ = 1e8  # If energy > 100 million kJ/mol, consider unstable

# Test parameters
TEMPERATURE = 310.0 * unit.kelvin
FRICTION = 2.0 / unit.picosecond
TIMESTEP = 1.0 * unit.femtoseconds
EQUIL_STEPS = 500000   # 500 ps equilibration (longer for stability test)
PROD_STEPS = 200000    # 200 ps production
ENERGY_INTERVAL = 1000
USE_NPT = True  # Keep barostat for realistic test


def regenerate_alchemical_system(base_dir, alpha_value, output_dir):
    """Regenerate alchemical system with specified softcore_alpha."""

    print(f"  Regenerating alchemical system with softcore_alpha={alpha_value}...")

    # Import openmmtools here to avoid issues
    try:
        from openmmtools.alchemy import AbsoluteAlchemicalFactory, AlchemicalRegion
    except ImportError:
        print("[FAIL] OpenMMTools not installed!")
        return None

    # Load original system
    system_path = os.path.join(base_dir, "system.xml")
    topology_path = os.path.join(base_dir, "topology.pdb")

    with open(system_path, "r") as f:
        system = XmlSerializer.deserialize(f.read())

    pdb = PDBFile(topology_path)
    topology = pdb.topology

    # Find ligand atoms (non-water, non-ion residues that aren't protein)
    EXCLUDE = {"HOH", "WAT", "TIP3", "SOL", "NA", "CL", "K", "MG", "CA", "ZN", "NA+", "CL-"}

    ligand_atom_indices = []
    for res in topology.residues():
        if res.name.upper() not in EXCLUDE:
            for atom in res.atoms():
                ligand_atom_indices.append(atom.index)

    if not ligand_atom_indices:
        print("  [FAIL] No ligand atoms found!")
        return None

    print(f"  Found {len(ligand_atom_indices)} ligand atoms")

    # Create alchemical region with specified alpha
    alchemical_region = AlchemicalRegion(
        alchemical_atoms=ligand_atom_indices,
        annihilate_electrostatics=True,
        annihilate_sterics=True,
        softcore_alpha=alpha_value,
        softcore_a=1,
        softcore_b=1,
        softcore_c=6
    )

    # Create alchemical factory
    factory = AbsoluteAlchemicalFactory(
        alchemical_pme_treatment='direct-space'
    )
    alchemical_system = factory.create_alchemical_system(system, alchemical_region)

    # Save to output directory
    os.makedirs(output_dir, exist_ok=True)

    alch_path = os.path.join(output_dir, "alchemical_system.xml")
    with open(alch_path, "w") as f:
        f.write(XmlSerializer.serialize(alchemical_system))

    # Copy other needed files
    for fname in ["topology.pdb", "positions.npy", "lambda_schedule.npy"]:
        src = os.path.join(base_dir, fname)
        dst = os.path.join(output_dir, fname)
        if os.path.exists(src):
            shutil.copy(src, dst)

    print(f"  [OK] Alchemical system saved with softcore_alpha={alpha_value}")
    return alch_path


def get_box_vectors_from_pdb(pdb_path):
    """Extract box vectors from PDB CRYST1 record."""
    with open(pdb_path, 'r') as f:
        for line in f:
            if line.startswith('CRYST1'):
                # CRYST1 format: a, b, c in Angstroms, angles in degrees
                a = float(line[6:15]) / 10.0  # Convert to nm
                b = float(line[15:24]) / 10.0
                c = float(line[24:33]) / 10.0
                # Assume orthorhombic box (90 degree angles)
                return (
                    Vec3(a, 0, 0) * unit.nanometer,
                    Vec3(0, b, 0) * unit.nanometer,
                    Vec3(0, 0, c) * unit.nanometer
                )
    return None


def run_canary_window(output_dir, window_idx=19):
    """Run canary test for window_19 (sterics=0.0)."""

    print(f"  Running window {window_idx} test...")

    # Load alchemical system
    alch_path = os.path.join(output_dir, "alchemical_system.xml")
    with open(alch_path, "r") as f:
        system = XmlSerializer.deserialize(f.read())

    # Optionally remove barostat
    if not USE_NPT:
        for i in range(system.getNumForces()):
            force = system.getForce(i)
            if isinstance(force, MonteCarloBarostat):
                system.removeForce(i)
                print("  Running NVT (barostat removed)")
                break
    else:
        print("  Running NPT (barostat kept)")

    # Load topology and positions
    pdb_path = os.path.join(output_dir, "topology.pdb")
    pdb = PDBFile(pdb_path)
    topology = pdb.topology

    # Load positions (stored in nanometers)
    positions = np.load(os.path.join(output_dir, "positions.npy"))
    positions = positions * unit.nanometer

    # Get box vectors from PDB
    box_vectors = get_box_vectors_from_pdb(pdb_path)
    if box_vectors is None:
        return False, "No CRYST1 record in topology.pdb"

    # Load lambda schedule
    lambda_schedule = np.load(os.path.join(output_dir, "lambda_schedule.npy"))
    n_windows = len(lambda_schedule)
    lam_e, lam_s, lam_r = lambda_schedule[window_idx]

    print(f"  Lambda: elec={lam_e:.4f}, sterics={lam_s:.4f}, restraints={lam_r:.4f}")

    # Create integrator
    integrator = LangevinMiddleIntegrator(TEMPERATURE, FRICTION, TIMESTEP)

    # Select platform
    try:
        platform = Platform.getPlatformByName('CUDA')
        properties = {'Precision': 'mixed'}
    except Exception:
        platform = Platform.getPlatformByName('CPU')
        properties = {}

    # Create simulation
    simulation = Simulation(topology, system, integrator, platform, properties)
    context = simulation.context

    # Set box vectors FIRST, then positions
    context.setPeriodicBoxVectors(*box_vectors)
    context.setPositions(positions)

    # Validate lambda parameters exist (STRICT - fail if missing)
    params = context.getParameters()
    if 'lambda_electrostatics' not in params:
        return False, "lambda_electrostatics parameter missing from system"
    if 'lambda_sterics' not in params:
        return False, "lambda_sterics parameter missing from system"

    # Set lambda parameters (solvent has no restraints)
    context.setParameter('lambda_electrostatics', float(lam_e))
    context.setParameter('lambda_sterics', float(lam_s))

    # Initialize velocities
    context.setVelocitiesToTemperature(TEMPERATURE)

    # Minimize
    print("  Minimizing...")
    try:
        simulation.minimizeEnergy(maxIterations=1000)
        state = context.getState(getEnergy=True)
        energy = state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
        print(f"  Energy after minimization: {energy:.1f} kJ/mol")
        if not np.isfinite(energy):
            return False, "NaN after minimization"
        if abs(energy) > MAX_ENERGY_KJ:
            return False, f"Energy too high after minimization: {energy:.1e} kJ/mol"
    except Exception as e:
        return False, f"Minimization failed: {e}"

    # Equilibration
    equil_time_ps = EQUIL_STEPS * TIMESTEP.value_in_unit(unit.picosecond)
    print(f"  Equilibrating ({EQUIL_STEPS} steps = {equil_time_ps:.0f} ps)...")
    try:
        chunk = 10000
        for i in range(0, EQUIL_STEPS, chunk):
            simulation.step(min(chunk, EQUIL_STEPS - i))
            state = context.getState(getEnergy=True)
            energy = state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
            if not np.isfinite(energy):
                return False, f"NaN during equilibration at step {i + chunk}"
            if abs(energy) > MAX_ENERGY_KJ:
                return False, f"Energy too high at step {i + chunk}: {energy:.1e} kJ/mol"
        print(f"  Equilibration complete, E = {energy:.1f} kJ/mol")
    except Exception as e:
        return False, f"Equilibration failed: {e}"

    # Test u_nk evaluation
    print("  Testing u_nk evaluation...")
    kT = (unit.MOLAR_GAS_CONSTANT_R * TEMPERATURE).value_in_unit(unit.kilojoules_per_mole)
    nan_count = 0
    energies = []

    for k in range(n_windows):
        lam_e_k, lam_s_k, lam_r_k = lambda_schedule[k]
        context.setParameter('lambda_electrostatics', float(lam_e_k))
        context.setParameter('lambda_sterics', float(lam_s_k))

        state = context.getState(getEnergy=True)
        e_kj = state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
        u = e_kj / kT
        energies.append(e_kj)

        if not np.isfinite(u):
            nan_count += 1
            print(f"    window {k}: NaN/Inf (sterics={lam_s_k:.2f})")

    # Restore lambda
    context.setParameter('lambda_electrostatics', float(lam_e))
    context.setParameter('lambda_sterics', float(lam_s))

    if nan_count > 0:
        return False, f"{nan_count}/{n_windows} NaN in u_nk evaluation"

    # Report energy range (large values are OK for MBAR)
    e_min, e_max = min(energies), max(energies)
    print(f"  u_nk: all {n_windows} states finite, E range: {e_min:.1f} to {e_max:.1f} kJ/mol")

    # Short production
    prod_time_ps = PROD_STEPS * TIMESTEP.value_in_unit(unit.picosecond)
    print(f"  Production ({PROD_STEPS} steps = {prod_time_ps:.0f} ps)...")
    try:
        n_samples = PROD_STEPS // ENERGY_INTERVAL
        for sample in range(n_samples):
            simulation.step(ENERGY_INTERVAL)
            state = context.getState(getEnergy=True)
            energy = state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
            if not np.isfinite(energy):
                return False, f"NaN during production at sample {sample + 1}"
            if abs(energy) > MAX_ENERGY_KJ:
                return False, f"Energy too high during production: {energy:.1e} kJ/mol"
        print(f"  Production complete, E = {energy:.1f} kJ/mol")
    except Exception as e:
        return False, f"Production failed: {e}"

    return True, "All tests passed"


def run_canary_test(alpha_value, base_dir, canary_dir):
    """Run complete canary test for one alpha value."""

    print(f"\n{'='*60}")
    print(f"CANARY TEST: softcore_alpha = {alpha_value}")
    print(f"{'='*60}")

    output_dir = os.path.join(canary_dir, f"alpha_{alpha_value:.2f}")

    # Regenerate alchemical system
    alch_path = regenerate_alchemical_system(base_dir, alpha_value, output_dir)
    if alch_path is None:
        return False, "Failed to regenerate alchemical system"

    # Run canary window
    success, message = run_canary_window(output_dir, window_idx=19)

    if success:
        print(f"\n[PASS] softcore_alpha={alpha_value}: {message}")
    else:
        print(f"\n[FAIL] softcore_alpha={alpha_value}: {message}")

    return success, message


def main():
    parser = argparse.ArgumentParser(description="Canary test for softcore_alpha")
    parser.add_argument("--alpha", type=float, help="Test single alpha value")
    parser.add_argument("--base", default="C:/Users/vasud/nod2-screening-data/fep_pmx/solvent",
                        help="Base solvent directory")
    parser.add_argument("--output", default="C:/Users/vasud/nod2-screening-data/fep_pmx_canary",
                        help="Output directory for canary tests")

    args = parser.parse_args()

    if args.alpha:
        # Test single alpha
        success, message = run_canary_test(args.alpha, args.base, args.output)
        sys.exit(0 if success else 1)
    else:
        # Test all alpha values
        alpha_values = [0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 1.00]
        results = {}

        print("=" * 60)
        print("CANARY TEST SUITE")
        print("=" * 60)
        print(f"Testing alpha values: {alpha_values}")
        print(f"Test window: 19 (sterics=0.0, fully decoupled)")
        print()

        for alpha in alpha_values:
            success, message = run_canary_test(alpha, args.base, args.output)
            results[alpha] = (success, message)

        # Summary
        print("\n" + "=" * 60)
        print("CANARY TEST RESULTS")
        print("=" * 60)

        min_passing = None
        for alpha in alpha_values:
            success, message = results[alpha]
            status = "PASS" if success else "FAIL"
            print(f"  alpha={alpha:.2f}: [{status}] {message}")
            if success and min_passing is None:
                min_passing = alpha

        print()
        if min_passing:
            print(f"RECOMMENDATION: Use softcore_alpha = {min_passing}")
            print(f"  This is the LOWEST alpha that passed all tests.")
        else:
            print("WARNING: All tests failed! May need alpha > 1.0 or other changes.")

        sys.exit(0 if min_passing else 1)


if __name__ == "__main__":
    main()
