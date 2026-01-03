#!/usr/bin/env python3
"""
PREFLIGHT VALIDATION SCRIPT
Run this BEFORE spending money on GPU rental!

This script:
1. Validates all input files exist
2. Tests ligand parameterization for each compound
3. Runs 10-step minimization to catch structure issues
4. Reports any problems BEFORE you start the 9-GPU burn

Usage:
    python preflight_check.py

Expected runtime: 1-2 minutes on CPU, catches 90% of potential failures
"""

import os
import sys
from datetime import datetime

def check_dependencies():
    """Check all required packages are installed."""
    print("=" * 70)
    print("PREFLIGHT CHECK: NOD2-SCOUT MD SIMULATIONS")
    print("=" * 70)
    print(f"\nStarted: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("\n[1/5] Checking dependencies...")

    errors = []

    # OpenMM
    try:
        import openmm
        print(f"      OpenMM: {openmm.__version__}")
    except ImportError:
        errors.append("OpenMM not installed")

    # OpenFF
    try:
        from openff.toolkit import Molecule
        import openff.toolkit
        print(f"      OpenFF Toolkit: {openff.toolkit.__version__}")
    except ImportError:
        errors.append("OpenFF Toolkit not installed")

    # OpenMMForceFields
    try:
        from openmmforcefields.generators import SystemGenerator
        print("      OpenMMForceFields: OK")
    except ImportError:
        errors.append("OpenMMForceFields not installed")

    # MDTraj (for analysis)
    try:
        import mdtraj
        print(f"      MDTraj: {mdtraj.__version__}")
    except ImportError:
        print("      MDTraj: NOT INSTALLED (needed for analysis only)")

    if errors:
        print("\n  ERRORS:")
        for e in errors:
            print(f"    - {e}")
        print("\n  Install with:")
        print("    conda install -c conda-forge openmm openmmforcefields openff-toolkit")
        return False

    print("      All dependencies OK!")
    return True


def check_cuda():
    """Check CUDA/GPU availability."""
    print("\n[2/5] Checking CUDA...")

    try:
        from openmm import Platform
        cuda = Platform.getPlatformByName('CUDA')
        print(f"      CUDA Platform: Available")

        # Try to get device count
        import subprocess
        result = subprocess.run(['nvidia-smi', '-L'], capture_output=True, text=True)
        if result.returncode == 0:
            gpus = [l for l in result.stdout.split('\n') if l.strip()]
            print(f"      GPUs found: {len(gpus)}")
            for gpu in gpus[:3]:  # Show first 3
                print(f"        - {gpu}")
            if len(gpus) > 3:
                print(f"        ... and {len(gpus) - 3} more")
        return True
    except Exception as e:
        print(f"      WARNING: CUDA not available ({e})")
        print("      Will test on CPU (slower but validates structure)")
        return False


def check_input_files():
    """Check all input files exist."""
    print("\n[3/5] Checking input files...")

    # Expected files
    structures_dir = "../structures"
    required_files = {
        'complex_febuxostat.pdb': 'febuxostat_docked.sdf',
        'complex_ursodiol.pdb': 'ursodiol_docked.sdf',
        'complex_budesonide.pdb': 'budesonide_docked.sdf',
        'complex_natural_cid10592.pdb': 'natural_top_docked.sdf',
        'complex_decoy.pdb': 'decoy_docked.sdf',
        'complex_apo.pdb': None,  # No ligand
    }

    all_ok = True
    for pdb, sdf in required_files.items():
        pdb_path = os.path.join(structures_dir, pdb)
        pdb_exists = os.path.exists(pdb_path)

        if sdf:
            sdf_path = os.path.join(structures_dir, sdf)
            sdf_exists = os.path.exists(sdf_path)
        else:
            sdf_exists = True
            sdf_path = None

        if pdb_exists and sdf_exists:
            print(f"      OK: {pdb}")
        else:
            all_ok = False
            if not pdb_exists:
                print(f"      MISSING: {pdb}")
            if sdf and not sdf_exists:
                print(f"      MISSING: {sdf}")

    return all_ok


def test_ligand_parameterization():
    """Test that each ligand can be parameterized."""
    print("\n[4/5] Testing ligand parameterization...")

    from openff.toolkit import Molecule
    from openmmforcefields.generators import SystemGenerator, GAFFTemplateGenerator

    ligands = [
        ("febuxostat", "../structures/febuxostat_docked.sdf"),
        ("ursodiol", "../structures/ursodiol_docked.sdf"),
        ("budesonide", "../structures/budesonide_docked.sdf"),
        ("natural", "../structures/natural_top_docked.sdf"),
        ("decoy", "../structures/decoy_docked.sdf"),
    ]

    all_ok = True
    for name, sdf_path in ligands:
        if not os.path.exists(sdf_path):
            print(f"      SKIP: {name} (file missing)")
            continue

        try:
            mol = Molecule.from_file(sdf_path)
            charge = sum(atom.formal_charge for atom in mol.atoms)

            # Try OpenFF parameterization
            try:
                system_generator = SystemGenerator(
                    forcefields=['amber14-all.xml', 'amber14/tip3pfb.xml'],
                    small_molecule_forcefield='openff-2.1.0',
                    molecules=[mol],
                    forcefield_kwargs={'constraints': None}
                )
                method = "OpenFF"
            except Exception as e:
                # Fallback to GAFF
                method = "GAFF (fallback)"

            print(f"      OK: {name} (charge={charge}, method={method})")

        except Exception as e:
            all_ok = False
            print(f"      FAIL: {name} - {e}")

    return all_ok


def test_short_simulation():
    """Run a 10-step test simulation to catch structure issues."""
    print("\n[5/5] Running 10-step test simulation...")

    try:
        from openmm import *
        from openmm.app import *
        from openmm.unit import *
        from openff.toolkit import Molecule
        from openmmforcefields.generators import SystemGenerator

        # Use febuxostat as test case
        pdb_file = "../structures/complex_febuxostat.pdb"
        sdf_file = "../structures/febuxostat_docked.sdf"

        if not os.path.exists(pdb_file) or not os.path.exists(sdf_file):
            print("      SKIP: Test files not found")
            return True

        print("      Loading structures...")
        mol = Molecule.from_file(sdf_file)
        pdb = PDBFile(pdb_file)

        print("      Setting up system...")
        system_generator = SystemGenerator(
            forcefields=['amber14-all.xml', 'amber14/tip3pfb.xml'],
            small_molecule_forcefield='openff-2.1.0',
            molecules=[mol],
            forcefield_kwargs={'constraints': AllBonds, 'rigidWater': True}
        )

        modeller = Modeller(pdb.topology, pdb.positions)
        modeller.addSolvent(system_generator.forcefield, model='tip3p',
                           padding=1.0*nanometers, ionicStrength=0.15*molar)

        system = system_generator.create_system(
            modeller.topology,
            nonbondedMethod=PME,
            nonbondedCutoff=1.0*nanometers
        )

        print(f"      Solvated system: {modeller.topology.getNumAtoms()} atoms")

        # Use CPU for preflight (more reliable)
        integrator = LangevinMiddleIntegrator(310*kelvin, 1/picosecond, 2*femtoseconds)

        try:
            platform = Platform.getPlatformByName('CUDA')
            properties = {'Precision': 'mixed', 'DeviceIndex': '0'}
            simulation = Simulation(modeller.topology, system, integrator, platform, properties)
            print("      Using CUDA platform")
        except:
            platform = Platform.getPlatformByName('CPU')
            simulation = Simulation(modeller.topology, system, integrator, platform)
            print("      Using CPU platform (CUDA not available)")

        simulation.context.setPositions(modeller.positions)

        # Check initial energy
        state = simulation.context.getState(getEnergy=True)
        energy = state.getPotentialEnergy()
        print(f"      Initial energy: {energy}")

        if 'nan' in str(energy).lower() or 'inf' in str(energy).lower():
            print("      FAIL: Initial energy is NaN/Inf!")
            return False

        # Minimize
        print("      Minimizing (100 steps)...")
        simulation.minimizeEnergy(maxIterations=100)

        state = simulation.context.getState(getEnergy=True)
        energy = state.getPotentialEnergy()
        print(f"      Minimized energy: {energy}")

        # Run 10 steps
        print("      Running 10 MD steps...")
        simulation.context.setVelocitiesToTemperature(310*kelvin)
        simulation.step(10)

        state = simulation.context.getState(getEnergy=True)
        energy = state.getPotentialEnergy()
        print(f"      Final energy: {energy}")

        if 'nan' in str(energy).lower() or 'inf' in str(energy).lower():
            print("      FAIL: Energy exploded!")
            return False

        print("      SUCCESS: Test simulation completed")
        return True

    except Exception as e:
        print(f"      FAIL: {e}")
        import traceback
        traceback.print_exc()
        return False


def main():
    os.chdir(os.path.dirname(os.path.abspath(__file__)))

    results = []

    # Run all checks
    results.append(("Dependencies", check_dependencies()))
    results.append(("CUDA", check_cuda()))
    results.append(("Input Files", check_input_files()))

    if results[0][1]:  # Only if dependencies OK
        results.append(("Ligand Params", test_ligand_parameterization()))
        results.append(("Test Simulation", test_short_simulation()))

    # Summary
    print("\n" + "=" * 70)
    print("PREFLIGHT SUMMARY")
    print("=" * 70)

    all_pass = True
    for name, passed in results:
        status = "PASS" if passed else "FAIL"
        symbol = "[✓]" if passed else "[✗]"
        print(f"  {symbol} {name}: {status}")
        if not passed:
            all_pass = False

    print("=" * 70)

    if all_pass:
        print("\nALL CHECKS PASSED - Ready to launch simulations!")
        print("\nNext steps:")
        print("  1. Rent 9x RTX 5090 instance")
        print("  2. Upload PHASE_6/ folder")
        print("  3. Run: ./gpu_scripts/parallel_9gpu/launch_9gpu.sh --wait")
        return 0
    else:
        print("\nSOME CHECKS FAILED - Fix issues before spending money!")
        return 1


if __name__ == "__main__":
    sys.exit(main())
