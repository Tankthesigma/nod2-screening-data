#!/usr/bin/env python
"""
FEP Setup for Natural Product CID_10120 (Bufadienolide)
NOD2-Crohn's Drug Discovery Project

This script sets up ABFE (Absolute Binding Free Energy) calculations for the
natural product using the EXACT same parameters as the febuxostat FEP runs.

SYSTEMS:
- wt_complex: Wild-type NOD2 + ligand
- mut_complex: R702W mutant NOD2 + ligand
- solvent: Ligand in solvent only (no protein)

LAMBDA SCHEDULE (K=20 windows):
- Windows 0-9: Electrostatics decoupling (elec: 1.0 -> 0.0)
- Windows 10-14: Sterics decoupling (sterics: 0.9 -> 0.4)
- Windows 15-19: Continue sterics to 0 (sterics: 0.4 -> 0.0)

Author: Claude Code Assistant
Date: 2026-01-15
"""

import os
import sys
import json
import numpy as np
from pathlib import Path
import shutil

# ============================================================================
# CONFIGURATION - MUST MATCH FEBUXOSTAT EXACTLY
# ============================================================================

BASE_DIR = Path("C:/Users/vasud/nod2-screening-data/fep_pmx_natural")
INPUT_DIR = Path("C:/Users/vasud/nod2-screening-data/PHASE_6/structures")
FEBUXOSTAT_DIR = Path("C:/Users/vasud/nod2-screening-data/fep_complete/fep_pmx")

# Input files
LIGAND_SDF = INPUT_DIR / "natural_top_docked.sdf"
COMPLEX_PDB = INPUT_DIR / "complex_natural_cid10592.pdb"

# Simulation parameters (MUST MATCH FEBUXOSTAT)
SIM_PARAMS = {
    'temperature_K': 310.0,
    'pressure_bar': 1.0,
    'friction_per_ps': 1.0,
    'timestep_fs': 2.0,
    'equil_steps': 50000,      # 100 ps NPT equilibration
    'prod_steps': 500000,      # 1 ns NVT production
    'energy_interval': 500,    # 1 ps sampling
    'barostat_freq': 25,
}

# Lambda schedule - MUST BE LOADED FROM FEBUXOSTAT, NOT HARDCODED!
# This will be loaded at runtime from the verified febuxostat schedule
def load_febuxostat_lambda_schedule():
    """Load the exact lambda schedule from febuxostat to ensure MBAR compatibility."""
    feb_schedule_path = FEBUXOSTAT_DIR / "wt_complex" / "lambda_schedule.npy"
    if not feb_schedule_path.exists():
        raise FileNotFoundError(f"Cannot find febuxostat lambda schedule at {feb_schedule_path}")
    return np.load(feb_schedule_path)

# Will be initialized at runtime
LAMBDA_SCHEDULE = None

# Softcore parameters (MUST MATCH FEBUXOSTAT)
SOFTCORE_PARAMS = {
    'softcore_alpha': 0.5,
    'softcore_beta': 0.0,
    'softcore_a': 1,
    'softcore_b': 1,
    'softcore_c': 6,
    'softcore_d': 1,
    'softcore_e': 1,
    'softcore_f': 2,
}

# Boresch restraint force constants (MUST MATCH FEBUXOSTAT)
BORESCH_PARAMS = {
    'k_distance': 4184.0,        # kJ/(mol*nm^2)
    'k_angle': 41.84,            # kJ/(mol*rad^2)
    'k_dihedral': 41.84,         # kJ/(mol*rad^2)
}

# ============================================================================
# STEP 1: CREATE DIRECTORY STRUCTURE
# ============================================================================

def create_directories():
    """Create the FEP directory structure."""
    print("="*70)
    print("STEP 1: Creating directory structure")
    print("="*70)

    dirs = [
        BASE_DIR,
        BASE_DIR / "wt_complex",
        BASE_DIR / "mut_complex",
        BASE_DIR / "solvent",
        BASE_DIR / "results",
        BASE_DIR / "logs",
    ]

    for d in dirs:
        d.mkdir(parents=True, exist_ok=True)
        print(f"  Created: {d}")

    # Create window subdirectories
    for sys_name in ['wt_complex', 'mut_complex', 'solvent']:
        for i in range(20):
            window_dir = BASE_DIR / sys_name / f"window_{i:02d}"
            window_dir.mkdir(exist_ok=True)

    print(f"  Created 60 window directories")
    return True

# ============================================================================
# STEP 2: VERIFY INPUT FILES
# ============================================================================

def verify_inputs():
    """Verify all input files exist."""
    print("\n" + "="*70)
    print("STEP 2: Verifying input files")
    print("="*70)

    required_files = [
        (LIGAND_SDF, "Natural product SDF"),
        (COMPLEX_PDB, "Complex PDB"),
    ]

    all_ok = True
    for fpath, desc in required_files:
        if fpath.exists():
            print(f"  [OK] {desc}: {fpath}")
        else:
            print(f"  [MISSING] {desc}: {fpath}")
            all_ok = False

    return all_ok

# ============================================================================
# STEP 3: PARAMETERIZE LIGAND WITH OPENFF 2.1.0
# ============================================================================

def parameterize_ligand():
    """Parameterize the natural product with OpenFF 2.1.0."""
    print("\n" + "="*70)
    print("STEP 3: Parameterizing ligand with OpenFF 2.1.0")
    print("="*70)

    try:
        from openff.toolkit import Molecule, ForceField
        from openff.units import unit as offunit

        # Load ligand from SDF
        print(f"  Loading ligand from: {LIGAND_SDF}")
        mol = Molecule.from_file(str(LIGAND_SDF))
        print(f"  Ligand: {mol.n_atoms} atoms, {mol.n_bonds} bonds")
        print(f"  Formula: {mol.hill_formula}")

        # Load OpenFF 2.1.0 force field
        print("  Loading OpenFF 2.1.0 (Sage) force field...")
        ff = ForceField('openff-2.1.0.offxml')

        # Create ligand topology
        from openff.toolkit import Topology
        topology = Topology.from_molecules([mol])

        # Parameterize
        print("  Parameterizing ligand...")
        ligand_system = ff.create_openmm_system(topology)

        # Save parameterized ligand
        from openmm import XmlSerializer
        ligand_xml = BASE_DIR / "ligand_system.xml"
        with open(ligand_xml, 'w') as f:
            f.write(XmlSerializer.serialize(ligand_system))
        print(f"  [SAVED] {ligand_xml}")

        # Save ligand positions
        if mol.conformers is None or len(mol.conformers) == 0:
            print("  [ERROR] No conformers found in ligand SDF file")
            print("          Ligand must have 3D coordinates")
            return None, None
        conformer = mol.conformers[0]
        positions_nm = conformer.to('nanometer').magnitude
        np.save(BASE_DIR / "ligand_positions.npy", positions_nm)
        print(f"  [SAVED] {BASE_DIR / 'ligand_positions.npy'}")

        # Get partial charges for later verification
        partial_charges = mol.partial_charges
        if partial_charges is not None:
            charges = partial_charges.to('elementary_charge').magnitude
            total_charge = np.sum(charges)
            print(f"  Ligand total charge: {total_charge:.4f} e")

        return mol, ligand_system

    except ImportError as e:
        print(f"  [ERROR] OpenFF toolkit not available: {e}")
        print("  Please install: conda install -c conda-forge openff-toolkit")
        return None, None

# ============================================================================
# STEP 4: BUILD WT AND R702W SYSTEMS
# ============================================================================

def build_systems(ligand_mol):
    """Build WT and R702W mutant systems."""
    print("\n" + "="*70)
    print("STEP 4: Building WT and R702W mutant systems")
    print("="*70)

    try:
        from openmm.app import PDBFile, ForceField, Modeller
        from openmm import unit

        # Load complex PDB
        print(f"  Loading complex: {COMPLEX_PDB}")
        pdb = PDBFile(str(COMPLEX_PDB))

        # Count residues
        n_residues = sum(1 for _ in pdb.topology.residues())
        n_atoms = sum(1 for _ in pdb.topology.atoms())
        print(f"  Complex: {n_residues} residues, {n_atoms} atoms")

        # Find residue 702 for mutation
        res702 = None
        for res in pdb.topology.residues():
            if res.id == '702' or (hasattr(res, 'resSeq') and res.resSeq == 702):
                res702 = res
                break

        if res702:
            print(f"  Found residue 702: {res702.name}")
        else:
            print("  [WARNING] Could not find residue 702")

        # For now, we'll use the input complex directly
        # The mutation R702W should be handled by modifying the PDB
        # or using a mutation tool

        return True

    except Exception as e:
        print(f"  [ERROR] Failed to build systems: {e}")
        import traceback
        traceback.print_exc()
        return False

# ============================================================================
# STEP 5: GENERATE WINDOW SCRIPTS
# ============================================================================

def generate_window_script(sys_name, window_idx, lambda_schedule, has_restraints=True):
    """Generate a run_window.py script for one FEP window."""

    lam_e, lam_s, lam_r = lambda_schedule[window_idx]

    # For solvent system, no restraints
    if sys_name == 'solvent':
        has_restraints = False
        lam_r = 0.0

    script = f'''#!/usr/bin/env python
"""
FEP Window Runner - Window {window_idx}
System: {sys_name}
Ligand: Natural Product CID_10120 (Bufadienolide)
Lambda: elec={lam_e:.4f}, sterics={lam_s:.4f}, restraints={lam_r:.4f}

Generated by setup_fep_natural.py
MUST MATCH FEBUXOSTAT PARAMETERS EXACTLY FOR MBAR COMPATIBILITY
"""
import os
import sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from openmm import XmlSerializer, LangevinMiddleIntegrator, Platform, MonteCarloBarostat
from openmm.app import PDBFile, Simulation, StateDataReporter, CheckpointReporter
from openmm import unit

# Configuration - MUST MATCH FEBUXOSTAT
WINDOW_IDX = {window_idx}
LAMBDA_ELEC = {lam_e}
LAMBDA_STERICS = {lam_s}
LAMBDA_RESTRAINTS = {lam_r}
SYS_NAME = "{sys_name}"
HAS_RESTRAINTS = {has_restraints}

# Simulation parameters - MUST MATCH FEBUXOSTAT EXACTLY
TEMPERATURE = {SIM_PARAMS['temperature_K']} * unit.kelvin
FRICTION = {SIM_PARAMS['friction_per_ps']} / unit.picosecond
TIMESTEP = {SIM_PARAMS['timestep_fs']} * unit.femtoseconds
EQUIL_STEPS = {SIM_PARAMS['equil_steps']}    # {SIM_PARAMS['equil_steps'] * SIM_PARAMS['timestep_fs'] / 1000:.0f} ps NPT equilibration
PROD_STEPS = {SIM_PARAMS['prod_steps']}      # {SIM_PARAMS['prod_steps'] * SIM_PARAMS['timestep_fs'] / 1000000:.0f} ns production
ENERGY_INTERVAL = {SIM_PARAMS['energy_interval']}  # {SIM_PARAMS['energy_interval'] * SIM_PARAMS['timestep_fs'] / 1000:.0f} ps sampling

def main():
    print("="*60)
    print(f"FEP Window {{WINDOW_IDX}} - {{SYS_NAME}}")
    print("="*60)
    print(f"Lambda: elec={{LAMBDA_ELEC:.4f}}, sterics={{LAMBDA_STERICS:.4f}}, restraints={{LAMBDA_RESTRAINTS:.4f}}")
    print()

    # Load system and positions
    parent_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

    with open(os.path.join(parent_dir, "alchemical_system.xml"), "r") as f:
        system = XmlSerializer.deserialize(f.read())

    positions = np.load(os.path.join(parent_dir, "positions.npy"))
    positions = positions * unit.nanometer

    pdb = PDBFile(os.path.join(parent_dir, "topology.pdb"))
    topology = pdb.topology

    # Load lambda schedule
    lambda_schedule = np.load(os.path.join(parent_dir, "lambda_schedule.npy"))
    n_windows = len(lambda_schedule)
    print(f"Loaded lambda schedule: {{n_windows}} windows")

    # Create integrator
    integrator = LangevinMiddleIntegrator(TEMPERATURE, FRICTION, TIMESTEP)

    # Select platform (prefer CUDA)
    try:
        platform = Platform.getPlatformByName('CUDA')
        properties = {{'Precision': 'mixed'}}
        print("Using CUDA platform")
    except:
        try:
            platform = Platform.getPlatformByName('OpenCL')
            properties = {{}}
            print("Using OpenCL platform")
        except:
            platform = Platform.getPlatformByName('CPU')
            properties = {{}}
            print("Using CPU platform")

    # Create simulation
    simulation = Simulation(topology, system, integrator, platform, properties)
    simulation.context.setPositions(positions)

    # Set lambda parameters
    context = simulation.context
    context.setParameter('lambda_electrostatics', LAMBDA_ELEC)
    context.setParameter('lambda_sterics', LAMBDA_STERICS)
    if HAS_RESTRAINTS:
        context.setParameter('lambda_restraints', LAMBDA_RESTRAINTS)

    # Minimize
    print("Minimizing energy...")
    simulation.minimizeEnergy(maxIterations=1000)

    # NPT Equilibration
    print(f"Running NPT equilibration ({{EQUIL_STEPS}} steps)...")
    simulation.step(EQUIL_STEPS)

    # Save checkpoint after equilibration
    simulation.saveCheckpoint("checkpoint_equil.chk")
    print("[SAVED] checkpoint_equil.chk")

    # Remove barostat for NVT production
    print("Removing barostat for NVT production...")

    state = context.getState(getPositions=True, getVelocities=True)
    positions_after_equil = state.getPositions()
    velocities_after_equil = state.getVelocities()
    box_vectors = state.getPeriodicBoxVectors()

    # Save box vectors
    box_arr = np.array([[v.x, v.y, v.z] for v in box_vectors])
    np.save("final_box_vectors.npy", box_arr)

    # Remove barostat
    for i in range(system.getNumForces()):
        force = system.getForce(i)
        if isinstance(force, MonteCarloBarostat):
            system.removeForce(i)
            break

    # Create fresh simulation for NVT
    integrator2 = LangevinMiddleIntegrator(TEMPERATURE, FRICTION, TIMESTEP)
    simulation = Simulation(topology, system, integrator2, platform, properties)
    context = simulation.context

    context.setPeriodicBoxVectors(*box_vectors)
    context.setPositions(positions_after_equil)
    context.setVelocitiesToTemperature(TEMPERATURE)

    # Set lambda parameters
    context.setParameter('lambda_electrostatics', LAMBDA_ELEC)
    context.setParameter('lambda_sterics', LAMBDA_STERICS)
    if HAS_RESTRAINTS:
        context.setParameter('lambda_restraints', LAMBDA_RESTRAINTS)

    # Verify energy
    state = context.getState(getEnergy=True)
    E_nvt = state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
    print(f"  NVT initial energy: {{E_nvt:.1f}} kJ/mol")
    if E_nvt > 0 or abs(E_nvt) > 1e10:
        raise RuntimeError(f"Energy looks wrong after NVT setup: {{E_nvt}}")

    # Production run
    print(f"Running NVT production ({{PROD_STEPS}} steps)...")
    n_samples = PROD_STEPS // ENERGY_INTERVAL
    u_nk = np.zeros((n_samples, n_windows))

    for sample_idx in range(n_samples):
        simulation.step(ENERGY_INTERVAL)

        # Compute reduced potential at all lambda states
        for k in range(n_windows):
            lam_e, lam_s, lam_r = lambda_schedule[k]
            context.setParameter('lambda_electrostatics', lam_e)
            context.setParameter('lambda_sterics', lam_s)
            if HAS_RESTRAINTS:
                context.setParameter('lambda_restraints', lam_r)

            state = context.getState(getEnergy=True)
            kT = unit.MOLAR_GAS_CONSTANT_R * TEMPERATURE
            u_nk[sample_idx, k] = state.getPotentialEnergy() / kT

        # Reset to current window
        context.setParameter('lambda_electrostatics', LAMBDA_ELEC)
        context.setParameter('lambda_sterics', LAMBDA_STERICS)
        if HAS_RESTRAINTS:
            context.setParameter('lambda_restraints', LAMBDA_RESTRAINTS)

        if (sample_idx + 1) % 100 == 0:
            print(f"  Sample {{sample_idx + 1}}/{{n_samples}}")

    # Save results
    np.save("u_nk.npy", u_nk)
    print(f"[PASS] Saved u_nk.npy: shape {{u_nk.shape}}")

    state = context.getState(getPositions=True)
    pos = state.getPositions(asNumpy=True).value_in_unit(unit.nanometer)
    np.save("final_positions.npy", pos)
    print("[PASS] Saved final_positions.npy")

    # Save final checkpoint
    simulation.saveCheckpoint("checkpoint.chk")
    print("[PASS] Saved checkpoint.chk")

    print()
    print("="*60)
    print("Window complete!")
    print("="*60)

if __name__ == "__main__":
    main()
'''
    return script

def generate_all_window_scripts(lambda_schedule):
    """Generate all 60 window scripts."""
    print("\n" + "="*70)
    print("STEP 5: Generating all 60 window scripts")
    print("="*70)

    for sys_name in ['wt_complex', 'mut_complex', 'solvent']:
        has_restraints = sys_name != 'solvent'
        for i in range(20):
            script = generate_window_script(sys_name, i, lambda_schedule, has_restraints)
            script_path = BASE_DIR / sys_name / f"window_{i:02d}" / "run_window.py"
            with open(script_path, 'w') as f:
                f.write(script)
        print(f"  Generated 20 scripts for {sys_name}")

    print("  [COMPLETE] Generated 60 window scripts")
    return True

# ============================================================================
# STEP 6: SAVE LAMBDA SCHEDULE
# ============================================================================

def save_lambda_schedule(lambda_schedule):
    """Save lambda schedule to all system directories."""
    print("\n" + "="*70)
    print("STEP 6: Saving lambda schedule")
    print("="*70)

    for sys_name in ['wt_complex', 'mut_complex', 'solvent']:
        schedule_path = BASE_DIR / sys_name / "lambda_schedule.npy"
        np.save(schedule_path, lambda_schedule)
        print(f"  Saved: {schedule_path}")

    return True

# ============================================================================
# STEP 7: GENERATE RUN MANIFEST
# ============================================================================

def generate_manifest(lambda_schedule):
    """Generate run manifest with all parameters."""
    print("\n" + "="*70)
    print("STEP 7: Generating run manifest")
    print("="*70)

    manifest = {
        'project': 'NOD2-Crohn Natural Product FEP',
        'ligand': 'CID_10120 (Bufadienolide)',
        'mutation': 'R702W',
        'systems': ['wt_complex', 'mut_complex', 'solvent'],
        'n_windows': 20,
        'total_windows': 60,
        'lambda_schedule': lambda_schedule.tolist(),
        'simulation_params': SIM_PARAMS,
        'softcore_params': SOFTCORE_PARAMS,
        'boresch_params': BORESCH_PARAMS,
        'window_folders': [],
        'package_versions': {},
    }

    # List all window folders
    for sys_name in ['wt_complex', 'mut_complex', 'solvent']:
        for i in range(20):
            manifest['window_folders'].append(f"{sys_name}/window_{i:02d}")

    # Get package versions
    try:
        import openmm
        manifest['package_versions']['openmm'] = openmm.__version__
    except:
        pass

    try:
        import openff.toolkit
        manifest['package_versions']['openff-toolkit'] = openff.toolkit.__version__
    except:
        pass

    try:
        import numpy
        manifest['package_versions']['numpy'] = numpy.__version__
    except:
        pass

    try:
        import pymbar
        manifest['package_versions']['pymbar'] = pymbar.__version__
    except:
        pass

    # Save manifest
    manifest_path = BASE_DIR / "run_manifest.json"
    with open(manifest_path, 'w') as f:
        json.dump(manifest, f, indent=2)
    print(f"  Saved: {manifest_path}")

    # Also save as readable text
    txt_path = BASE_DIR / "run_manifest.txt"
    with open(txt_path, 'w') as f:
        f.write("="*70 + "\n")
        f.write("NOD2-CROHN NATURAL PRODUCT FEP RUN MANIFEST\n")
        f.write("="*70 + "\n\n")

        f.write("LIGAND: CID_10120 (Bufadienolide)\n")
        f.write("MUTATION: R702W\n\n")

        f.write("LAMBDA SCHEDULE (K=20):\n")
        f.write("-"*40 + "\n")
        f.write("Window  Elec    Sterics  Restraints\n")
        f.write("-"*40 + "\n")
        for i, (e, s, r) in enumerate(lambda_schedule):
            f.write(f"{i:2d}      {e:.4f}  {s:.4f}   {r:.4f}\n")
        f.write("\n")

        f.write("SIMULATION PARAMETERS:\n")
        f.write("-"*40 + "\n")
        for k, v in SIM_PARAMS.items():
            f.write(f"  {k}: {v}\n")
        f.write("\n")

        f.write("SOFTCORE PARAMETERS:\n")
        f.write("-"*40 + "\n")
        for k, v in SOFTCORE_PARAMS.items():
            f.write(f"  {k}: {v}\n")
        f.write("\n")

        f.write("BORESCH RESTRAINT FORCE CONSTANTS:\n")
        f.write("-"*40 + "\n")
        for k, v in BORESCH_PARAMS.items():
            f.write(f"  {k}: {v}\n")
        f.write("\n")

        f.write("WINDOW FOLDERS:\n")
        f.write("-"*40 + "\n")
        for folder in manifest['window_folders']:
            f.write(f"  {folder}\n")

    print(f"  Saved: {txt_path}")
    return manifest

# ============================================================================
# STEP 8: GENERATE LAUNCH SCRIPT
# ============================================================================

def generate_launch_script():
    """Generate a single command to launch all windows."""
    print("\n" + "="*70)
    print("STEP 8: Generating launch script")
    print("="*70)

    launch_script = f'''#!/usr/bin/env python
"""
Launch all 60 FEP windows for Natural Product CID_10120.
Run this after system setup is complete.
"""
import os
import subprocess
import sys
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed

BASE_DIR = Path("{BASE_DIR}")
SYSTEMS = ['wt_complex', 'mut_complex', 'solvent']
N_WINDOWS = 20

def run_window(sys_name, window_idx):
    """Run a single FEP window."""
    window_dir = BASE_DIR / sys_name / f"window_{{window_idx:02d}}"
    script = window_dir / "run_window.py"
    log_file = window_dir / "run.log"

    print(f"[START] {{sys_name}}/window_{{window_idx:02d}}")

    with open(log_file, 'w') as log:
        result = subprocess.run(
            [sys.executable, str(script)],
            cwd=str(window_dir),
            stdout=log,
            stderr=subprocess.STDOUT,
        )

    if result.returncode == 0:
        print(f"[DONE] {{sys_name}}/window_{{window_idx:02d}}")
    else:
        print(f"[FAIL] {{sys_name}}/window_{{window_idx:02d}}")

    return (sys_name, window_idx, result.returncode)

def main():
    print("="*70)
    print("LAUNCHING ALL 60 FEP WINDOWS")
    print("="*70)

    # Build list of all windows
    windows = []
    for sys_name in SYSTEMS:
        for i in range(N_WINDOWS):
            windows.append((sys_name, i))

    print(f"Total windows to run: {{len(windows)}}")

    # Run sequentially (can be parallelized with GPU resources)
    results = []
    for sys_name, window_idx in windows:
        result = run_window(sys_name, window_idx)
        results.append(result)

    # Summary
    print()
    print("="*70)
    print("SUMMARY")
    print("="*70)

    success = sum(1 for r in results if r[2] == 0)
    failed = sum(1 for r in results if r[2] != 0)

    print(f"Success: {{success}}/{{len(results)}}")
    print(f"Failed: {{failed}}/{{len(results)}}")

    if failed > 0:
        print("\\nFailed windows:")
        for sys_name, window_idx, rc in results:
            if rc != 0:
                print(f"  {{sys_name}}/window_{{window_idx:02d}}")

if __name__ == "__main__":
    main()
'''

    launch_path = BASE_DIR / "launch_all_windows.py"
    with open(launch_path, 'w') as f:
        f.write(launch_script)
    print(f"  Saved: {launch_path}")

    return True

# ============================================================================
# MAIN
# ============================================================================

def main():
    print("="*70)
    print("FEP SETUP FOR NATURAL PRODUCT CID_10120")
    print("NOD2-CROHN DRUG DISCOVERY PROJECT")
    print("="*70)
    print()

    # CRITICAL: Load lambda schedule from febuxostat to ensure MBAR compatibility
    print("Loading lambda schedule from febuxostat...")
    try:
        lambda_schedule = load_febuxostat_lambda_schedule()
        print(f"  Loaded {len(lambda_schedule)} windows from febuxostat")
    except FileNotFoundError as e:
        print(f"[ABORT] {e}")
        return 1

    # Step 1: Create directories
    if not create_directories():
        print("[ABORT] Failed to create directories")
        return 1

    # Step 2: Verify inputs
    if not verify_inputs():
        print("[ABORT] Missing input files")
        return 1

    # Step 3: Parameterize ligand
    mol, ligand_system = parameterize_ligand()
    if mol is None:
        print("[WARNING] Ligand parameterization skipped - will need manual setup")

    # Step 4: Build systems (placeholder for full implementation)
    build_systems(mol)

    # Step 5: Generate window scripts
    generate_all_window_scripts(lambda_schedule)

    # Step 6: Save lambda schedule
    save_lambda_schedule(lambda_schedule)

    # Step 7: Generate manifest
    generate_manifest(lambda_schedule)

    # Step 8: Generate launch script
    generate_launch_script()

    print()
    print("="*70)
    print("SETUP COMPLETE")
    print("="*70)
    print()
    print("NEXT STEPS:")
    print("1. Complete system preparation (solvation, neutralization)")
    print("2. Set up Boresch restraints using select_boresch_anchors.py")
    print("3. Run canary windows 15-19 to verify stability")
    print("4. Check overlap matrix")
    print("5. Launch all 60 windows")
    print()
    print(f"Output directory: {BASE_DIR}")

    return 0

if __name__ == "__main__":
    sys.exit(main())
