#!/usr/bin/env python3
"""
LOCAL RTX 4060 Ti - Phase A1 Mutant MD Simulations (4 sims, 80ns total)
Sequential execution, NO shutdown.

Simulations:
1. R702W + Febuxostat rep1 (20ns)
2. R702W + Febuxostat rep2 (20ns)
3. R702W + Febuxostat rep3 (20ns)
4. G908R + Febuxostat rep1 (20ns)

Expected runtime: ~8-10 hours on RTX 4060 Ti with HMR (~200 ns/day)
"""

from openmm import *
from openmm.app import *
from openmm.unit import *
import sys
import os
import gc
import math
from datetime import datetime
import traceback

try:
    from openmmforcefields.generators import SystemGenerator, GAFFTemplateGenerator
    from openff.toolkit.topology import Molecule
    OPENFF_AVAILABLE = True
except ImportError:
    print("FATAL: OpenFF/openmmforcefields not installed!")
    print("Run: conda install -c conda-forge openmm openmmforcefields openff-toolkit -y")
    sys.exit(1)

try:
    from pdbfixer import PDBFixer
    PDBFIXER_AVAILABLE = True
except ImportError:
    print("WARNING: PDBFixer not installed - will skip structure fixing")
    print("Run: conda install -c conda-forge pdbfixer -y")
    PDBFIXER_AVAILABLE = False

# ============================================================
# CONFIGURATION
# ============================================================
GPU_ID = 0

STRUCTURES_DIR = "structures"
LIGANDS_DIR = "ligands"
OUTPUT_DIR = "trajectories"
LOG_DIR = "logs"
CHECKPOINT_DIR = "checkpoints"

TEMPERATURE = 310.15 * kelvin
PRESSURE = 1.0 * atmospheres
TIMESTEP = 4.0 * femtoseconds
FRICTION = 1.0 / picoseconds
USE_HMR = True
HYDROGEN_MASS = 3.0 * amu

MINIMIZATION_STEPS = 5000
EQUILIBRATION_NVT_PS = 100
EQUILIBRATION_NPT_PS = 500
PRODUCTION_NS = 20

REPORT_INTERVAL = 2500
TRAJECTORY_INTERVAL = 2500
CHECKPOINT_INTERVAL = 250000

# LOCAL 4060 SET: R702W+Febuxostat x3, G908R+Febuxostat x1
SIMULATIONS = [
    ("R702W_febuxostat", 1, "NOD2_R702W.pdb", "febuxostat_docked.sdf"),
    ("R702W_febuxostat", 2, "NOD2_R702W.pdb", "febuxostat_docked.sdf"),
    ("R702W_febuxostat", 3, "NOD2_R702W.pdb", "febuxostat_docked.sdf"),
    ("G908R_febuxostat", 1, "NOD2_G908R.pdb", "febuxostat_docked.sdf"),
]

# ============================================================
# HELPER FUNCTIONS
# ============================================================
def setup_directories():
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    os.makedirs(LOG_DIR, exist_ok=True)
    os.makedirs(CHECKPOINT_DIR, exist_ok=True)

def fix_pdb_structure(pdb_path):
    """Fix PDB structure using PDBFixer - adds missing atoms/residues."""
    if not PDBFIXER_AVAILABLE:
        print("      PDBFixer not available, loading as-is")
        return PDBFile(pdb_path)

    print("      Running PDBFixer...")
    fixer = PDBFixer(filename=pdb_path)
    fixer.findMissingResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(7.0)  # pH 7.0

    print(f"      Fixed: {fixer.topology.getNumAtoms()} atoms")

    class FixedPDB:
        def __init__(self, topology, positions):
            self.topology = topology
            self.positions = positions

    return FixedPDB(fixer.topology, fixer.positions)

def check_energy(simulation, stage=""):
    """Check if energy is valid (not NaN/Inf). Returns False if exploded."""
    state = simulation.context.getState(getEnergy=True)
    energy = state.getPotentialEnergy()
    energy_val = energy.value_in_unit(kilojoules_per_mole)

    if math.isnan(energy_val) or math.isinf(energy_val):
        print(f"\n!!! CRITICAL: Energy exploded during {stage}!")
        print(f"!!! Energy = {energy}")
        print("!!! Check input PDB/SDF for clashes or bad geometry")
        return False, energy
    return True, energy

def create_system_with_fallback(pdb, ligand_mol, forcefield_kwargs):
    """Create system with OpenFF, falling back to GAFF if needed."""
    modeller = Modeller(pdb.topology, pdb.positions)
    modeller.deleteWater()
    print("      Deleted existing water molecules")

    # ADD LIGAND TO MODELLER
    if ligand_mol is not None:
        print("      Adding ligand to system...")
        ligand_topology = ligand_mol.to_topology().to_openmm()
        ligand_positions = ligand_mol.conformers[0].to_openmm()
        modeller.add(ligand_topology, ligand_positions)
        print(f"      Ligand added: {ligand_mol.n_atoms} atoms")

    if OPENFF_AVAILABLE:
        try:
            print("      Trying OpenFF...")
            system_generator = SystemGenerator(
                forcefields=['amber14-all.xml', 'amber14/tip3pfb.xml'],
                small_molecule_forcefield='openff-2.1.0',
                molecules=[ligand_mol],
                forcefield_kwargs=forcefield_kwargs,
                periodic_forcefield_kwargs={
                    'nonbondedMethod': PME,
                    'nonbondedCutoff': 1.0*nanometers
                }
            )
            modeller.addSolvent(system_generator.forcefield, model='tip3p',
                              padding=1.0*nanometers, ionicStrength=0.15*molar)
            system = system_generator.create_system(modeller.topology)
            for force in system.getForces():
                if isinstance(force, NonbondedForce):
                    force.setUseDispersionCorrection(True)
            print("      SUCCESS: OpenFF (with HMR + dispersion correction)")
            return system, modeller
        except Exception as e:
            print(f"      OpenFF failed: {e}, trying GAFF...")
            modeller = Modeller(pdb.topology, pdb.positions)
            modeller.deleteWater()
            if ligand_mol is not None:
                ligand_topology = ligand_mol.to_topology().to_openmm()
                ligand_positions = ligand_mol.conformers[0].to_openmm()
                modeller.add(ligand_topology, ligand_positions)

    print("      Using GAFF...")
    forcefield = ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')
    gaff_generator = GAFFTemplateGenerator(molecules=[ligand_mol])
    forcefield.registerTemplateGenerator(gaff_generator.generator)
    modeller.addSolvent(forcefield, model='tip3p',
                      padding=1.0*nanometers, ionicStrength=0.15*molar)
    system = forcefield.createSystem(
        modeller.topology, nonbondedMethod=PME, nonbondedCutoff=1.0*nanometers,
        constraints=HBonds, hydrogenMass=HYDROGEN_MASS if USE_HMR else None
    )
    for force in system.getForces():
        if isinstance(force, NonbondedForce):
            force.setUseDispersionCorrection(True)
    print("      SUCCESS: GAFF (with HMR + dispersion correction)")
    return system, modeller

# ============================================================
# SIMULATION FUNCTION
# ============================================================
def run_simulation(name, replicate, pdb_file, sdf_file, sim_num, total_sims):
    """Run a single MD simulation with full error checking."""
    # Seeds: 100000 base for local script, well separated from Vast scripts
    random_seed = 100000 + sim_num * 1000 + replicate
    output_prefix = f"{name}_rep{replicate}"

    pdb_path = f"{STRUCTURES_DIR}/{pdb_file}"
    sdf_path = f"{LIGANDS_DIR}/{sdf_file}"

    print("\n" + "=" * 70)
    print(f"[{sim_num}/{total_sims}] {name.upper()} - Replicate {replicate}")
    print(f"Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("=" * 70)

    # Validate inputs
    if not os.path.exists(pdb_path):
        print(f"ERROR: PDB not found: {pdb_path}")
        return False
    if not os.path.exists(sdf_path):
        print(f"ERROR: SDF not found: {sdf_path}")
        return False

    # Load ligand
    print("\n[1/8] Loading ligand...")
    ligand_mol = Molecule.from_file(sdf_path)
    if isinstance(ligand_mol, list):
        ligand_mol = ligand_mol[0]
    print(f"      Atoms: {ligand_mol.n_atoms}")

    # Load and fix structure
    print("\n[2/8] Loading and fixing structure...")
    pdb = fix_pdb_structure(pdb_path)
    modeller_h = Modeller(pdb.topology, pdb.positions)
    print(f"      Ready: {modeller_h.topology.getNumAtoms()} atoms")

    # Create system
    print("\n[3/8] Setting up force fields...")
    forcefield_kwargs = {
        'constraints': HBonds,
        'rigidWater': True,
        'removeCMMotion': True,
        'hydrogenMass': HYDROGEN_MASS if USE_HMR else None
    }
    system, modeller = create_system_with_fallback(modeller_h, ligand_mol, forcefield_kwargs)
    print(f"      Solvated atoms: {modeller.topology.getNumAtoms()}")

    # Save solvated topology
    solvated_pdb = f"{OUTPUT_DIR}/{output_prefix}_solvated.pdb"
    with open(solvated_pdb, 'w') as f:
        PDBFile.writeFile(modeller.topology, modeller.positions, f)
    print(f"      Saved: {solvated_pdb}")

    # Setup simulation
    print(f"\n[4/8] Setting up CUDA (GPU {GPU_ID})...")
    integrator = LangevinMiddleIntegrator(TEMPERATURE, FRICTION, TIMESTEP)
    integrator.setRandomNumberSeed(random_seed)
    integrator.setConstraintTolerance(1e-6)
    try:
        platform = Platform.getPlatformByName('CUDA')
        properties = {'Precision': 'mixed', 'DeviceIndex': str(GPU_ID)}
        simulation = Simulation(modeller.topology, system, integrator, platform, properties)
    except Exception as e:
        print(f"      CUDA failed: {e}, trying OpenCL...")
        platform = Platform.getPlatformByName('OpenCL')
        properties = {'Precision': 'mixed', 'OpenCLPlatformIndex': '0', 'DeviceIndex': str(GPU_ID)}
        simulation = Simulation(modeller.topology, system, integrator, platform, properties)
    simulation.context.setPositions(modeller.positions)

    # Minimization
    print("\n[5/8] MINIMIZATION")
    ok, energy = check_energy(simulation, "pre-minimization")
    print(f"      Initial: {energy}")

    simulation.minimizeEnergy(maxIterations=MINIMIZATION_STEPS)

    ok, energy = check_energy(simulation, "post-minimization")
    if not ok:
        return False
    print(f"      Final: {energy}")

    # NVT with GRADUAL HEATING
    print(f"\n[6/8] NVT EQUILIBRATION ({EQUILIBRATION_NVT_PS} ps) - GRADUAL HEATING")

    heating_stages = [
        (50*kelvin, 5),
        (100*kelvin, 10),
        (150*kelvin, 10),
        (200*kelvin, 15),
        (250*kelvin, 20),
        (TEMPERATURE, 40),
    ]

    for temp, duration_ps in heating_stages:
        print(f"      Heating: {duration_ps}ps at {temp}")
        simulation.context.setVelocitiesToTemperature(temp, random_seed)
        integrator.setTemperature(temp)
        steps = int(duration_ps * 1000 / TIMESTEP.value_in_unit(femtoseconds))
        simulation.step(steps)

        ok, energy = check_energy(simulation, f"NVT-{temp}")
        if not ok:
            return False

    print("      NVT complete")

    # NPT
    print(f"\n[7/8] NPT EQUILIBRATION ({EQUILIBRATION_NPT_PS} ps)")
    system.addForce(MonteCarloBarostat(PRESSURE, TEMPERATURE, 100))
    simulation.context.reinitialize(preserveState=True)
    npt_steps = int(EQUILIBRATION_NPT_PS * 1000 / TIMESTEP.value_in_unit(femtoseconds))
    simulation.step(npt_steps)

    ok, energy = check_energy(simulation, "NPT")
    if not ok:
        return False

    # Production
    print(f"\n[8/8] PRODUCTION ({PRODUCTION_NS} ns)")
    production_steps = int(PRODUCTION_NS * 1e6 / TIMESTEP.value_in_unit(femtoseconds))

    simulation.reporters.append(DCDReporter(
        f'{OUTPUT_DIR}/{output_prefix}.dcd', TRAJECTORY_INTERVAL
    ))
    simulation.reporters.append(StateDataReporter(
        f'{LOG_DIR}/{output_prefix}.log', REPORT_INTERVAL,
        step=True, time=True, potentialEnergy=True, kineticEnergy=True,
        temperature=True, speed=True, progress=True, remainingTime=True,
        totalSteps=production_steps
    ))
    simulation.reporters.append(CheckpointReporter(
        f'{CHECKPOINT_DIR}/{output_prefix}.chk', CHECKPOINT_INTERVAL
    ))

    simulation.step(production_steps)
    simulation.saveState(f'{CHECKPOINT_DIR}/{output_prefix}_final.xml')

    ok, energy = check_energy(simulation, "production")
    if not ok:
        print("WARNING: Final energy check failed but trajectory saved")

    print(f"\n>>> COMPLETE: {output_prefix}")
    print(f"    Finished: {datetime.now().strftime('%H:%M:%S')}")

    del simulation, system, integrator, modeller
    gc.collect()

    return True

# ============================================================
# MAIN
# ============================================================
def preflight_check():
    """Check ALL input files exist before starting."""
    print("\n[PREFLIGHT] Checking input files...")
    missing = []
    for name, rep, pdb, sdf in SIMULATIONS:
        pdb_path = f"{STRUCTURES_DIR}/{pdb}"
        sdf_path = f"{LIGANDS_DIR}/{sdf}"
        if not os.path.exists(pdb_path):
            missing.append(pdb_path)
        if not os.path.exists(sdf_path):
            missing.append(sdf_path)

    if missing:
        print("\n!!! FATAL: Missing input files !!!")
        for f in missing:
            print(f"    - {f}")
        sys.exit(1)

    print(f"    All input files found.")

    # Validate ligands
    print("\n[PREFLIGHT] Validating ligand files...")
    unique_sdfs = set(sdf for _, _, _, sdf in SIMULATIONS)
    for sdf in unique_sdfs:
        sdf_path = f"{LIGANDS_DIR}/{sdf}"
        try:
            mol = Molecule.from_file(sdf_path)
            if isinstance(mol, list):
                mol = mol[0]
            print(f"    OK {sdf}: {mol.n_atoms} atoms, {mol.n_bonds} bonds")
        except Exception as e:
            print(f"\n!!! FATAL: Cannot load ligand {sdf} !!!")
            print(f"    Error: {e}")
            sys.exit(1)

    print("    All ligands validated successfully.")

def main():
    print("=" * 70)
    print("LOCAL RTX 4060 Ti - PHASE A1 MUTANT MD SIMULATIONS")
    print(f"Running {len(SIMULATIONS)} simulations ({len(SIMULATIONS) * PRODUCTION_NS}ns total)")
    print(f"Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("=" * 70)

    preflight_check()
    setup_directories()

    results = []
    for i, (name, rep, pdb, sdf) in enumerate(SIMULATIONS, 1):
        try:
            success = run_simulation(name, rep, pdb, sdf, i, len(SIMULATIONS))
            results.append((f"{name}_rep{rep}", success))
        except Exception as e:
            print(f"\nERROR in {name}_rep{rep}: {e}")
            traceback.print_exc()
            results.append((f"{name}_rep{rep}", False))
            print("\nContinuing to next simulation...")
            gc.collect()

    # Summary
    print("\n" + "=" * 70)
    print("FINAL SUMMARY")
    print("=" * 70)
    success_count = len([s for _, s in results if s])
    for name, success in results:
        status = "COMPLETE" if success else "FAILED"
        print(f"  {name}: {status}")
    print(f"\nTotal: {success_count}/{len(results)} successful")
    print(f"Finished: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("=" * 70)

if __name__ == "__main__":
    script_dir = os.path.dirname(os.path.abspath(__file__))
    phase_a1_dir = os.path.join(script_dir, "..")
    os.chdir(phase_a1_dir)
    print(f"Working directory: {os.getcwd()}")

    try:
        main()
    except Exception as e:
        print(f"\n!!! FATAL ERROR: {e}")
        traceback.print_exc()
        sys.exit(1)
