#!/usr/bin/env python3
"""
VAST.AI RTX 5090 - Heavy Load (10 simulations, 200ns total)
Sequential execution with auto-shutdown on completion.

Simulations:
1. Febuxostat rep1-3 (60ns)
2. Ursodiol rep1-3 (60ns)
3. Budesonide rep1-3 (60ns)
4. Natural_top rep1 (20ns)

Expected runtime: ~5 hours on RTX 5090 with HMR
"""

from openmm import *
from openmm.app import *
from openmm.unit import *
import sys
import os
import gc
from datetime import datetime
import traceback

try:
    from openmmforcefields.generators import SystemGenerator, GAFFTemplateGenerator
    from openff.toolkit import Molecule
    OPENFF_AVAILABLE = True
except ImportError:
    OPENFF_AVAILABLE = False
    print("WARNING: OpenFF not available, using GAFF only")

# ============================================================
# CONFIGURATION
# ============================================================
GPU_ID = 0  # Single GPU on rented instance

# Directories (relative to PHASE_6/)
STRUCTURES_DIR = "structures"
OUTPUT_DIR = "trajectories"
LOG_DIR = "logs"
CHECKPOINT_DIR = "checkpoints"

# Physics
TEMPERATURE = 310.15 * kelvin
PRESSURE = 1.0 * atmospheres
TIMESTEP = 4.0 * femtoseconds
FRICTION = 1.0 / picoseconds
USE_HMR = True
HYDROGEN_MASS = 3.0 * amu  # Proper HMR for 4fs

# Protocol
MINIMIZATION_STEPS = 5000
EQUILIBRATION_NVT_PS = 100
EQUILIBRATION_NPT_PS = 500
PRODUCTION_NS = 20

# Output intervals
REPORT_INTERVAL = 2500       # 10ps
TRAJECTORY_INTERVAL = 2500   # 10ps
CHECKPOINT_INTERVAL = 250000 # 1ns

# ============================================================
# SIMULATIONS TO RUN (10 total = 200ns)
# ============================================================
SIMULATIONS = [
    # (name, replicate, pdb_file, sdf_file)
    ("febuxostat", 1, "complex_febuxostat.pdb", "febuxostat_docked.sdf"),
    ("febuxostat", 2, "complex_febuxostat.pdb", "febuxostat_docked.sdf"),
    ("febuxostat", 3, "complex_febuxostat.pdb", "febuxostat_docked.sdf"),
    ("ursodiol", 1, "complex_ursodiol.pdb", "ursodiol_docked.sdf"),
    ("ursodiol", 2, "complex_ursodiol.pdb", "ursodiol_docked.sdf"),
    ("ursodiol", 3, "complex_ursodiol.pdb", "ursodiol_docked.sdf"),
    ("budesonide", 1, "complex_budesonide.pdb", "budesonide_docked.sdf"),
    ("budesonide", 2, "complex_budesonide.pdb", "budesonide_docked.sdf"),
    ("budesonide", 3, "complex_budesonide.pdb", "budesonide_docked.sdf"),
    ("natural_top", 1, "complex_natural_cid10592.pdb", "natural_top_docked.sdf"),
]

# ============================================================
# HELPER FUNCTIONS
# ============================================================
def setup_directories():
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    os.makedirs(LOG_DIR, exist_ok=True)
    os.makedirs(CHECKPOINT_DIR, exist_ok=True)

def create_system_with_fallback(pdb, ligand_mol, forcefield_kwargs):
    """Create system with OpenFF, falling back to GAFF if needed."""
    modeller = Modeller(pdb.topology, pdb.positions)

    # Remove existing water to prevent doubling
    modeller.deleteWater()
    print("      Deleted existing water molecules")

    if OPENFF_AVAILABLE:
        try:
            print("      Trying OpenFF...")
            system_generator = SystemGenerator(
                forcefields=['amber14-all.xml', 'amber14/tip3pfb.xml'],
                small_molecule_forcefield='openff-2.1.0',
                molecules=[ligand_mol],
                forcefield_kwargs=forcefield_kwargs
            )
            modeller.addSolvent(system_generator.forcefield, model='tip3p',
                              padding=1.0*nanometers, ionicStrength=0.15*molar)
            system = system_generator.create_system(
                modeller.topology, nonbondedMethod=PME, nonbondedCutoff=1.0*nanometers
            )
            print("      SUCCESS: OpenFF")
            return system, modeller
        except Exception as e:
            print(f"      OpenFF failed: {e}, trying GAFF...")
            modeller = Modeller(pdb.topology, pdb.positions)
            modeller.deleteWater()

    # Fallback to GAFF
    print("      Using GAFF...")
    forcefield = ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')
    gaff_generator = GAFFTemplateGenerator(molecules=[ligand_mol])
    forcefield.registerTemplateGenerator(gaff_generator.generator)
    modeller.addSolvent(forcefield, model='tip3p',
                      padding=1.0*nanometers, ionicStrength=0.15*molar)
    system = forcefield.createSystem(
        modeller.topology, nonbondedMethod=PME, nonbondedCutoff=1.0*nanometers,
        constraints=AllBonds, hydrogenMass=HYDROGEN_MASS if USE_HMR else None
    )
    print("      SUCCESS: GAFF")
    return system, modeller

# ============================================================
# SIMULATION FUNCTION
# ============================================================
def run_simulation(name, replicate, pdb_file, sdf_file, sim_num, total_sims):
    """Run a single MD simulation."""
    random_seed = sim_num * 1000 + replicate
    output_prefix = f"{name}_rep{replicate}"

    pdb_path = f"{STRUCTURES_DIR}/{pdb_file}"
    sdf_path = f"{STRUCTURES_DIR}/{sdf_file}"

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
    print(f"      Atoms: {ligand_mol.n_atoms}")

    # Load structure
    print("\n[2/8] Loading structure...")
    pdb = PDBFile(pdb_path)
    print(f"      Atoms: {pdb.topology.getNumAtoms()}")

    # Create system
    print("\n[3/8] Setting up force fields...")
    forcefield_kwargs = {
        'constraints': AllBonds,
        'rigidWater': True,
        'removeCMMotion': True,
        'hydrogenMass': HYDROGEN_MASS if USE_HMR else None
    }
    system, modeller = create_system_with_fallback(pdb, ligand_mol, forcefield_kwargs)
    print(f"      Solvated atoms: {modeller.topology.getNumAtoms()}")

    # Setup simulation
    print(f"\n[4/8] Setting up CUDA (GPU {GPU_ID})...")
    integrator = LangevinMiddleIntegrator(TEMPERATURE, FRICTION, TIMESTEP)
    integrator.setRandomNumberSeed(random_seed)
    platform = Platform.getPlatformByName('CUDA')
    properties = {'Precision': 'mixed', 'DeviceIndex': str(GPU_ID)}
    simulation = Simulation(modeller.topology, system, integrator, platform, properties)
    simulation.context.setPositions(modeller.positions)

    # Minimization
    print("\n[5/8] MINIMIZATION")
    state = simulation.context.getState(getEnergy=True)
    print(f"      Initial: {state.getPotentialEnergy()}")
    simulation.minimizeEnergy(maxIterations=MINIMIZATION_STEPS)
    state = simulation.context.getState(getEnergy=True)
    print(f"      Final: {state.getPotentialEnergy()}")

    # NVT
    print(f"\n[6/8] NVT EQUILIBRATION ({EQUILIBRATION_NVT_PS} ps)")
    simulation.context.setVelocitiesToTemperature(TEMPERATURE, random_seed)
    nvt_steps = int(EQUILIBRATION_NVT_PS * 1000 / TIMESTEP.value_in_unit(femtoseconds))
    simulation.step(nvt_steps)

    # NPT
    print(f"\n[7/8] NPT EQUILIBRATION ({EQUILIBRATION_NPT_PS} ps)")
    system.addForce(MonteCarloBarostat(PRESSURE, TEMPERATURE, 25))
    simulation.context.reinitialize(preserveState=True)
    npt_steps = int(EQUILIBRATION_NPT_PS * 1000 / TIMESTEP.value_in_unit(femtoseconds))
    simulation.step(npt_steps)

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

    print(f"\n>>> COMPLETE: {output_prefix}")
    print(f"    Finished: {datetime.now().strftime('%H:%M:%S')}")

    # Cleanup
    del simulation, system, integrator, modeller
    gc.collect()

    return True

# ============================================================
# MAIN
# ============================================================
def main():
    print("=" * 70)
    print("VAST.AI RTX 5090 - NOD2-SCOUT MD SIMULATIONS")
    print(f"Running {len(SIMULATIONS)} simulations ({len(SIMULATIONS) * PRODUCTION_NS}ns total)")
    print(f"Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("=" * 70)

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
    success_count = sum(1 for _, s in results if s)
    for name, success in results:
        status = "COMPLETE" if success else "FAILED"
        print(f"  {name}: {status}")
    print(f"\nTotal: {success_count}/{len(results)} successful")
    print(f"Finished: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("=" * 70)

    # AUTO-SHUTDOWN
    print("\n" + "!" * 70)
    print("ALL SIMULATIONS COMPLETE - SHUTTING DOWN IN 60 SECONDS")
    print("!" * 70)

    import time
    time.sleep(60)

    print("SHUTTING DOWN NOW...")
    os.system('sudo shutdown now')

if __name__ == "__main__":
    # Change to PHASE_6 directory
    script_dir = os.path.dirname(os.path.abspath(__file__))
    phase6_dir = os.path.join(script_dir, "../..")
    os.chdir(phase6_dir)
    print(f"Working directory: {os.getcwd()}")

    main()
