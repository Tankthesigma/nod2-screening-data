#!/usr/bin/env python3
"""
GPU 8: REMAINING SIMULATIONS (Sequential with Error Isolation)
- budesonide_rep3
- natural_cid10592_rep1, rep2, rep3
- apo_rep1 (NEGATIVE CONTROL - no ligand)
- decoy_rep1 (NEGATIVE CONTROL - weak binder)

FULLY FIXED VERSION:
- Try/except around each simulation (error isolation)
- Proper output directories
- GAFF fallback for parameterization
- Memory cleanup between simulations
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
    print("WARNING: OpenFF not available")

# ============================================================
# CONFIGURATION
# ============================================================
GPU_ID = 8

# Output directories
OUTPUT_DIR = "../trajectories"
LOG_DIR = "../logs"
CHECKPOINT_DIR = "../checkpoints"

# Physics
TEMPERATURE = 310.15 * kelvin
PRESSURE = 1.0 * atmospheres
TIMESTEP = 4.0 * femtoseconds
FRICTION = 1.0 / picoseconds
USE_HMR = True
HYDROGEN_MASS = 1.5 * amu

# Protocol
MINIMIZATION_STEPS = 5000
EQUILIBRATION_NVT_PS = 100
EQUILIBRATION_NPT_PS = 500
PRODUCTION_NS = 20

# Output intervals (optimized for I/O)
REPORT_INTERVAL = 2500       # 10ps
TRAJECTORY_INTERVAL = 2500   # 10ps
CHECKPOINT_INTERVAL = 250000 # 1ns

# Simulations to run sequentially
# Format: (system_name, replicate, pdb_file, sdf_file_or_None)
SIMULATIONS = [
    ("budesonide", 3, "../structures/complex_budesonide.pdb", "../structures/budesonide_docked.sdf"),
    ("natural_cid10592", 1, "../structures/complex_natural_cid10592.pdb", "../structures/natural_top_docked.sdf"),
    ("natural_cid10592", 2, "../structures/complex_natural_cid10592.pdb", "../structures/natural_top_docked.sdf"),
    ("natural_cid10592", 3, "../structures/complex_natural_cid10592.pdb", "../structures/natural_top_docked.sdf"),
    ("apo", 1, "../structures/complex_apo.pdb", None),
    ("decoy", 1, "../structures/complex_decoy.pdb", "../structures/decoy_docked.sdf"),
]

# ============================================================
# HELPER FUNCTIONS
# ============================================================
def setup_directories():
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    os.makedirs(LOG_DIR, exist_ok=True)
    os.makedirs(CHECKPOINT_DIR, exist_ok=True)

def create_system_with_fallback(pdb, ligand_mol, forcefield_kwargs, is_apo=False):
    """Create system with OpenFF, falling back to GAFF if needed."""
    modeller = Modeller(pdb.topology, pdb.positions)

    if is_apo:
        # APO system - no ligand, use standard ForceField
        forcefield = ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')
        modeller.addSolvent(forcefield, model='tip3p',
                          padding=1.0*nanometers, ionicStrength=0.15*molar)
        system = forcefield.createSystem(
            modeller.topology,
            nonbondedMethod=PME,
            nonbondedCutoff=1.0*nanometers,
            constraints=AllBonds,
            hydrogenMass=HYDROGEN_MASS if USE_HMR else None
        )
        return system, modeller

    # Try OpenFF first
    if OPENFF_AVAILABLE:
        try:
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
            return system, modeller
        except Exception as e:
            print(f"      OpenFF failed: {e}, trying GAFF...")

    # Fallback to GAFF
    forcefield = ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')
    gaff_generator = GAFFTemplateGenerator(molecules=[ligand_mol])
    forcefield.registerTemplateGenerator(gaff_generator.generator)
    modeller = Modeller(pdb.topology, pdb.positions)
    modeller.addSolvent(forcefield, model='tip3p',
                      padding=1.0*nanometers, ionicStrength=0.15*molar)
    system = forcefield.createSystem(
        modeller.topology, nonbondedMethod=PME, nonbondedCutoff=1.0*nanometers,
        constraints=AllBonds, hydrogenMass=HYDROGEN_MASS if USE_HMR else None
    )
    return system, modeller

# ============================================================
# SIMULATION FUNCTION (with error isolation)
# ============================================================
def run_simulation(system_name, replicate, pdb_file, sdf_file, sim_num, total_sims):
    """
    Run a single MD simulation with full error isolation.
    If this simulation fails, the loop continues to the next one.
    """
    random_seed = 8000 + sim_num * 100 + replicate
    is_apo = sdf_file is None
    output_prefix = f"{system_name}_rep{replicate}"

    print("\n" + "=" * 70)
    print(f"GPU {GPU_ID} | SIMULATION {sim_num}/{total_sims}")
    print(f"{system_name.upper()} - Replicate {replicate}")
    print(f"Started: {datetime.now().strftime('%H:%M:%S')}")
    print(f"Random seed: {random_seed}")
    print("=" * 70)

    # Check input files
    if not os.path.exists(pdb_file):
        print(f"ERROR: PDB file not found: {pdb_file}")
        return False

    if not is_apo and not os.path.exists(sdf_file):
        print(f"ERROR: SDF file not found: {sdf_file}")
        return False

    # Load ligand
    if not is_apo:
        print("\n[1/8] Loading ligand...")
        ligand_mol = Molecule.from_file(sdf_file)
        print(f"      Atoms: {ligand_mol.n_atoms}")
    else:
        print("\n[1/8] APO system - no ligand")
        ligand_mol = None

    # Load structure
    print("\n[2/8] Loading structure...")
    pdb = PDBFile(pdb_file)
    print(f"      Atoms: {pdb.topology.getNumAtoms()}")

    # Create system
    print("\n[3/8] Setting up force fields...")
    forcefield_kwargs = {
        'constraints': AllBonds,
        'rigidWater': True,
        'removeCMMotion': True,
        'hydrogenMass': HYDROGEN_MASS if USE_HMR else None
    }
    system, modeller = create_system_with_fallback(pdb, ligand_mol, forcefield_kwargs, is_apo)
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

    simulation.reporters.append(DCDReporter(f'{OUTPUT_DIR}/{output_prefix}.dcd', TRAJECTORY_INTERVAL))
    simulation.reporters.append(StateDataReporter(
        f'{LOG_DIR}/{output_prefix}.log', REPORT_INTERVAL,
        step=True, time=True, potentialEnergy=True, kineticEnergy=True,
        temperature=True, speed=True, progress=True, remainingTime=True,
        totalSteps=production_steps
    ))
    simulation.reporters.append(CheckpointReporter(f'{CHECKPOINT_DIR}/{output_prefix}.chk', CHECKPOINT_INTERVAL))

    simulation.step(production_steps)
    simulation.saveState(f'{CHECKPOINT_DIR}/{output_prefix}_final.xml')

    print(f"\n>>> COMPLETE: {output_prefix}")
    print(f"    Finished: {datetime.now().strftime('%H:%M:%S')}")

    # Cleanup
    del simulation
    del system
    del integrator
    del modeller
    gc.collect()

    return True

# ============================================================
# MAIN
# ============================================================
def main():
    print("=" * 70)
    print("GPU 8: SEQUENTIAL SIMULATIONS (with error isolation)")
    print(f"Running {len(SIMULATIONS)} simulations")
    print(f"Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("=" * 70)

    setup_directories()

    results = []
    for i, (name, rep, pdb, sdf) in enumerate(SIMULATIONS, 1):
        try:
            success = run_simulation(name, rep, pdb, sdf, i, len(SIMULATIONS))
            results.append((f"{name}_rep{rep}", success))
        except Exception as e:
            print(f"\n{'='*70}")
            print(f"ERROR in {name}_rep{rep}: {e}")
            print(f"{'='*70}")
            traceback.print_exc()
            results.append((f"{name}_rep{rep}", False))
            # Continue to next simulation instead of crashing
            print("\nContinuing to next simulation...")
            gc.collect()  # Clean up before next sim

    # Summary
    print("\n" + "=" * 70)
    print("GPU 8 SUMMARY")
    print("=" * 70)
    success_count = sum(1 for _, s in results if s)
    for name, success in results:
        status = "COMPLETE" if success else "FAILED"
        print(f"  {name}: {status}")
    print(f"\nTotal: {success_count}/{len(results)} successful")
    print(f"Finished: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("=" * 70)

if __name__ == "__main__":
    os.chdir(os.path.dirname(os.path.abspath(__file__)) + "/..")
    main()
