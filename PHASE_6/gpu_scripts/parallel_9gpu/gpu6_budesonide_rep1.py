#!/usr/bin/env python3
"""
GPU 6: BUDESONIDE Replicate 1
RTX 5090 + HMR (4fs timestep) - ~30 min per 20ns

FULLY FIXED VERSION:
- OpenFF SystemGenerator with GAFF fallback
- Proper output directories
- Ligand-PDB sanity checks
- Optimized I/O intervals
"""

from openmm import *
from openmm.app import *
from openmm.unit import *
import sys
import os
from datetime import datetime
import traceback
import gc

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
GPU_ID = 6
SYSTEM_NAME = "budesonide"
REPLICATE = 1
PDB_FILE = "../structures/complex_budesonide.pdb"
SDF_FILE = "../structures/budesonide_docked.sdf"

OUTPUT_DIR = "../trajectories"
LOG_DIR = "../logs"
CHECKPOINT_DIR = "../checkpoints"

RANDOM_SEED = 1000 * GPU_ID + REPLICATE

TEMPERATURE = 310.15 * kelvin
PRESSURE = 1.0 * atmospheres
TIMESTEP = 4.0 * femtoseconds
FRICTION = 1.0 / picoseconds
USE_HMR = True
HYDROGEN_MASS = 1.5 * amu

MINIMIZATION_STEPS = 5000
EQUILIBRATION_NVT_PS = 100
EQUILIBRATION_NPT_PS = 500
PRODUCTION_NS = 20

REPORT_INTERVAL = 2500       # 10ps
TRAJECTORY_INTERVAL = 2500   # 10ps
CHECKPOINT_INTERVAL = 250000 # 1ns

def validate_inputs():
    errors = []
    if not os.path.exists(PDB_FILE):
        errors.append(f"PDB not found: {PDB_FILE}")
    if not os.path.exists(SDF_FILE):
        errors.append(f"SDF not found: {SDF_FILE}")
    if errors:
        for e in errors:
            print(f"ERROR: {e}")
        sys.exit(1)

def setup_directories():
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    os.makedirs(LOG_DIR, exist_ok=True)
    os.makedirs(CHECKPOINT_DIR, exist_ok=True)

def create_system_with_fallback(pdb, ligand_mol, forcefield_kwargs):
    modeller = Modeller(pdb.topology, pdb.positions)

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
            system = system_generator.create_system(modeller.topology,
                nonbondedMethod=PME, nonbondedCutoff=1.0*nanometers)
            print("      SUCCESS: OpenFF")
            return system, modeller
        except Exception as e:
            print(f"      OpenFF failed: {e}, trying GAFF...")

    try:
        forcefield = ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')
        gaff_generator = GAFFTemplateGenerator(molecules=[ligand_mol])
        forcefield.registerTemplateGenerator(gaff_generator.generator)
        modeller = Modeller(pdb.topology, pdb.positions)
        modeller.addSolvent(forcefield, model='tip3p',
                          padding=1.0*nanometers, ionicStrength=0.15*molar)
        system = forcefield.createSystem(modeller.topology,
            nonbondedMethod=PME, nonbondedCutoff=1.0*nanometers,
            constraints=AllBonds, hydrogenMass=HYDROGEN_MASS if USE_HMR else None)
        print("      SUCCESS: GAFF")
        return system, modeller
    except Exception as e:
        raise RuntimeError(f"Parameterization failed: {e}")

def run():
    print("=" * 70)
    print(f"GPU {GPU_ID}: {SYSTEM_NAME.upper()} - Replicate {REPLICATE}")
    print(f"Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("=" * 70)

    validate_inputs()
    setup_directories()
    output_prefix = f"{SYSTEM_NAME}_rep{REPLICATE}"

    print("\n[1/8] Loading ligand...")
    ligand_mol = Molecule.from_file(SDF_FILE)
    print(f"      Atoms: {ligand_mol.n_atoms}")

    print("\n[2/8] Loading structure...")
    pdb = PDBFile(PDB_FILE)
    print(f"      Atoms: {pdb.topology.getNumAtoms()}")

    print("\n[3/8] Setting up force fields...")
    forcefield_kwargs = {'constraints': AllBonds, 'rigidWater': True,
        'removeCMMotion': True, 'hydrogenMass': HYDROGEN_MASS if USE_HMR else None}
    system, modeller = create_system_with_fallback(pdb, ligand_mol, forcefield_kwargs)
    print(f"      Solvated: {modeller.topology.getNumAtoms()}")

    print(f"\n[4/8] Setting up CUDA (GPU {GPU_ID})...")
    integrator = LangevinMiddleIntegrator(TEMPERATURE, FRICTION, TIMESTEP)
    integrator.setRandomNumberSeed(RANDOM_SEED)
    platform = Platform.getPlatformByName('CUDA')
    properties = {'Precision': 'mixed', 'DeviceIndex': str(GPU_ID)}
    simulation = Simulation(modeller.topology, system, integrator, platform, properties)
    simulation.context.setPositions(modeller.positions)

    print("\n[5/8] MINIMIZATION")
    simulation.minimizeEnergy(maxIterations=MINIMIZATION_STEPS)

    print(f"\n[6/8] NVT EQUILIBRATION ({EQUILIBRATION_NVT_PS} ps)")
    simulation.context.setVelocitiesToTemperature(TEMPERATURE, RANDOM_SEED)
    simulation.step(int(EQUILIBRATION_NVT_PS * 1000 / TIMESTEP.value_in_unit(femtoseconds)))

    print(f"\n[7/8] NPT EQUILIBRATION ({EQUILIBRATION_NPT_PS} ps)")
    system.addForce(MonteCarloBarostat(PRESSURE, TEMPERATURE, 25))
    simulation.context.reinitialize(preserveState=True)
    simulation.step(int(EQUILIBRATION_NPT_PS * 1000 / TIMESTEP.value_in_unit(femtoseconds)))

    print(f"\n[8/8] PRODUCTION ({PRODUCTION_NS} ns)")
    production_steps = int(PRODUCTION_NS * 1e6 / TIMESTEP.value_in_unit(femtoseconds))
    simulation.reporters.append(DCDReporter(f'{OUTPUT_DIR}/{output_prefix}.dcd', TRAJECTORY_INTERVAL))
    simulation.reporters.append(StateDataReporter(f'{LOG_DIR}/{output_prefix}.log', REPORT_INTERVAL,
        step=True, time=True, potentialEnergy=True, kineticEnergy=True,
        temperature=True, speed=True, progress=True, remainingTime=True, totalSteps=production_steps))
    simulation.reporters.append(CheckpointReporter(f'{CHECKPOINT_DIR}/{output_prefix}.chk', CHECKPOINT_INTERVAL))
    simulation.step(production_steps)
    simulation.saveState(f'{CHECKPOINT_DIR}/{output_prefix}_final.xml')

    print(f"\nCOMPLETE: {output_prefix} at {datetime.now().strftime('%H:%M:%S')}")

if __name__ == "__main__":
    os.chdir(os.path.dirname(os.path.abspath(__file__)) + "/..")
    try:
        run()
    except Exception as e:
        print(f"FATAL ERROR: {e}")
        traceback.print_exc()
        sys.exit(1)
