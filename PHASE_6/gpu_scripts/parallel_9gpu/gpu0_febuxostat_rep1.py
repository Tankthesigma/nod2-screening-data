#!/usr/bin/env python3
"""
GPU 0: FEBUXOSTAT Replicate 1
RTX 5090 + HMR (4fs timestep) - ~30 min per 20ns

FIXED VERSION:
- OpenFF SystemGenerator for ligand parameterization
- Proper NVT→NPT equilibration (barostat added after NVT)
- Random seeds for reproducibility
- Corrected interval comments
"""

from openmm import *
from openmm.app import *
from openmm.unit import *
import sys
import os
from datetime import datetime

# OpenFF for ligand parameterization
from openmmforcefields.generators import SystemGenerator
from openff.toolkit import Molecule

# ============================================================
# CONFIGURATION
# ============================================================
GPU_ID = 0
SYSTEM_NAME = "febuxostat"
REPLICATE = 1
PDB_FILE = "../structures/complex_febuxostat.pdb"
SDF_FILE = "../structures/febuxostat_docked.sdf"  # For ligand parameters

# Random seed for reproducibility
RANDOM_SEED = 1000 * GPU_ID + REPLICATE  # Unique per GPU/replicate

# Physics
TEMPERATURE = 310.15 * kelvin
PRESSURE = 1.0 * atmospheres
TIMESTEP = 4.0 * femtoseconds  # 4fs with HMR
FRICTION = 1.0 / picoseconds

# HMR - Hydrogen Mass Repartitioning
USE_HMR = True
HYDROGEN_MASS = 1.5 * amu

# Protocol
MINIMIZATION_STEPS = 5000
EQUILIBRATION_NVT_PS = 100
EQUILIBRATION_NPT_PS = 500
PRODUCTION_NS = 20

# Output intervals (at 4fs timestep):
REPORT_INTERVAL = 2500      # Every 10ps (2500 * 4fs = 10ps)
TRAJECTORY_INTERVAL = 1250  # Every 5ps (1250 * 4fs = 5ps)
CHECKPOINT_INTERVAL = 25000 # Every 100ps (25000 * 4fs = 100ps)

# ============================================================
# MAIN SIMULATION
# ============================================================
def run():
    print("=" * 60)
    print(f"GPU {GPU_ID}: {SYSTEM_NAME.upper()} - Replicate {REPLICATE}")
    print(f"Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"Random seed: {RANDOM_SEED}")
    print("=" * 60)

    output_prefix = f"{SYSTEM_NAME}_rep{REPLICATE}"

    # --------------------------------------------------------
    # 1. Load ligand molecule for parameterization
    # --------------------------------------------------------
    print("\n[1/7] Loading ligand from SDF for parameterization...")
    ligand_mol = Molecule.from_file(SDF_FILE)
    print(f"      Ligand: {ligand_mol.name or SYSTEM_NAME}")
    print(f"      Atoms: {ligand_mol.n_atoms}")

    # --------------------------------------------------------
    # 2. Load complex structure
    # --------------------------------------------------------
    print("\n[2/7] Loading complex structure...")
    pdb = PDBFile(PDB_FILE)
    print(f"      Complex atoms: {pdb.topology.getNumAtoms()}")

    # --------------------------------------------------------
    # 3. Create SystemGenerator with OpenFF for ligand
    # --------------------------------------------------------
    print("\n[3/7] Setting up force fields (Amber14 + OpenFF)...")

    # Force field kwargs
    forcefield_kwargs = {
        'constraints': AllBonds,
        'rigidWater': True,
        'removeCMMotion': True,
        'hydrogenMass': HYDROGEN_MASS if USE_HMR else None
    }

    system_generator = SystemGenerator(
        forcefields=['amber14-all.xml', 'amber14/tip3pfb.xml'],
        small_molecule_forcefield='openff-2.1.0',  # OpenFF Sage
        molecules=[ligand_mol],
        forcefield_kwargs=forcefield_kwargs
    )
    print("      Protein: Amber14")
    print("      Ligand: OpenFF 2.1.0 (Sage)")
    print("      Water: TIP3P-FB")

    # --------------------------------------------------------
    # 4. Solvate the system
    # --------------------------------------------------------
    print("\n[4/7] Adding solvent and ions...")
    modeller = Modeller(pdb.topology, pdb.positions)

    # Need to use the system_generator's forcefield for solvation
    modeller.addSolvent(
        system_generator.forcefield,
        model='tip3p',
        padding=1.0*nanometers,
        ionicStrength=0.15*molar
    )
    print(f"      Solvated atoms: {modeller.topology.getNumAtoms()}")

    # --------------------------------------------------------
    # 5. Create system (WITHOUT barostat - added after NVT)
    # --------------------------------------------------------
    print("\n[5/7] Creating system (HMR enabled)...")
    system = system_generator.create_system(
        modeller.topology,
        nonbondedMethod=PME,
        nonbondedCutoff=1.0*nanometers
    )
    print(f"      Forces: {system.getNumForces()}")

    # --------------------------------------------------------
    # 6. Setup simulation
    # --------------------------------------------------------
    print(f"\n[6/7] Setting up CUDA (GPU {GPU_ID})...")
    integrator = LangevinMiddleIntegrator(TEMPERATURE, FRICTION, TIMESTEP)
    integrator.setRandomNumberSeed(RANDOM_SEED)

    platform = Platform.getPlatformByName('CUDA')
    properties = {'Precision': 'mixed', 'DeviceIndex': str(GPU_ID)}

    simulation = Simulation(modeller.topology, system, integrator, platform, properties)
    simulation.context.setPositions(modeller.positions)

    # --------------------------------------------------------
    # MINIMIZATION
    # --------------------------------------------------------
    print("\n>>> MINIMIZATION")
    state = simulation.context.getState(getEnergy=True)
    print(f"    Initial energy: {state.getPotentialEnergy()}")
    simulation.minimizeEnergy(maxIterations=MINIMIZATION_STEPS)
    state = simulation.context.getState(getEnergy=True)
    print(f"    Final energy: {state.getPotentialEnergy()}")

    # --------------------------------------------------------
    # NVT EQUILIBRATION (no barostat = true NVT)
    # --------------------------------------------------------
    print(f"\n>>> NVT EQUILIBRATION ({EQUILIBRATION_NVT_PS} ps)")
    simulation.context.setVelocitiesToTemperature(TEMPERATURE, RANDOM_SEED)
    nvt_steps = int(EQUILIBRATION_NVT_PS * 1000 / TIMESTEP.value_in_unit(femtoseconds))
    simulation.step(nvt_steps)
    state = simulation.context.getState(getEnergy=True)
    print(f"    Final energy: {state.getPotentialEnergy()}")

    # --------------------------------------------------------
    # ADD BAROSTAT FOR NPT (after NVT)
    # --------------------------------------------------------
    print("\n>>> Adding barostat for NPT...")
    system.addForce(MonteCarloBarostat(PRESSURE, TEMPERATURE, 25))
    simulation.context.reinitialize(preserveState=True)

    # --------------------------------------------------------
    # NPT EQUILIBRATION
    # --------------------------------------------------------
    print(f"\n>>> NPT EQUILIBRATION ({EQUILIBRATION_NPT_PS} ps)")
    npt_steps = int(EQUILIBRATION_NPT_PS * 1000 / TIMESTEP.value_in_unit(femtoseconds))
    simulation.step(npt_steps)
    state = simulation.context.getState(getEnergy=True)
    print(f"    Final energy: {state.getPotentialEnergy()}")

    # --------------------------------------------------------
    # PRODUCTION
    # --------------------------------------------------------
    print(f"\n>>> PRODUCTION ({PRODUCTION_NS} ns)")
    production_steps = int(PRODUCTION_NS * 1e6 / TIMESTEP.value_in_unit(femtoseconds))

    simulation.reporters.append(DCDReporter(f'{output_prefix}.dcd', TRAJECTORY_INTERVAL))
    simulation.reporters.append(StateDataReporter(
        f'{output_prefix}.log', REPORT_INTERVAL,
        step=True, time=True, potentialEnergy=True, kineticEnergy=True,
        temperature=True, speed=True, progress=True, remainingTime=True,
        totalSteps=production_steps
    ))
    simulation.reporters.append(CheckpointReporter(f'{output_prefix}.chk', CHECKPOINT_INTERVAL))

    simulation.step(production_steps)
    simulation.saveState(f'{output_prefix}_final.xml')

    print("\n" + "=" * 60)
    print(f"COMPLETE: {output_prefix}")
    print(f"Finished: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("=" * 60)

if __name__ == "__main__":
    os.chdir(os.path.dirname(os.path.abspath(__file__)) + "/..")
    run()
