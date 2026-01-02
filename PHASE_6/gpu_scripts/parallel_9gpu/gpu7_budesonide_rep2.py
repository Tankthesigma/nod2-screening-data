#!/usr/bin/env python3
"""
GPU 7: BUDESONIDE Replicate 2 (POSITIVE CONTROL)
RTX 5090 + HMR (4fs timestep) - ~30 min per 20ns

FIXED: OpenFF parameterization, proper NVT→NPT, random seeds
"""

from openmm import *
from openmm.app import *
from openmm.unit import *
import sys
import os
from datetime import datetime
from openmmforcefields.generators import SystemGenerator
from openff.toolkit import Molecule

GPU_ID = 7
SYSTEM_NAME = "budesonide"
REPLICATE = 2
PDB_FILE = "../structures/complex_budesonide.pdb"
SDF_FILE = "../structures/budesonide_docked.sdf"
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

REPORT_INTERVAL = 2500
TRAJECTORY_INTERVAL = 1250
CHECKPOINT_INTERVAL = 25000

def run():
    print("=" * 60)
    print(f"GPU {GPU_ID}: {SYSTEM_NAME.upper()} - Replicate {REPLICATE}")
    print(f"[POSITIVE CONTROL]")
    print(f"Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("=" * 60)

    output_prefix = f"{SYSTEM_NAME}_rep{REPLICATE}"
    ligand_mol = Molecule.from_file(SDF_FILE)
    pdb = PDBFile(PDB_FILE)

    forcefield_kwargs = {'constraints': AllBonds, 'rigidWater': True,
        'removeCMMotion': True, 'hydrogenMass': HYDROGEN_MASS if USE_HMR else None}
    system_generator = SystemGenerator(
        forcefields=['amber14-all.xml', 'amber14/tip3pfb.xml'],
        small_molecule_forcefield='openff-2.1.0',
        molecules=[ligand_mol], forcefield_kwargs=forcefield_kwargs)

    modeller = Modeller(pdb.topology, pdb.positions)
    modeller.addSolvent(system_generator.forcefield, model='tip3p',
                        padding=1.0*nanometers, ionicStrength=0.15*molar)
    print(f"Solvated atoms: {modeller.topology.getNumAtoms()}")

    system = system_generator.create_system(modeller.topology,
        nonbondedMethod=PME, nonbondedCutoff=1.0*nanometers)

    integrator = LangevinMiddleIntegrator(TEMPERATURE, FRICTION, TIMESTEP)
    integrator.setRandomNumberSeed(RANDOM_SEED)
    platform = Platform.getPlatformByName('CUDA')
    properties = {'Precision': 'mixed', 'DeviceIndex': str(GPU_ID)}
    simulation = Simulation(modeller.topology, system, integrator, platform, properties)
    simulation.context.setPositions(modeller.positions)

    print(">>> MINIMIZATION")
    simulation.minimizeEnergy(maxIterations=MINIMIZATION_STEPS)

    print(f">>> NVT EQUILIBRATION ({EQUILIBRATION_NVT_PS} ps)")
    simulation.context.setVelocitiesToTemperature(TEMPERATURE, RANDOM_SEED)
    simulation.step(int(EQUILIBRATION_NVT_PS * 1000 / TIMESTEP.value_in_unit(femtoseconds)))

    system.addForce(MonteCarloBarostat(PRESSURE, TEMPERATURE, 25))
    simulation.context.reinitialize(preserveState=True)

    print(f">>> NPT EQUILIBRATION ({EQUILIBRATION_NPT_PS} ps)")
    simulation.step(int(EQUILIBRATION_NPT_PS * 1000 / TIMESTEP.value_in_unit(femtoseconds)))

    print(f">>> PRODUCTION ({PRODUCTION_NS} ns)")
    production_steps = int(PRODUCTION_NS * 1e6 / TIMESTEP.value_in_unit(femtoseconds))
    simulation.reporters.append(DCDReporter(f'{output_prefix}.dcd', TRAJECTORY_INTERVAL))
    simulation.reporters.append(StateDataReporter(f'{output_prefix}.log', REPORT_INTERVAL,
        step=True, time=True, potentialEnergy=True, kineticEnergy=True,
        temperature=True, speed=True, progress=True, remainingTime=True, totalSteps=production_steps))
    simulation.reporters.append(CheckpointReporter(f'{output_prefix}.chk', CHECKPOINT_INTERVAL))
    simulation.step(production_steps)
    simulation.saveState(f'{output_prefix}_final.xml')

    print(f"COMPLETE: {output_prefix} at {datetime.now().strftime('%H:%M:%S')}")

if __name__ == "__main__":
    os.chdir(os.path.dirname(os.path.abspath(__file__)) + "/..")
    run()
