#!/usr/bin/env python3
"""
GPU 2: FEBUXOSTAT Replicate 3
RTX 5090 + HMR (4fs timestep) - ~30 min per 20ns
"""

from openmm import *
from openmm.app import *
from openmm.unit import *
import sys
import os

# ============================================================
# CONFIGURATION
# ============================================================
GPU_ID = 2
SYSTEM_NAME = "febuxostat"
REPLICATE = 3
PDB_FILE = "../structures/complex_febuxostat.pdb"

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

# Output
REPORT_INTERVAL = 5000      # Every 10ps
TRAJECTORY_INTERVAL = 2500  # Every 5ps
CHECKPOINT_INTERVAL = 25000 # Every 50ps

# ============================================================
# MAIN SIMULATION
# ============================================================
def run():
    print("=" * 60)
    print(f"GPU {GPU_ID}: {SYSTEM_NAME.upper()} - Replicate {REPLICATE}")
    print("=" * 60)

    output_prefix = f"{SYSTEM_NAME}_rep{REPLICATE}"

    print("\n[1/6] Loading structure...")
    pdb = PDBFile(PDB_FILE)
    print(f"      Atoms: {pdb.topology.getNumAtoms()}")

    print("[2/6] Setting up force field...")
    forcefield = ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')

    print("[3/6] Adding solvent...")
    modeller = Modeller(pdb.topology, pdb.positions)
    modeller.addSolvent(forcefield, model='tip3p',
                        padding=1.0*nanometers,
                        ionicStrength=0.15*molar)
    print(f"      Solvated atoms: {modeller.topology.getNumAtoms()}")

    print("[4/6] Creating system (HMR enabled)...")
    system = forcefield.createSystem(
        modeller.topology,
        nonbondedMethod=PME,
        nonbondedCutoff=1.0*nanometers,
        constraints=AllBonds,
        hydrogenMass=HYDROGEN_MASS if USE_HMR else None
    )
    system.addForce(MonteCarloBarostat(PRESSURE, TEMPERATURE))

    integrator = LangevinMiddleIntegrator(TEMPERATURE, FRICTION, TIMESTEP)

    print(f"[5/6] Setting up CUDA (GPU {GPU_ID})...")
    platform = Platform.getPlatformByName('CUDA')
    properties = {'Precision': 'mixed', 'DeviceIndex': str(GPU_ID)}

    simulation = Simulation(modeller.topology, system, integrator, platform, properties)
    simulation.context.setPositions(modeller.positions)

    print("\n>>> MINIMIZATION")
    state = simulation.context.getState(getEnergy=True)
    print(f"    Initial: {state.getPotentialEnergy()}")
    simulation.minimizeEnergy(maxIterations=MINIMIZATION_STEPS)
    state = simulation.context.getState(getEnergy=True)
    print(f"    Final: {state.getPotentialEnergy()}")

    print(f"\n>>> NVT EQUILIBRATION ({EQUILIBRATION_NVT_PS} ps)")
    simulation.context.setVelocitiesToTemperature(TEMPERATURE)
    nvt_steps = int(EQUILIBRATION_NVT_PS * 1000 / TIMESTEP.value_in_unit(femtoseconds))
    simulation.step(nvt_steps)

    print(f"\n>>> NPT EQUILIBRATION ({EQUILIBRATION_NPT_PS} ps)")
    npt_steps = int(EQUILIBRATION_NPT_PS * 1000 / TIMESTEP.value_in_unit(femtoseconds))
    simulation.step(npt_steps)

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
    print("=" * 60)

if __name__ == "__main__":
    os.chdir(os.path.dirname(os.path.abspath(__file__)) + "/..")
    run()
