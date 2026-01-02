#!/usr/bin/env python3
"""
GPU 8: REMAINING SIMULATIONS (Sequential)
- budesonide_rep3
- natural_cid10592_rep1, rep2, rep3
- apo_rep1
- decoy_rep1

RTX 5090 + HMR (4fs timestep) - ~30 min each = ~3 hours total
"""

from openmm import *
from openmm.app import *
from openmm.unit import *
import sys
import os
from datetime import datetime

# ============================================================
# CONFIGURATION
# ============================================================
GPU_ID = 8

# Physics
TEMPERATURE = 310.15 * kelvin
PRESSURE = 1.0 * atmospheres
TIMESTEP = 4.0 * femtoseconds  # 4fs with HMR
FRICTION = 1.0 / picoseconds

# HMR
USE_HMR = True
HYDROGEN_MASS = 1.5 * amu

# Protocol
MINIMIZATION_STEPS = 5000
EQUILIBRATION_NVT_PS = 100
EQUILIBRATION_NPT_PS = 500
PRODUCTION_NS = 20

# Output
REPORT_INTERVAL = 5000
TRAJECTORY_INTERVAL = 2500
CHECKPOINT_INTERVAL = 25000

# Simulations to run sequentially
SIMULATIONS = [
    ("budesonide", 3, "../structures/complex_budesonide.pdb"),
    ("natural_cid10592", 1, "../structures/complex_natural_cid10592.pdb"),
    ("natural_cid10592", 2, "../structures/complex_natural_cid10592.pdb"),
    ("natural_cid10592", 3, "../structures/complex_natural_cid10592.pdb"),
    ("apo", 1, "../structures/complex_apo.pdb"),
    ("decoy", 1, "../structures/complex_decoy.pdb"),
]

# ============================================================
# SIMULATION FUNCTION
# ============================================================
def run_simulation(system_name, replicate, pdb_file, sim_num, total_sims):
    print("\n" + "=" * 70)
    print(f"GPU {GPU_ID} | SIMULATION {sim_num}/{total_sims}")
    print(f"{system_name.upper()} - Replicate {replicate}")
    print(f"Started: {datetime.now().strftime('%H:%M:%S')}")
    print("=" * 70)

    output_prefix = f"{system_name}_rep{replicate}"

    # Load structure
    print("\n[1/6] Loading structure...")
    pdb = PDBFile(pdb_file)
    print(f"      Atoms: {pdb.topology.getNumAtoms()}")

    # Force field
    print("[2/6] Setting up force field...")
    forcefield = ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')

    # Add solvent
    print("[3/6] Adding solvent...")
    modeller = Modeller(pdb.topology, pdb.positions)
    modeller.addSolvent(forcefield, model='tip3p',
                        padding=1.0*nanometers,
                        ionicStrength=0.15*molar)
    print(f"      Solvated atoms: {modeller.topology.getNumAtoms()}")

    # Create system with HMR
    print("[4/6] Creating system (HMR enabled)...")
    system = forcefield.createSystem(
        modeller.topology,
        nonbondedMethod=PME,
        nonbondedCutoff=1.0*nanometers,
        constraints=AllBonds,
        hydrogenMass=HYDROGEN_MASS if USE_HMR else None
    )
    system.addForce(MonteCarloBarostat(PRESSURE, TEMPERATURE))

    # Integrator
    integrator = LangevinMiddleIntegrator(TEMPERATURE, FRICTION, TIMESTEP)

    # Platform
    print(f"[5/6] Setting up CUDA (GPU {GPU_ID})...")
    platform = Platform.getPlatformByName('CUDA')
    properties = {'Precision': 'mixed', 'DeviceIndex': str(GPU_ID)}

    simulation = Simulation(modeller.topology, system, integrator, platform, properties)
    simulation.context.setPositions(modeller.positions)

    # Minimization
    print("\n>>> MINIMIZATION")
    state = simulation.context.getState(getEnergy=True)
    print(f"    Initial: {state.getPotentialEnergy()}")
    simulation.minimizeEnergy(maxIterations=MINIMIZATION_STEPS)
    state = simulation.context.getState(getEnergy=True)
    print(f"    Final: {state.getPotentialEnergy()}")

    # NVT
    print(f"\n>>> NVT EQUILIBRATION ({EQUILIBRATION_NVT_PS} ps)")
    simulation.context.setVelocitiesToTemperature(TEMPERATURE)
    nvt_steps = int(EQUILIBRATION_NVT_PS * 1000 / TIMESTEP.value_in_unit(femtoseconds))
    simulation.step(nvt_steps)

    # NPT
    print(f"\n>>> NPT EQUILIBRATION ({EQUILIBRATION_NPT_PS} ps)")
    npt_steps = int(EQUILIBRATION_NPT_PS * 1000 / TIMESTEP.value_in_unit(femtoseconds))
    simulation.step(npt_steps)

    # Production
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

    print(f"\n>>> COMPLETE: {output_prefix}")
    print(f"    Finished: {datetime.now().strftime('%H:%M:%S')}")

    # Clear memory
    del simulation, system, integrator, modeller
    import gc
    gc.collect()

# ============================================================
# MAIN
# ============================================================
def main():
    print("=" * 70)
    print("GPU 8: SEQUENTIAL SIMULATIONS")
    print(f"Running {len(SIMULATIONS)} simulations")
    print(f"Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("=" * 70)

    for i, (name, rep, pdb) in enumerate(SIMULATIONS, 1):
        run_simulation(name, rep, pdb, i, len(SIMULATIONS))

    print("\n" + "=" * 70)
    print("ALL GPU 8 SIMULATIONS COMPLETE!")
    print(f"Finished: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("=" * 70)

if __name__ == "__main__":
    os.chdir(os.path.dirname(os.path.abspath(__file__)) + "/..")
    main()
