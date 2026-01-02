#!/usr/bin/env python3
"""
MD Simulation Script for ursodiol
Run on GPU instance with OpenMM installed

Protocol:
- Energy minimization: 5000 steps
- NVT equilibration: 100 ps
- NPT equilibration: 500 ps
- Production: 3 x 20 ns

Estimated time on RTX 4090: ~2-4 hours per system
"""

from openmm import *
from openmm.app import *
from openmm.unit import *
import sys

# Configuration
SYSTEM_NAME = "ursodiol"
TEMPERATURE = 310.15  # K
PRESSURE = 1.0  # atm
TIMESTEP = 2.0  # fs
FRICTION = 1.0  # 1/ps

MINIMIZATION_STEPS = 5000
EQUILIBRATION_NVT_PS = 100
EQUILIBRATION_NPT_PS = 500
PRODUCTION_NS = 20
N_REPLICATES = 3

REPORT_INTERVAL = 10000
TRAJECTORY_INTERVAL = 5000
CHECKPOINT_INTERVAL = 50000


def run_simulation(complex_pdb, replicate):
    """Run full MD simulation"""

    print(f"\n============================================================")
    print(f"MD SIMULATION: {SYSTEM_NAME} - Replicate {replicate+1}")
    print(f"============================================================\n")

    # Load structure
    print("Loading structure...")
    pdb = PDBFile(complex_pdb)

    # Force field
    print("Setting up force field...")
    forcefield = ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')

    # Create system
    print("Creating system...")
    system = forcefield.createSystem(
        pdb.topology,
        nonbondedMethod=PME,
        nonbondedCutoff=1.0*nanometers,
        constraints=HBonds
    )

    # Integrator
    integrator = LangevinMiddleIntegrator(
        TEMPERATURE*kelvin,
        FRICTION/picoseconds,
        TIMESTEP*femtoseconds
    )

    # Platform (CUDA preferred)
    print("Selecting platform...")
    try:
        platform = Platform.getPlatformByName('CUDA')
        properties = {'CudaPrecision': 'mixed'}
        print("  Using CUDA platform")
    except:
        try:
            platform = Platform.getPlatformByName('OpenCL')
            properties = {}
            print("  Using OpenCL platform")
        except:
            platform = Platform.getPlatformByName('CPU')
            properties = {}
            print("  WARNING: Using CPU platform (slow!)")

    # Create simulation
    simulation = Simulation(pdb.topology, system, integrator, platform, properties)
    simulation.context.setPositions(pdb.positions)

    # MINIMIZATION
    print("\nEnergy minimization...")
    state = simulation.context.getState(getEnergy=True)
    print(f"  Initial energy: {state.getPotentialEnergy()}")

    simulation.minimizeEnergy(maxIterations=MINIMIZATION_STEPS)

    state = simulation.context.getState(getEnergy=True)
    print(f"  Final energy: {state.getPotentialEnergy()}")

    # NVT EQUILIBRATION
    print(f"\nNVT equilibration ({EQUILIBRATION_NVT_PS} ps)...")
    simulation.context.setVelocitiesToTemperature(TEMPERATURE*kelvin)

    nvt_steps = int(EQUILIBRATION_NVT_PS * 1000 / TIMESTEP)
    simulation.step(nvt_steps)
    print("  NVT complete")

    # Add barostat for NPT
    system.addForce(MonteCarloBarostat(PRESSURE*atmospheres, TEMPERATURE*kelvin))
    simulation.context.reinitialize(preserveState=True)

    # NPT EQUILIBRATION
    print(f"\nNPT equilibration ({EQUILIBRATION_NPT_PS} ps)...")
    npt_steps = int(EQUILIBRATION_NPT_PS * 1000 / TIMESTEP)
    simulation.step(npt_steps)
    print("  NPT complete")

    # PRODUCTION
    print(f"\nProduction MD ({PRODUCTION_NS} ns)...")

    # Output files
    traj_file = f"{SYSTEM_NAME}_rep{replicate+1}.dcd"
    log_file = f"{SYSTEM_NAME}_rep{replicate+1}.log"
    chk_file = f"{SYSTEM_NAME}_rep{replicate+1}.chk"

    production_steps = int(PRODUCTION_NS * 1e6 / TIMESTEP)

    # Add reporters
    simulation.reporters.append(DCDReporter(traj_file, TRAJECTORY_INTERVAL))
    simulation.reporters.append(StateDataReporter(
        log_file, REPORT_INTERVAL,
        step=True, time=True,
        potentialEnergy=True, kineticEnergy=True,
        temperature=True, progress=True,
        remainingTime=True, speed=True,
        totalSteps=production_steps
    ))
    simulation.reporters.append(CheckpointReporter(chk_file, CHECKPOINT_INTERVAL))

    # Run production
    simulation.step(production_steps)

    print(f"\nProduction complete!")
    print(f"  Trajectory: {traj_file}")
    print(f"  Log: {log_file}")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Run MD simulation")
    parser.add_argument("--pdb", required=True, help="Input complex PDB file")
    parser.add_argument("--replicate", type=int, default=0, help="Replicate number (0-indexed)")

    args = parser.parse_args()

    run_simulation(args.pdb, args.replicate)
    print("\nSimulation finished successfully!")
