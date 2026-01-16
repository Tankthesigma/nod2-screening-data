#!/usr/bin/env python
"""Quick canary test - run window 15 briefly to check stability."""
import os
import numpy as np
from pathlib import Path
from openmm import *
from openmm.app import *
from openmm import unit

BASE_DIR = Path("C:/Users/vasud/nod2-screening-data/fep_pmx_natural")
SYS_NAME = "wt_complex"

# Lambda values for window 15
LAMBDA_ELECTROSTATICS = 0.0
LAMBDA_STERICS = 0.4

# Short test - 1000 steps = 2 ps (CPU is slower)
TEST_STEPS = 1000
TIMESTEP = 2.0 * unit.femtoseconds
TEMPERATURE = 310.0 * unit.kelvin

def main():
    print("="*60)
    print("QUICK CANARY TEST - Window 15")
    print("="*60)

    sys_dir = BASE_DIR / SYS_NAME

    # Load system
    print("\nLoading alchemical system...")
    with open(sys_dir / "alchemical_system.xml", 'r') as f:
        system = XmlSerializer.deserialize(f.read())

    # Load topology
    print("Loading topology...")
    pdb = PDBFile(str(sys_dir / "topology.pdb"))

    # Load positions
    print("Loading positions...")
    positions_nm = np.load(sys_dir / "positions.npy")
    positions = [Vec3(x, y, z) * unit.nanometers for x, y, z in positions_nm]

    # Set lambda values
    print(f"\nSetting lambda: electrostatics={LAMBDA_ELECTROSTATICS}, sterics={LAMBDA_STERICS}")
    for i in range(system.getNumForces()):
        force = system.getForce(i)
        if hasattr(force, 'getNumGlobalParameters'):
            for j in range(force.getNumGlobalParameters()):
                name = force.getGlobalParameterName(j)
                if name == 'lambda_electrostatics':
                    force.setGlobalParameterDefaultValue(j, LAMBDA_ELECTROSTATICS)
                elif name == 'lambda_sterics':
                    force.setGlobalParameterDefaultValue(j, LAMBDA_STERICS)

    # Create integrator
    integrator = LangevinMiddleIntegrator(TEMPERATURE, 1.0/unit.picoseconds, TIMESTEP)

    # Use CPU platform (CUDA PTX issue on Windows)
    print("\nInitializing simulation...")
    platform = Platform.getPlatformByName('CPU')
    properties = {}
    print("  Using CPU platform")

    simulation = Simulation(pdb.topology, system, integrator, platform, properties)
    simulation.context.setPositions(positions)

    # Minimize
    print("\nMinimizing energy...")
    state = simulation.context.getState(getEnergy=True)
    print(f"  Initial energy: {state.getPotentialEnergy()}")

    simulation.minimizeEnergy(maxIterations=1000)

    state = simulation.context.getState(getEnergy=True)
    print(f"  After minimization: {state.getPotentialEnergy()}")

    # Set velocities
    simulation.context.setVelocitiesToTemperature(TEMPERATURE)

    # Run short test
    print(f"\nRunning {TEST_STEPS} steps ({TEST_STEPS * 0.002} ps)...")
    import time
    t0 = time.time()

    for i in range(5):
        simulation.step(TEST_STEPS // 5)
        state = simulation.context.getState(getEnergy=True, getPositions=True)
        pe = state.getPotentialEnergy()
        print(f"  Step {(i+1)*TEST_STEPS//5}: PE = {pe}")

        # Check for NaN
        if np.isnan(pe.value_in_unit(unit.kilojoules_per_mole)):
            print("  [FAIL] NaN energy detected!")
            return False

    elapsed = time.time() - t0
    print(f"\nCompleted in {elapsed:.1f}s ({TEST_STEPS/elapsed:.0f} steps/s)")

    # Final check
    state = simulation.context.getState(getEnergy=True)
    final_pe = state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)

    if np.isnan(final_pe) or abs(final_pe) > 1e10:
        print("\n[FAIL] Unstable simulation!")
        return False

    print("\n" + "="*60)
    print("[PASS] Quick canary test completed successfully")
    print("="*60)
    return True

if __name__ == "__main__":
    success = main()
    exit(0 if success else 1)
