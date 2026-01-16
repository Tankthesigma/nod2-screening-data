#!/usr/bin/env python
"""Simple canary test for windows 15-19 - starts from initial positions."""
import numpy as np
import time
from pathlib import Path
from openmm import *
from openmm.app import *
from openmm import unit

BASE_DIR = Path("C:/Users/vasud/nod2-screening-data/fep_pmx_natural")

# Lambda schedule (stored in wt_complex folder from febuxostat)
LAMBDA_SCHEDULE = np.load(BASE_DIR / "wt_complex" / "lambda_schedule.npy")

# Simulation parameters (reduced for canary testing)
EQUIL_STEPS = 5000     # 10 ps equilibration
PROD_STEPS = 10000     # 20 ps production (quick test)
TIMESTEP = 2.0 * unit.femtoseconds
TEMPERATURE = 310.0 * unit.kelvin


def run_window_canary(sys_name, window_idx):
    """Run a single window canary test."""
    sys_dir = BASE_DIR / sys_name

    # Get lambda values
    lambda_elec = float(LAMBDA_SCHEDULE[window_idx, 0])
    lambda_sterics = float(LAMBDA_SCHEDULE[window_idx, 1])

    print(f"\n  Window {window_idx}: lambda_elec={lambda_elec:.2f}, lambda_sterics={lambda_sterics:.2f}")

    # Load system
    with open(sys_dir / "alchemical_system.xml", 'r') as f:
        system = XmlSerializer.deserialize(f.read())

    # Load topology and positions
    pdb = PDBFile(str(sys_dir / "topology.pdb"))
    positions_nm = np.load(sys_dir / "positions.npy")
    positions = [Vec3(x, y, z) * unit.nanometers for x, y, z in positions_nm]

    # Set lambda values in system
    for i in range(system.getNumForces()):
        force = system.getForce(i)
        if hasattr(force, 'getNumGlobalParameters'):
            for j in range(force.getNumGlobalParameters()):
                name = force.getGlobalParameterName(j)
                if name == 'lambda_electrostatics':
                    force.setGlobalParameterDefaultValue(j, lambda_elec)
                elif name == 'lambda_sterics':
                    force.setGlobalParameterDefaultValue(j, lambda_sterics)

    # Create integrator and simulation
    integrator = LangevinMiddleIntegrator(TEMPERATURE, 1.0/unit.picoseconds, TIMESTEP)
    platform = Platform.getPlatformByName('CPU')

    simulation = Simulation(pdb.topology, system, integrator, platform, {})
    simulation.context.setPositions(positions)

    # Also set lambda in context
    simulation.context.setParameter('lambda_electrostatics', lambda_elec)
    simulation.context.setParameter('lambda_sterics', lambda_sterics)

    # Minimization
    print(f"    Minimizing...", end='', flush=True)
    simulation.minimizeEnergy(maxIterations=500)

    state = simulation.context.getState(getEnergy=True)
    min_energy = state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
    print(f" {min_energy:.0f} kJ/mol")

    if np.isnan(min_energy) or abs(min_energy) > 1e10:
        print(f"    [FAIL] Bad energy after minimization!")
        return False

    # Set velocities
    simulation.context.setVelocitiesToTemperature(TEMPERATURE)

    # Equilibration
    print(f"    Equilibrating ({EQUIL_STEPS} steps)...", end='', flush=True)
    t0 = time.time()

    try:
        simulation.step(EQUIL_STEPS)
    except Exception as e:
        print(f"\n    [FAIL] {e}")
        return False

    print(f" done ({time.time()-t0:.1f}s)")

    state = simulation.context.getState(getEnergy=True)
    eq_energy = state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)

    if np.isnan(eq_energy):
        print(f"    [FAIL] NaN energy after equilibration!")
        return False

    # Production
    print(f"    Production ({PROD_STEPS} steps)...", end='', flush=True)
    t0 = time.time()

    energies = []
    n_chunks = 5
    chunk_size = PROD_STEPS // n_chunks

    for i in range(n_chunks):
        try:
            simulation.step(chunk_size)
        except Exception as e:
            print(f"\n    [FAIL] {e}")
            return False

        state = simulation.context.getState(getEnergy=True)
        pe = state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
        energies.append(pe)

        if np.isnan(pe):
            print(f"\n    [FAIL] NaN energy during production!")
            return False

    elapsed = time.time() - t0
    print(f" done ({elapsed:.1f}s)")

    avg_energy = np.mean(energies)
    std_energy = np.std(energies)
    print(f"    Final: PE = {avg_energy:.0f} +/- {std_energy:.0f} kJ/mol")

    if std_energy > 1e6:
        print(f"    [FAIL] Unstable simulation!")
        return False

    print(f"    [PASS]")
    return True


def main():
    print("="*60)
    print("CANARY TESTS - Windows 15-19 (Simple)")
    print("="*60)

    results = {}

    for sys_name in ['wt_complex', 'mut_complex']:
        print(f"\n{'='*60}")
        print(f"Testing {sys_name}")
        print("="*60)

        results[sys_name] = {}

        for window_idx in range(15, 20):
            try:
                success = run_window_canary(sys_name, window_idx)
                results[sys_name][window_idx] = success
            except Exception as e:
                print(f"    [ERROR] {e}")
                results[sys_name][window_idx] = False

    # Summary
    print("\n" + "="*60)
    print("CANARY TEST SUMMARY")
    print("="*60)

    all_pass = True
    for sys_name in ['wt_complex', 'mut_complex']:
        print(f"\n{sys_name}:")
        for window_idx in range(15, 20):
            status = "PASS" if results[sys_name][window_idx] else "FAIL"
            print(f"  Window {window_idx}: {status}")
            if not results[sys_name][window_idx]:
                all_pass = False

    print("\n" + "="*60)
    if all_pass:
        print("[ALL CANARY TESTS PASSED]")
    else:
        print("[SOME CANARY TESTS FAILED]")
    print("="*60)

    return all_pass


if __name__ == "__main__":
    success = main()
    exit(0 if success else 1)
