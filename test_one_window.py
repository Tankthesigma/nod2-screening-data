#!/usr/bin/env python
"""
Quick test of one FEP window to verify setup works.
Runs minimal steps just to confirm everything loads and runs.
"""
import os
import sys
import numpy as np
import time

from openmm import XmlSerializer, LangevinMiddleIntegrator, Platform, MonteCarloBarostat
from openmm.app import PDBFile, Simulation
from openmm import unit

# Test configuration
SYS_NAME = "wt_complex"  # or "mut_complex" or "solvent"
WINDOW_IDX = 0
TEST_STEPS = 1000  # Just 2 ps for quick test

def main():
    print("="*60)
    print("FEP WINDOW TEST")
    print("="*60)

    sys_dir = f"fep_pmx/{SYS_NAME}"

    # Check files exist
    print(f"\n[1] Checking files for {SYS_NAME}...")
    required_files = [
        "alchemical_system.xml",
        "topology.pdb",
        "positions.npy",
        "lambda_schedule.npy"
    ]
    for f in required_files:
        path = os.path.join(sys_dir, f)
        if os.path.exists(path):
            print(f"  [PASS] {f}")
        else:
            print(f"  [FAIL] {f} not found!")
            sys.exit(1)

    # Load system
    print(f"\n[2] Loading alchemical system...")
    t0 = time.time()
    with open(os.path.join(sys_dir, "alchemical_system.xml"), "r") as f:
        system = XmlSerializer.deserialize(f.read())
    print(f"  Loaded in {time.time()-t0:.1f}s")
    print(f"  Forces: {system.getNumForces()}")
    print(f"  Particles: {system.getNumParticles()}")

    # Check for barostat
    has_barostat = False
    for i in range(system.getNumForces()):
        force = system.getForce(i)
        if isinstance(force, MonteCarloBarostat):
            has_barostat = True
            print(f"  [PASS] Barostat found")
            break
    if not has_barostat:
        print(f"  [FAIL] No barostat!")
        sys.exit(1)

    # Load positions
    print(f"\n[3] Loading positions...")
    positions = np.load(os.path.join(sys_dir, "positions.npy"))
    positions = positions * unit.nanometer
    print(f"  Shape: {positions.shape}")

    # Load topology
    pdb = PDBFile(os.path.join(sys_dir, "topology.pdb"))
    topology = pdb.topology

    # Load lambda schedule
    lambda_schedule = np.load(os.path.join(sys_dir, "lambda_schedule.npy"))
    n_windows = len(lambda_schedule)
    print(f"  Lambda windows: {n_windows}")

    # Get lambda for this window
    lam_e, lam_s, lam_r = lambda_schedule[WINDOW_IDX]
    print(f"\n[4] Window {WINDOW_IDX} lambda values:")
    print(f"  electrostatics: {lam_e:.4f}")
    print(f"  sterics: {lam_s:.4f}")
    print(f"  restraints: {lam_r:.4f}")

    # Create integrator
    print(f"\n[5] Creating simulation...")
    integrator = LangevinMiddleIntegrator(
        310 * unit.kelvin,
        1.0 / unit.picosecond,
        2.0 * unit.femtoseconds
    )

    # Select platform - prefer CUDA
    try:
        platform = Platform.getPlatformByName('CUDA')
        properties = {'Precision': 'mixed'}
        print(f"  Platform: CUDA")
    except Exception as e:
        print(f"  CUDA failed: {e}")
        platform = Platform.getPlatformByName('CPU')
        properties = {}
        print(f"  Platform: CPU (fallback)")

    # Create simulation
    simulation = Simulation(topology, system, integrator, platform, properties)
    simulation.context.setPositions(positions)

    # Set lambda parameters
    context = simulation.context
    context.setParameter('lambda_electrostatics', lam_e)
    context.setParameter('lambda_sterics', lam_s)
    if SYS_NAME != 'solvent':
        context.setParameter('lambda_restraints', lam_r)
    print(f"  [PASS] Lambda parameters set")

    # Get initial energy
    state = context.getState(getEnergy=True)
    E0 = state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
    print(f"  Initial energy: {E0:.1f} kJ/mol")

    # Minimize
    print(f"\n[6] Minimizing...")
    t0 = time.time()
    simulation.minimizeEnergy(maxIterations=100)
    state = context.getState(getEnergy=True)
    E1 = state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
    print(f"  Minimized in {time.time()-t0:.1f}s")
    print(f"  Energy: {E0:.1f} -> {E1:.1f} kJ/mol")

    # Run short test
    print(f"\n[7] Running {TEST_STEPS} steps...")
    t0 = time.time()
    simulation.step(TEST_STEPS)
    elapsed = time.time() - t0
    ns_per_day = (TEST_STEPS * 0.002 / 1000) / (elapsed / 86400)
    print(f"  Completed in {elapsed:.1f}s")
    print(f"  Performance: ~{ns_per_day:.1f} ns/day")

    # Test u_nk computation
    print(f"\n[8] Testing u_nk computation...")
    kT = unit.MOLAR_GAS_CONSTANT_R * 310 * unit.kelvin
    u_nk_test = np.zeros(n_windows)

    for k in range(n_windows):
        lam_e_k, lam_s_k, lam_r_k = lambda_schedule[k]
        context.setParameter('lambda_electrostatics', lam_e_k)
        context.setParameter('lambda_sterics', lam_s_k)
        if SYS_NAME != 'solvent':
            context.setParameter('lambda_restraints', lam_r_k)

        state = context.getState(getEnergy=True)
        u_nk_test[k] = state.getPotentialEnergy() / kT

    # Reset lambda
    context.setParameter('lambda_electrostatics', lam_e)
    context.setParameter('lambda_sterics', lam_s)
    if SYS_NAME != 'solvent':
        context.setParameter('lambda_restraints', lam_r)

    print(f"  u_nk shape: ({n_windows},)")
    print(f"  u_nk range: [{u_nk_test.min():.1f}, {u_nk_test.max():.1f}]")
    print(f"  [PASS] u_nk computation works")

    print()
    print("="*60)
    print("TEST PASSED! Ready for production runs.")
    print("="*60)
    print()
    print("Next steps:")
    print(f"  1. Run all 20 windows for {SYS_NAME}")
    print(f"  2. Run all 20 windows for other systems")
    print(f"  3. Collect u_nk.npy files and run MBAR")
    print()
    print("To run one window:")
    print(f"  cd fep_pmx/{SYS_NAME}/window_00")
    print(f"  python run_window.py")

if __name__ == "__main__":
    main()
