#!/usr/bin/env python3
"""
GPU 8: REMAINING SIMULATIONS (Sequential)
- budesonide_rep3
- natural_cid10592_rep1, rep2, rep3
- apo_rep1 (NEGATIVE CONTROL - no ligand)
- decoy_rep1 (NEGATIVE CONTROL - weak binder)

RTX 5090 + HMR (4fs timestep) - ~30 min each = ~3 hours total

FIXED: OpenFF parameterization, proper NVT→NPT, random seeds
"""

from openmm import *
from openmm.app import *
from openmm.unit import *
import sys
import os
import gc
from datetime import datetime
from openmmforcefields.generators import SystemGenerator
from openff.toolkit import Molecule

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

# Output intervals (at 4fs timestep)
REPORT_INTERVAL = 2500      # 10ps at 4fs
TRAJECTORY_INTERVAL = 1250  # 5ps at 4fs
CHECKPOINT_INTERVAL = 25000 # 100ps at 4fs

# Simulations to run sequentially
# Format: (system_name, replicate, pdb_file, sdf_file_or_None)
SIMULATIONS = [
    ("budesonide", 3, "../structures/complex_budesonide.pdb", "../structures/budesonide_docked.sdf"),
    ("natural_cid10592", 1, "../structures/complex_natural_cid10592.pdb", "../structures/natural_top_docked.sdf"),
    ("natural_cid10592", 2, "../structures/complex_natural_cid10592.pdb", "../structures/natural_top_docked.sdf"),
    ("natural_cid10592", 3, "../structures/complex_natural_cid10592.pdb", "../structures/natural_top_docked.sdf"),
    ("apo", 1, "../structures/complex_apo.pdb", None),  # No ligand
    ("decoy", 1, "../structures/complex_decoy.pdb", "../structures/decoy_docked.sdf"),
]

# ============================================================
# SIMULATION FUNCTION
# ============================================================
def run_simulation(system_name, replicate, pdb_file, sdf_file, sim_num, total_sims):
    """Run a single MD simulation with proper parameterization."""

    random_seed = 8000 + sim_num * 100 + replicate  # Unique seed

    print("\n" + "=" * 70)
    print(f"GPU {GPU_ID} | SIMULATION {sim_num}/{total_sims}")
    print(f"{system_name.upper()} - Replicate {replicate}")
    print(f"Started: {datetime.now().strftime('%H:%M:%S')}")
    print(f"Random seed: {random_seed}")
    print("=" * 70)

    output_prefix = f"{system_name}_rep{replicate}"
    is_apo = sdf_file is None

    # --------------------------------------------------------
    # 1. Load ligand (if present)
    # --------------------------------------------------------
    if not is_apo:
        print("\n[1/7] Loading ligand from SDF...")
        ligand_mol = Molecule.from_file(sdf_file)
        print(f"      Ligand atoms: {ligand_mol.n_atoms}")
    else:
        print("\n[1/7] APO system - no ligand")
        ligand_mol = None

    # --------------------------------------------------------
    # 2. Load complex structure
    # --------------------------------------------------------
    print("\n[2/7] Loading structure...")
    pdb = PDBFile(pdb_file)
    print(f"      Atoms: {pdb.topology.getNumAtoms()}")

    # --------------------------------------------------------
    # 3. Create SystemGenerator
    # --------------------------------------------------------
    print("\n[3/7] Setting up force fields...")

    forcefield_kwargs = {
        'constraints': AllBonds,
        'rigidWater': True,
        'removeCMMotion': True,
        'hydrogenMass': HYDROGEN_MASS if USE_HMR else None
    }

    if not is_apo:
        system_generator = SystemGenerator(
            forcefields=['amber14-all.xml', 'amber14/tip3pfb.xml'],
            small_molecule_forcefield='openff-2.1.0',
            molecules=[ligand_mol],
            forcefield_kwargs=forcefield_kwargs
        )
        print("      Protein: Amber14 | Ligand: OpenFF 2.1.0")
    else:
        # APO system - no ligand, use standard ForceField
        system_generator = SystemGenerator(
            forcefields=['amber14-all.xml', 'amber14/tip3pfb.xml'],
            small_molecule_forcefield='openff-2.1.0',
            molecules=[],
            forcefield_kwargs=forcefield_kwargs
        )
        print("      Protein: Amber14 | No ligand (APO)")

    # --------------------------------------------------------
    # 4. Solvate
    # --------------------------------------------------------
    print("\n[4/7] Adding solvent...")
    modeller = Modeller(pdb.topology, pdb.positions)
    modeller.addSolvent(
        system_generator.forcefield,
        model='tip3p',
        padding=1.0*nanometers,
        ionicStrength=0.15*molar
    )
    print(f"      Solvated atoms: {modeller.topology.getNumAtoms()}")

    # --------------------------------------------------------
    # 5. Create system
    # --------------------------------------------------------
    print("\n[5/7] Creating system (HMR enabled)...")
    system = system_generator.create_system(
        modeller.topology,
        nonbondedMethod=PME,
        nonbondedCutoff=1.0*nanometers
    )

    # --------------------------------------------------------
    # 6. Setup simulation
    # --------------------------------------------------------
    print(f"\n[6/7] Setting up CUDA (GPU {GPU_ID})...")
    integrator = LangevinMiddleIntegrator(TEMPERATURE, FRICTION, TIMESTEP)
    integrator.setRandomNumberSeed(random_seed)

    platform = Platform.getPlatformByName('CUDA')
    properties = {'Precision': 'mixed', 'DeviceIndex': str(GPU_ID)}

    simulation = Simulation(modeller.topology, system, integrator, platform, properties)
    simulation.context.setPositions(modeller.positions)

    # --------------------------------------------------------
    # MINIMIZATION
    # --------------------------------------------------------
    print("\n>>> MINIMIZATION")
    state = simulation.context.getState(getEnergy=True)
    print(f"    Initial: {state.getPotentialEnergy()}")
    simulation.minimizeEnergy(maxIterations=MINIMIZATION_STEPS)
    state = simulation.context.getState(getEnergy=True)
    print(f"    Final: {state.getPotentialEnergy()}")

    # --------------------------------------------------------
    # NVT EQUILIBRATION (true NVT - no barostat yet)
    # --------------------------------------------------------
    print(f"\n>>> NVT EQUILIBRATION ({EQUILIBRATION_NVT_PS} ps)")
    simulation.context.setVelocitiesToTemperature(TEMPERATURE, random_seed)
    nvt_steps = int(EQUILIBRATION_NVT_PS * 1000 / TIMESTEP.value_in_unit(femtoseconds))
    simulation.step(nvt_steps)

    # --------------------------------------------------------
    # ADD BAROSTAT FOR NPT
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

    print(f"\n>>> COMPLETE: {output_prefix}")
    print(f"    Finished: {datetime.now().strftime('%H:%M:%S')}")

    # --------------------------------------------------------
    # CLEANUP (prevent memory leaks between simulations)
    # --------------------------------------------------------
    del simulation
    del system
    del integrator
    del modeller
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

    for i, (name, rep, pdb, sdf) in enumerate(SIMULATIONS, 1):
        run_simulation(name, rep, pdb, sdf, i, len(SIMULATIONS))

    print("\n" + "=" * 70)
    print("ALL GPU 8 SIMULATIONS COMPLETE!")
    print(f"Finished: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("=" * 70)

if __name__ == "__main__":
    os.chdir(os.path.dirname(os.path.abspath(__file__)) + "/..")
    main()
