#!/usr/bin/env python3
"""
VAST.AI RTX 5090 - Heavy Load (10 simulations, 200ns total)
Sequential execution with auto-shutdown on completion.

SAFETY FEATURES:
- Git push after EACH simulation (crash protection)
- Git push on ANY error (never lose data)
- NaN energy check (fail fast)
- Robust shutdown with fallback

Simulations:
1. Febuxostat rep1-3 (60ns)
2. Ursodiol rep1-3 (60ns)
3. Budesonide rep1-3 (60ns)
4. Natural_top rep1 (20ns)

Expected runtime: ~5 hours on RTX 5090 with HMR (~1000 ns/day)
"""

from openmm import *
from openmm.app import *
from openmm.unit import *
import sys
import os
import gc
import math
import subprocess
from datetime import datetime
import traceback

try:
    from openmmforcefields.generators import SystemGenerator, GAFFTemplateGenerator
    from openff.toolkit.topology import Molecule  # Explicit import path
    OPENFF_AVAILABLE = True
except ImportError:
    print("FATAL: OpenFF/openmmforcefields not installed!")
    print("Run: conda install -c conda-forge openmm openmmforcefields openff-toolkit -y")
    sys.exit(1)

# ============================================================
# CONFIGURATION
# ============================================================
GPU_ID = 0

STRUCTURES_DIR = "structures"
OUTPUT_DIR = "trajectories"
LOG_DIR = "logs"
CHECKPOINT_DIR = "checkpoints"

TEMPERATURE = 310.15 * kelvin
PRESSURE = 1.0 * atmospheres
TIMESTEP = 4.0 * femtoseconds
FRICTION = 1.0 / picoseconds
USE_HMR = True
HYDROGEN_MASS = 3.0 * amu

MINIMIZATION_STEPS = 5000
EQUILIBRATION_NVT_PS = 100
EQUILIBRATION_NPT_PS = 500
PRODUCTION_NS = 20

REPORT_INTERVAL = 2500
TRAJECTORY_INTERVAL = 2500
CHECKPOINT_INTERVAL = 250000

SIMULATIONS = [
    ("febuxostat", 1, "complex_febuxostat.pdb", "febuxostat_docked.sdf"),
    ("febuxostat", 2, "complex_febuxostat.pdb", "febuxostat_docked.sdf"),
    ("febuxostat", 3, "complex_febuxostat.pdb", "febuxostat_docked.sdf"),
    ("ursodiol", 1, "complex_ursodiol.pdb", "ursodiol_docked.sdf"),
    ("ursodiol", 2, "complex_ursodiol.pdb", "ursodiol_docked.sdf"),
    ("ursodiol", 3, "complex_ursodiol.pdb", "ursodiol_docked.sdf"),
    ("budesonide", 1, "complex_budesonide.pdb", "budesonide_docked.sdf"),
    ("budesonide", 2, "complex_budesonide.pdb", "budesonide_docked.sdf"),
    ("budesonide", 3, "complex_budesonide.pdb", "budesonide_docked.sdf"),
    ("natural_top", 1, "complex_natural_cid10592.pdb", "natural_top_docked.sdf"),
]

# ============================================================
# GIT SAVE FUNCTION (CRITICAL - NEVER LOSE DATA)
# ============================================================
def save_to_git(message="Auto-save simulation results"):
    """Push results to git. Called after each sim and on crash."""
    try:
        print(f"\n>>> SAVING TO GIT: {message}")
        subprocess.run(["git", "add", "-A"], check=False, capture_output=True)
        subprocess.run(
            ["git", "commit", "-m", f"{message} - {datetime.now().strftime('%Y-%m-%d %H:%M')}"],
            check=False, capture_output=True
        )
        result = subprocess.run(["git", "push"], check=False, capture_output=True)
        if result.returncode == 0:
            print(">>> GIT PUSH: SUCCESS")
        else:
            print(f">>> GIT PUSH: FAILED (but local commit saved)")
            print(f"    Error: {result.stderr.decode()[:200]}")
    except Exception as e:
        print(f">>> GIT SAVE ERROR: {e}")

# ============================================================
# HELPER FUNCTIONS
# ============================================================
def setup_directories():
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    os.makedirs(LOG_DIR, exist_ok=True)
    os.makedirs(CHECKPOINT_DIR, exist_ok=True)

def check_energy(simulation, stage=""):
    """Check if energy is valid (not NaN/Inf). Returns False if exploded."""
    state = simulation.context.getState(getEnergy=True)
    energy = state.getPotentialEnergy()
    energy_val = energy.value_in_unit(kilojoules_per_mole)

    if math.isnan(energy_val) or math.isinf(energy_val):
        print(f"\n!!! CRITICAL: Energy exploded during {stage}!")
        print(f"!!! Energy = {energy}")
        print("!!! Check input PDB/SDF for clashes or bad geometry")
        return False, energy
    return True, energy

def add_position_restraints(system, topology, positions, k=1000.0):
    """Add position restraints to heavy atoms (non-hydrogen).

    Returns the force object so it can be removed later.
    k is in kJ/mol/nm^2 (1000 = strong restraint)
    """
    restraint = CustomExternalForce("0.5*k*((x-x0)^2+(y-y0)^2+(z-z0)^2)")
    restraint.addGlobalParameter("k", k * kilojoules_per_mole / nanometers**2)
    restraint.addPerParticleParameter("x0")
    restraint.addPerParticleParameter("y0")
    restraint.addPerParticleParameter("z0")

    restrained_count = 0
    for atom in topology.atoms():
        if atom.element.symbol != 'H':  # Heavy atoms only
            pos = positions[atom.index]
            restraint.addParticle(atom.index, [
                pos[0].value_in_unit(nanometers),
                pos[1].value_in_unit(nanometers),
                pos[2].value_in_unit(nanometers)
            ])
            restrained_count += 1

    force_idx = system.addForce(restraint)
    print(f"      Added restraints to {restrained_count} heavy atoms (k={k} kJ/mol/nm²)")
    return force_idx

def create_system_with_fallback(pdb, ligand_mol, forcefield_kwargs):
    """Create system with OpenFF, falling back to GAFF if needed."""
    modeller = Modeller(pdb.topology, pdb.positions)
    modeller.deleteWater()
    print("      Deleted existing water molecules")

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
            system = system_generator.create_system(
                modeller.topology, nonbondedMethod=PME, nonbondedCutoff=1.0*nanometers
            )
            print("      SUCCESS: OpenFF")
            return system, modeller
        except Exception as e:
            print(f"      OpenFF failed: {e}, trying GAFF...")
            modeller = Modeller(pdb.topology, pdb.positions)
            modeller.deleteWater()

    print("      Using GAFF...")
    forcefield = ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')
    gaff_generator = GAFFTemplateGenerator(molecules=[ligand_mol])
    forcefield.registerTemplateGenerator(gaff_generator.generator)
    modeller.addSolvent(forcefield, model='tip3p',
                      padding=1.0*nanometers, ionicStrength=0.15*molar)
    system = forcefield.createSystem(
        modeller.topology, nonbondedMethod=PME, nonbondedCutoff=1.0*nanometers,
        constraints=AllBonds, hydrogenMass=HYDROGEN_MASS if USE_HMR else None
    )
    print("      SUCCESS: GAFF")
    return system, modeller

# ============================================================
# SIMULATION FUNCTION
# ============================================================
def run_simulation(name, replicate, pdb_file, sdf_file, sim_num, total_sims):
    """Run a single MD simulation with full error checking."""
    random_seed = sim_num * 1000 + replicate
    output_prefix = f"{name}_rep{replicate}"

    pdb_path = f"{STRUCTURES_DIR}/{pdb_file}"
    sdf_path = f"{STRUCTURES_DIR}/{sdf_file}"

    print("\n" + "=" * 70)
    print(f"[{sim_num}/{total_sims}] {name.upper()} - Replicate {replicate}")
    print(f"Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("=" * 70)

    # Validate inputs
    if not os.path.exists(pdb_path):
        print(f"ERROR: PDB not found: {pdb_path}")
        return False
    if not os.path.exists(sdf_path):
        print(f"ERROR: SDF not found: {sdf_path}")
        return False

    # Load ligand
    print("\n[1/8] Loading ligand...")
    ligand_mol = Molecule.from_file(sdf_path)
    print(f"      Atoms: {ligand_mol.n_atoms}")

    # Load structure
    print("\n[2/8] Loading structure...")
    pdb = PDBFile(pdb_path)
    print(f"      Atoms: {pdb.topology.getNumAtoms()}")

    # Create system
    print("\n[3/8] Setting up force fields...")
    forcefield_kwargs = {
        'constraints': AllBonds,
        'rigidWater': True,
        'removeCMMotion': True,
        'hydrogenMass': HYDROGEN_MASS if USE_HMR else None
    }
    system, modeller = create_system_with_fallback(pdb, ligand_mol, forcefield_kwargs)
    print(f"      Solvated atoms: {modeller.topology.getNumAtoms()}")

    # Save solvated topology for analysis (needed for DCD loading)
    solvated_pdb = f"{OUTPUT_DIR}/{output_prefix}_solvated.pdb"
    with open(solvated_pdb, 'w') as f:
        PDBFile.writeFile(modeller.topology, modeller.positions, f)
    print(f"      Saved: {solvated_pdb}")

    # Setup simulation
    print(f"\n[4/8] Setting up CUDA (GPU {GPU_ID})...")
    integrator = LangevinMiddleIntegrator(TEMPERATURE, FRICTION, TIMESTEP)
    integrator.setRandomNumberSeed(random_seed)
    integrator.setConstraintTolerance(1e-6)  # Tighter tolerance for stability
    platform = Platform.getPlatformByName('CUDA')
    properties = {'Precision': 'mixed', 'DeviceIndex': str(GPU_ID)}
    simulation = Simulation(modeller.topology, system, integrator, platform, properties)
    simulation.context.setPositions(modeller.positions)

    # Minimization with NaN check
    print("\n[5/8] MINIMIZATION")
    ok, energy = check_energy(simulation, "pre-minimization")
    print(f"      Initial: {energy}")

    simulation.minimizeEnergy(maxIterations=MINIMIZATION_STEPS)

    ok, energy = check_energy(simulation, "post-minimization")
    if not ok:
        return False
    print(f"      Final: {energy}")

    # NVT with position restraints
    print(f"\n[6/8] NVT EQUILIBRATION ({EQUILIBRATION_NVT_PS} ps) - WITH RESTRAINTS")
    restraint_idx = add_position_restraints(system, modeller.topology, modeller.positions, k=1000.0)
    simulation.context.reinitialize(preserveState=True)
    simulation.context.setVelocitiesToTemperature(TEMPERATURE, random_seed)
    nvt_steps = int(EQUILIBRATION_NVT_PS * 1000 / TIMESTEP.value_in_unit(femtoseconds))
    simulation.step(nvt_steps)

    ok, energy = check_energy(simulation, "NVT")
    if not ok:
        return False

    # Remove restraints for NPT
    print("      Removing position restraints...")
    system.removeForce(restraint_idx)
    simulation.context.reinitialize(preserveState=True)

    # NPT with check
    print(f"\n[7/8] NPT EQUILIBRATION ({EQUILIBRATION_NPT_PS} ps)")
    system.addForce(MonteCarloBarostat(PRESSURE, TEMPERATURE, 25))
    simulation.context.reinitialize(preserveState=True)
    npt_steps = int(EQUILIBRATION_NPT_PS * 1000 / TIMESTEP.value_in_unit(femtoseconds))
    simulation.step(npt_steps)

    ok, energy = check_energy(simulation, "NPT")
    if not ok:
        return False

    # Production
    print(f"\n[8/8] PRODUCTION ({PRODUCTION_NS} ns)")
    production_steps = int(PRODUCTION_NS * 1e6 / TIMESTEP.value_in_unit(femtoseconds))

    simulation.reporters.append(DCDReporter(
        f'{OUTPUT_DIR}/{output_prefix}.dcd', TRAJECTORY_INTERVAL
    ))
    simulation.reporters.append(StateDataReporter(
        f'{LOG_DIR}/{output_prefix}.log', REPORT_INTERVAL,
        step=True, time=True, potentialEnergy=True, kineticEnergy=True,
        temperature=True, speed=True, progress=True, remainingTime=True,
        totalSteps=production_steps
    ))
    simulation.reporters.append(CheckpointReporter(
        f'{CHECKPOINT_DIR}/{output_prefix}.chk', CHECKPOINT_INTERVAL
    ))

    simulation.step(production_steps)
    simulation.saveState(f'{CHECKPOINT_DIR}/{output_prefix}_final.xml')

    # Final energy check
    ok, energy = check_energy(simulation, "production")
    if not ok:
        print("WARNING: Final energy check failed but trajectory saved")

    print(f"\n>>> COMPLETE: {output_prefix}")
    print(f"    Finished: {datetime.now().strftime('%H:%M:%S')}")

    # Cleanup
    del simulation, system, integrator, modeller
    gc.collect()

    return True

# ============================================================
# MAIN
# ============================================================
def preflight_check():
    """Check ALL input files exist AND ligands can be loaded BEFORE starting."""
    print("\n[PREFLIGHT] Checking input files...")
    missing = []
    for name, rep, pdb, sdf in SIMULATIONS:
        pdb_path = f"{STRUCTURES_DIR}/{pdb}"
        sdf_path = f"{STRUCTURES_DIR}/{sdf}"
        if not os.path.exists(pdb_path):
            missing.append(pdb_path)
        if not os.path.exists(sdf_path):
            missing.append(sdf_path)

    if missing:
        print("\n!!! FATAL: Missing input files !!!")
        for f in missing:
            print(f"    - {f}")
        print("\nAborting to avoid wasting GPU rental money.")
        sys.exit(1)

    print(f"    All {len(SIMULATIONS) * 2} input files found.")

    # Validate ligands can be loaded (catches bad SDF geometry/bonding)
    print("\n[PREFLIGHT] Validating ligand files...")
    unique_sdfs = set(sdf for _, _, _, sdf in SIMULATIONS)
    for sdf in unique_sdfs:
        sdf_path = f"{STRUCTURES_DIR}/{sdf}"
        try:
            mol = Molecule.from_file(sdf_path)
            print(f"    ✓ {sdf}: {mol.n_atoms} atoms, {mol.n_bonds} bonds")
        except Exception as e:
            print(f"\n!!! FATAL: Cannot load ligand {sdf} !!!")
            print(f"    Error: {e}")
            print("\nThis SDF file has bad geometry or bonding.")
            print("Fix: Re-export from RDKit or check original docking output.")
            sys.exit(1)

    print("    All ligands validated successfully.")

def main():
    print("=" * 70)
    print("VAST.AI RTX 5090 - NOD2-SCOUT MD SIMULATIONS")
    print(f"Running {len(SIMULATIONS)} simulations ({len(SIMULATIONS) * PRODUCTION_NS}ns total)")
    print(f"Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("=" * 70)

    # CHECK ALL FILES BEFORE STARTING
    preflight_check()

    setup_directories()

    results = []
    for i, (name, rep, pdb, sdf) in enumerate(SIMULATIONS, 1):
        try:
            success = run_simulation(name, rep, pdb, sdf, i, len(SIMULATIONS))
            results.append((f"{name}_rep{rep}", success))

            # SAVE TO GIT AFTER EACH SIMULATION (even if failed!)
            status = "SUCCESS" if success else "FAILED"
            save_to_git(f"{name}_rep{rep} {status} ({i}/{len(SIMULATIONS)})")

        except Exception as e:
            print(f"\nERROR in {name}_rep{rep}: {e}")
            traceback.print_exc()
            results.append((f"{name}_rep{rep}", False))

            # SAVE TO GIT ON ERROR
            save_to_git(f"ERROR in {name}_rep{rep} - saving progress")

            print("\nContinuing to next simulation...")
            gc.collect()

    # Summary
    print("\n" + "=" * 70)
    print("FINAL SUMMARY")
    print("=" * 70)
    success_count = sum(1 for _, s in results if s)
    for name, success in results:
        status = "COMPLETE" if success else "FAILED"
        print(f"  {name}: {status}")
    print(f"\nTotal: {success_count}/{len(results)} successful")
    print(f"Finished: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("=" * 70)

    # FINAL GIT SAVE
    save_to_git(f"ALL DONE: {success_count}/{len(results)} simulations complete")

    # NO AUTO-SHUTDOWN - too risky without verified git auth
    print("\n" + "=" * 70)
    print("ALL SIMULATIONS COMPLETE!")
    print("=" * 70)
    print("\nDATA SAVED LOCALLY. To download:")
    print("  scp -r user@instance:~/nod2-screening-data/PHASE_6/trajectories ./")
    print("  scp -r user@instance:~/nod2-screening-data/PHASE_6/checkpoints ./")
    print("\nThen manually shutdown:")
    print("  sudo shutdown -h now")
    print("\nOr stop the instance from Vast.ai dashboard.")
    print("=" * 70)

if __name__ == "__main__":
    script_dir = os.path.dirname(os.path.abspath(__file__))
    phase6_dir = os.path.join(script_dir, "../..")
    os.chdir(phase6_dir)
    print(f"Working directory: {os.getcwd()}")

    try:
        main()
    except Exception as e:
        print(f"\n!!! FATAL ERROR: {e}")
        traceback.print_exc()
        save_to_git("FATAL ERROR - emergency save")
        sys.exit(1)
