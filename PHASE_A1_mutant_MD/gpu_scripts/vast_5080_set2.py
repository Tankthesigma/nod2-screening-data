#!/usr/bin/env python3
"""
VAST.AI RTX 5080 #2 - Phase A1 Mutant MD Simulations (4 sims, 80ns total)
Sequential execution with auto-shutdown on completion.

SAFETY FEATURES:
- Git push after EACH simulation (crash protection)
- Git push on ANY error (never lose data)
- NaN energy check (fail fast)
- DCD upload to transfer.sh (backup)

Simulations:
1. R702W + Natural rep3 (20ns)
2. G908R + Natural rep1 (20ns)
3. G908R + Natural rep2 (20ns)
4. G908R + Natural rep3 (20ns)

Expected runtime: ~2-3 hours on RTX 5080 with HMR (~600-800 ns/day)
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
    from openff.toolkit.topology import Molecule
    OPENFF_AVAILABLE = True
except ImportError:
    print("FATAL: OpenFF/openmmforcefields not installed!")
    print("Run: conda install -c conda-forge openmm openmmforcefields openff-toolkit -y")
    sys.exit(1)

try:
    from pdbfixer import PDBFixer
    PDBFIXER_AVAILABLE = True
except ImportError:
    print("WARNING: PDBFixer not installed - will skip structure fixing")
    print("Run: conda install -c conda-forge pdbfixer -y")
    PDBFIXER_AVAILABLE = False

# ============================================================
# CONFIGURATION
# ============================================================
GPU_ID = 0

STRUCTURES_DIR = "structures"
LIGANDS_DIR = "ligands"
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

# VAST 5080 SET #2: R702W+Natural rep3 + G908R+Natural rep1,2,3
SIMULATIONS = [
    ("R702W_natural", 3, "NOD2_R702W.pdb", "natural_top_docked.sdf"),
    ("G908R_natural", 1, "NOD2_G908R.pdb", "natural_top_docked.sdf"),
    ("G908R_natural", 2, "NOD2_G908R.pdb", "natural_top_docked.sdf"),
    ("G908R_natural", 3, "NOD2_G908R.pdb", "natural_top_docked.sdf"),
]

# ============================================================
# GIT SAVE FUNCTION (CRITICAL - NEVER LOSE DATA)
# ============================================================
def setup_git_auth():
    """Configure git to use token authentication from environment."""
    token = os.environ.get("GITHUB_TOKEN")
    if not token:
        print(">>> WARNING: GITHUB_TOKEN not set - git push may fail")
        print(">>> Set it with: export GITHUB_TOKEN=ghp_your_token_here")
        return False

    remote_url = f"https://Tankthesigma:{token}@github.com/Tankthesigma/nod2-screening-data.git"
    subprocess.run(["git", "remote", "set-url", "origin", remote_url], check=False, capture_output=True)
    subprocess.run(["git", "config", "user.email", "nod2scout@isef.local"], check=False, capture_output=True)
    subprocess.run(["git", "config", "user.name", "NOD2-Scout"], check=False, capture_output=True)
    print(">>> Git authentication configured")
    return True

def upload_dcd(dcd_path):
    """Upload DCD file to transfer.sh (14 day storage, up to 10GB)."""
    if not os.path.exists(dcd_path):
        print(f">>> DCD not found: {dcd_path}")
        return None

    print(f"\n>>> UPLOADING DCD: {dcd_path}")
    file_size = os.path.getsize(dcd_path) / (1024**3)
    print(f"    Size: {file_size:.2f} GB")

    try:
        result = subprocess.run(
            ["curl", "--upload-file", dcd_path, f"https://transfer.sh/{os.path.basename(dcd_path)}"],
            capture_output=True, text=True, timeout=3600
        )
        if result.returncode == 0:
            url = result.stdout.strip()
            print(f">>> UPLOAD SUCCESS: {url}")
            with open("trajectories/DOWNLOAD_URLS.txt", "a") as f:
                f.write(f"{os.path.basename(dcd_path)}: {url}\n")
            return url
        else:
            print(f">>> UPLOAD FAILED: {result.stderr[:200]}")
            return None
    except Exception as e:
        print(f">>> UPLOAD ERROR: {e}")
        return None

def save_to_git(message="Auto-save simulation results"):
    """Push results to git. Called after each sim and on crash."""
    try:
        print(f"\n>>> SAVING TO GIT: {message}")
        subprocess.run(["git", "add", "-A"], check=False, capture_output=True)
        subprocess.run(
            ["git", "commit", "-m", f"{message} - {datetime.now().strftime('%Y-%m-%d %H:%M')}"],
            check=False, capture_output=True
        )
        # Pull with rebase to avoid conflicts with other vast instance
        subprocess.run(["git", "pull", "--rebase", "origin", "main"], check=False, capture_output=True)
        result = subprocess.run(["git", "push"], check=False, capture_output=True)
        if result.returncode == 0:
            print(">>> GIT PUSH: SUCCESS")
            return True
        else:
            print(f">>> GIT PUSH: FAILED (but local commit saved)")
            print(f"    Error: {result.stderr.decode()[:200]}")
            return False
    except Exception as e:
        print(f">>> GIT SAVE ERROR: {e}")
        return False

# ============================================================
# HELPER FUNCTIONS
# ============================================================
def setup_directories():
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    os.makedirs(LOG_DIR, exist_ok=True)
    os.makedirs(CHECKPOINT_DIR, exist_ok=True)

def fix_pdb_structure(pdb_path):
    """Fix PDB structure using PDBFixer - adds missing atoms/residues."""
    if not PDBFIXER_AVAILABLE:
        print("      PDBFixer not available, loading as-is")
        return PDBFile(pdb_path)

    print("      Running PDBFixer...")
    fixer = PDBFixer(filename=pdb_path)
    fixer.findMissingResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(7.0)

    print(f"      Fixed: {fixer.topology.getNumAtoms()} atoms")

    class FixedPDB:
        def __init__(self, topology, positions):
            self.topology = topology
            self.positions = positions

    return FixedPDB(fixer.topology, fixer.positions)

def check_energy(simulation, stage=""):
    """Check if energy is valid (not NaN/Inf). Returns False if exploded."""
    state = simulation.context.getState(getEnergy=True)
    energy = state.getPotentialEnergy()
    energy_val = energy.value_in_unit(kilojoules_per_mole)

    if math.isnan(energy_val) or math.isinf(energy_val):
        print(f"\n!!! CRITICAL: Energy exploded during {stage}!")
        print(f"!!! Energy = {energy}")
        return False, energy
    return True, energy

def create_system_with_fallback(pdb, ligand_mol, forcefield_kwargs):
    """Create system with OpenFF, falling back to GAFF if needed."""
    modeller = Modeller(pdb.topology, pdb.positions)
    modeller.deleteWater()
    print("      Deleted existing water molecules")

    if ligand_mol is not None:
        print("      Adding ligand to system...")
        ligand_topology = ligand_mol.to_topology().to_openmm()
        ligand_positions = ligand_mol.conformers[0].to_openmm()
        modeller.add(ligand_topology, ligand_positions)
        print(f"      Ligand added: {ligand_mol.n_atoms} atoms")

    if OPENFF_AVAILABLE:
        try:
            print("      Trying OpenFF...")
            system_generator = SystemGenerator(
                forcefields=['amber14-all.xml', 'amber14/tip3pfb.xml'],
                small_molecule_forcefield='openff-2.1.0',
                molecules=[ligand_mol],
                forcefield_kwargs=forcefield_kwargs,
                periodic_forcefield_kwargs={
                    'nonbondedMethod': PME,
                    'nonbondedCutoff': 1.0*nanometers
                }
            )
            modeller.addSolvent(system_generator.forcefield, model='tip3p',
                              padding=1.0*nanometers, ionicStrength=0.15*molar)
            system = system_generator.create_system(modeller.topology)
            for force in system.getForces():
                if isinstance(force, NonbondedForce):
                    force.setUseDispersionCorrection(True)
            print("      SUCCESS: OpenFF (with HMR + dispersion correction)")
            return system, modeller
        except Exception as e:
            print(f"      OpenFF failed: {e}, trying GAFF...")
            modeller = Modeller(pdb.topology, pdb.positions)
            modeller.deleteWater()
            if ligand_mol is not None:
                ligand_topology = ligand_mol.to_topology().to_openmm()
                ligand_positions = ligand_mol.conformers[0].to_openmm()
                modeller.add(ligand_topology, ligand_positions)

    print("      Using GAFF...")
    forcefield = ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')
    gaff_generator = GAFFTemplateGenerator(molecules=[ligand_mol])
    forcefield.registerTemplateGenerator(gaff_generator.generator)
    modeller.addSolvent(forcefield, model='tip3p',
                      padding=1.0*nanometers, ionicStrength=0.15*molar)
    system = forcefield.createSystem(
        modeller.topology, nonbondedMethod=PME, nonbondedCutoff=1.0*nanometers,
        constraints=HBonds, hydrogenMass=HYDROGEN_MASS if USE_HMR else None
    )
    for force in system.getForces():
        if isinstance(force, NonbondedForce):
            force.setUseDispersionCorrection(True)
    print("      SUCCESS: GAFF (with HMR + dispersion correction)")
    return system, modeller

# ============================================================
# SIMULATION FUNCTION
# ============================================================
def run_simulation(name, replicate, pdb_file, sdf_file, sim_num, total_sims):
    """Run a single MD simulation with full error checking."""
    # Seeds: 300000 base for Vast Set2, well separated from other scripts
    random_seed = 300000 + sim_num * 1000 + replicate
    output_prefix = f"{name}_rep{replicate}"

    pdb_path = f"{STRUCTURES_DIR}/{pdb_file}"
    sdf_path = f"{LIGANDS_DIR}/{sdf_file}"

    print("\n" + "=" * 70)
    print(f"[{sim_num}/{total_sims}] {name.upper()} - Replicate {replicate}")
    print(f"Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("=" * 70)

    if not os.path.exists(pdb_path):
        print(f"ERROR: PDB not found: {pdb_path}")
        return False
    if not os.path.exists(sdf_path):
        print(f"ERROR: SDF not found: {sdf_path}")
        return False

    print("\n[1/8] Loading ligand...")
    ligand_mol = Molecule.from_file(sdf_path)
    if isinstance(ligand_mol, list):
        ligand_mol = ligand_mol[0]
    print(f"      Atoms: {ligand_mol.n_atoms}")

    print("\n[2/8] Loading and fixing structure...")
    pdb = fix_pdb_structure(pdb_path)
    modeller_h = Modeller(pdb.topology, pdb.positions)
    print(f"      Ready: {modeller_h.topology.getNumAtoms()} atoms")

    print("\n[3/8] Setting up force fields...")
    forcefield_kwargs = {
        'constraints': HBonds,
        'rigidWater': True,
        'removeCMMotion': True,
        'hydrogenMass': HYDROGEN_MASS if USE_HMR else None
    }
    system, modeller = create_system_with_fallback(modeller_h, ligand_mol, forcefield_kwargs)
    print(f"      Solvated atoms: {modeller.topology.getNumAtoms()}")

    solvated_pdb = f"{OUTPUT_DIR}/{output_prefix}_solvated.pdb"
    with open(solvated_pdb, 'w') as f:
        PDBFile.writeFile(modeller.topology, modeller.positions, f)
    print(f"      Saved: {solvated_pdb}")

    print(f"\n[4/8] Setting up CUDA (GPU {GPU_ID})...")
    integrator = LangevinMiddleIntegrator(TEMPERATURE, FRICTION, TIMESTEP)
    integrator.setRandomNumberSeed(random_seed)
    integrator.setConstraintTolerance(1e-6)
    try:
        platform = Platform.getPlatformByName('CUDA')
        properties = {'Precision': 'mixed', 'DeviceIndex': str(GPU_ID)}
        simulation = Simulation(modeller.topology, system, integrator, platform, properties)
    except Exception as e:
        print(f"      CUDA failed: {e}, trying OpenCL...")
        platform = Platform.getPlatformByName('OpenCL')
        properties = {'Precision': 'mixed', 'OpenCLPlatformIndex': '0', 'DeviceIndex': str(GPU_ID)}
        simulation = Simulation(modeller.topology, system, integrator, platform, properties)
    simulation.context.setPositions(modeller.positions)

    print("\n[5/8] MINIMIZATION")
    ok, energy = check_energy(simulation, "pre-minimization")
    print(f"      Initial: {energy}")

    simulation.minimizeEnergy(maxIterations=MINIMIZATION_STEPS)

    ok, energy = check_energy(simulation, "post-minimization")
    if not ok:
        return False
    print(f"      Final: {energy}")

    print(f"\n[6/8] NVT EQUILIBRATION ({EQUILIBRATION_NVT_PS} ps) - GRADUAL HEATING")

    heating_stages = [
        (50*kelvin, 5),
        (100*kelvin, 10),
        (150*kelvin, 10),
        (200*kelvin, 15),
        (250*kelvin, 20),
        (TEMPERATURE, 40),
    ]

    for temp, duration_ps in heating_stages:
        print(f"      Heating: {duration_ps}ps at {temp}")
        simulation.context.setVelocitiesToTemperature(temp, random_seed)
        integrator.setTemperature(temp)
        steps = int(duration_ps * 1000 / TIMESTEP.value_in_unit(femtoseconds))
        simulation.step(steps)

        ok, energy = check_energy(simulation, f"NVT-{temp}")
        if not ok:
            return False

    print("      NVT complete")

    print(f"\n[7/8] NPT EQUILIBRATION ({EQUILIBRATION_NPT_PS} ps)")
    system.addForce(MonteCarloBarostat(PRESSURE, TEMPERATURE, 100))
    simulation.context.reinitialize(preserveState=True)
    npt_steps = int(EQUILIBRATION_NPT_PS * 1000 / TIMESTEP.value_in_unit(femtoseconds))
    simulation.step(npt_steps)

    ok, energy = check_energy(simulation, "NPT")
    if not ok:
        return False

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

    ok, energy = check_energy(simulation, "production")
    if not ok:
        print("WARNING: Final energy check failed but trajectory saved")

    print(f"\n>>> COMPLETE: {output_prefix}")
    print(f"    Finished: {datetime.now().strftime('%H:%M:%S')}")

    del simulation, system, integrator, modeller
    gc.collect()

    return True

# ============================================================
# MAIN
# ============================================================
def preflight_check():
    """Check ALL input files exist before starting."""
    print("\n[PREFLIGHT] Checking input files...")
    missing = []
    for name, rep, pdb, sdf in SIMULATIONS:
        pdb_path = f"{STRUCTURES_DIR}/{pdb}"
        sdf_path = f"{LIGANDS_DIR}/{sdf}"
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

    print(f"    All input files found.")

    print("\n[PREFLIGHT] Validating ligand files...")
    unique_sdfs = set(sdf for _, _, _, sdf in SIMULATIONS)
    for sdf in unique_sdfs:
        sdf_path = f"{LIGANDS_DIR}/{sdf}"
        try:
            mol = Molecule.from_file(sdf_path)
            if isinstance(mol, list):
                mol = mol[0]
            print(f"    OK {sdf}: {mol.n_atoms} atoms, {mol.n_bonds} bonds")
        except Exception as e:
            print(f"\n!!! FATAL: Cannot load ligand {sdf} !!!")
            print(f"    Error: {e}")
            sys.exit(1)

    print("    All ligands validated successfully.")

def main():
    print("=" * 70)
    print("VAST.AI RTX 5080 #2 - PHASE A1 MUTANT MD SIMULATIONS")
    print(f"Running {len(SIMULATIONS)} simulations ({len(SIMULATIONS) * PRODUCTION_NS}ns total)")
    print(f"Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("=" * 70)

    setup_git_auth()
    preflight_check()
    setup_directories()

    results = []
    for i, (name, rep, pdb, sdf) in enumerate(SIMULATIONS, 1):
        try:
            success = run_simulation(name, rep, pdb, sdf, i, len(SIMULATIONS))
            results.append((f"{name}_rep{rep}", success))

            if success:
                dcd_path = f"{OUTPUT_DIR}/{name}_rep{rep}.dcd"
                upload_dcd(dcd_path)

            status = "SUCCESS" if success else "FAILED"
            save_to_git(f"PHASE_A1 {name}_rep{rep} {status} ({i}/{len(SIMULATIONS)})")

        except Exception as e:
            print(f"\nERROR in {name}_rep{rep}: {e}")
            traceback.print_exc()
            results.append((f"{name}_rep{rep}", False))
            save_to_git(f"PHASE_A1 ERROR in {name}_rep{rep} - saving progress")
            print("\nContinuing to next simulation...")
            gc.collect()

    # Summary
    print("\n" + "=" * 70)
    print("FINAL SUMMARY")
    print("=" * 70)
    success_count = len([s for _, s in results if s])
    for name, success in results:
        status = "COMPLETE" if success else "FAILED"
        print(f"  {name}: {status}")
    print(f"\nTotal: {success_count}/{len(results)} successful")
    print(f"Finished: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("=" * 70)

    push_success = save_to_git(f"PHASE_A1 DONE: {success_count}/{len(results)} simulations complete")

    print("\n" + "=" * 70)
    print("ALL SIMULATIONS COMPLETE!")
    print("=" * 70)

    if push_success:
        print("\n>>> ALL DATA PUSHED TO GITHUB SUCCESSFULLY")
        print(">>> SHUTTING DOWN IN 60 SECONDS...")
        print(">>> (Stop instance from Vast.ai dashboard to cancel)")

        import time
        time.sleep(60)

        print(">>> SHUTTING DOWN NOW...")
        result = subprocess.run(["sudo", "-n", "shutdown", "-h", "now"], capture_output=True)
        if result.returncode != 0:
            os.system('sudo shutdown -h now')
    else:
        print("\n!!! GIT PUSH FAILED - NOT SHUTTING DOWN !!!")
        print("\nData saved locally. Download manually before stopping instance.")

    print("=" * 70)

if __name__ == "__main__":
    script_dir = os.path.dirname(os.path.abspath(__file__))
    phase_a1_dir = os.path.join(script_dir, "..")
    os.chdir(phase_a1_dir)
    print(f"Working directory: {os.getcwd()}")

    try:
        main()
    except Exception as e:
        print(f"\n!!! FATAL ERROR: {e}")
        traceback.print_exc()
        save_to_git("PHASE_A1 FATAL ERROR - emergency save")
        sys.exit(1)
