#!/usr/bin/env python
"""
FEP System Setup with Boresch Restraints (Open Source - No OpenEye Required)

This script prepares solvated systems and alchemical topologies for
Free Energy Perturbation calculations of the R702W mutation effect
on Febuxostat binding to NOD2.

R702W mutation: ARG (charged +1) -> TRP (neutral 0)

Builds solvated systems:
  - wt_complex: WT protein + ligand (with Boresch restraints)
  - mut_complex: R702W mutant protein + ligand (with Boresch restraints)
  - solvent: ligand only in water (no restraints)

Creates alchemical systems with 19-window staged lambda schedule:
  - Stage 1 (windows 0-9): Turn off electrostatics, turn ON restraints
  - Stage 2 (windows 10-18): Turn off sterics, restraints stay ON

Boresch restraints prevent ligand drift at low lambda and enable
quantitative absolute binding free energy calculations. The standard
state correction must be applied after MBAR analysis.

NOTE: This is SETUP ONLY. After running this script:
  1. Run equilibration MD for each system
  2. Run alchemical FEP sampling at each lambda window
  3. Collect u_nk energies and analyze with MBAR
  4. Apply standard state correction to final ΔG_bind

Usage:
  python setup_fep_pmx.py
"""
import os
import sys
import glob
import shutil
import numpy as np

# Boresch restraint force constants (kcal/mol units)
# These must be consistent between restraint creation and standard state correction
K_DISTANCE = 10.0   # kcal/mol/Å²
K_ANGLE = 10.0      # kcal/mol/rad²
K_DIHEDRAL = 10.0   # kcal/mol/rad²
TEMPERATURE_K = 310.0  # Kelvin - must match production temperature


def normalize_resid(val):
    """Normalize residue ID to int for comparison."""
    if val is None:
        return None
    try:
        return int(str(val).rstrip('ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz'))
    except (ValueError, TypeError):
        return None


def get_ligand_name_from_sdf(sdf_path):
    """Extract ligand name from SDF file."""
    from rdkit import Chem
    suppl = Chem.SDMolSupplier(sdf_path, removeHs=False)
    for mol in suppl:
        if mol is not None:
            name = mol.GetProp("_Name") if mol.HasProp("_Name") else None
            if not name:
                name = os.path.splitext(os.path.basename(sdf_path))[0]
            return name, mol
    return None, None


def select_boresch_anchors(topology, positions, ligand_atom_indices, cutoff_nm=1.0):
    """
    Select anchor atoms for Boresch restraints.

    Returns:
        dict with keys:
            'protein': [P1, P2, P3] - 3 protein CA atom indices
            'ligand': [L1, L2, L3] - 3 ligand heavy atom indices
            'values': dict of equilibrium values for restraints
    """
    from openmm import unit

    # Get positions as numpy array in nm
    pos_unit = positions.value_in_unit(unit.nanometer)
    try:
        pos_nm = np.array([[p.x, p.y, p.z] for p in pos_unit])
    except (AttributeError, TypeError):
        pos_nm = np.asarray(pos_unit)

    # Calculate ligand centroid
    ligand_coords = pos_nm[ligand_atom_indices]
    ligand_centroid = np.mean(ligand_coords, axis=0)

    # Find protein CA atoms within cutoff of ligand centroid
    protein_ca_indices = []
    for atom in topology.atoms():
        if atom.name == 'CA' and atom.residue.name not in ['HOH', 'WAT', 'NA', 'CL', 'Na+', 'Cl-']:
            if atom.index not in ligand_atom_indices:
                dist = np.linalg.norm(pos_nm[atom.index] - ligand_centroid)
                if dist < cutoff_nm:
                    protein_ca_indices.append((atom.index, dist, atom.residue.name, atom.residue.id))

    if len(protein_ca_indices) < 3:
        return None, "Not enough protein CA atoms near ligand"

    # Sort by distance and select 3 well-separated CA atoms
    protein_ca_indices.sort(key=lambda x: x[1])

    # Select P1 (closest), then P2 and P3 to form a triangle
    P1 = protein_ca_indices[0][0]
    P1_pos = pos_nm[P1]

    # Find P2: next closest that's at least 0.5 nm from P1
    P2 = None
    for idx, dist, resname, resid in protein_ca_indices[1:]:
        if np.linalg.norm(pos_nm[idx] - P1_pos) > 0.5:
            P2 = idx
            break
    if P2 is None:
        P2 = protein_ca_indices[1][0]
    P2_pos = pos_nm[P2]

    # Find P3: next that's at least 0.5 nm from both P1 and P2
    P3 = None
    for idx, dist, resname, resid in protein_ca_indices[2:]:
        if (np.linalg.norm(pos_nm[idx] - P1_pos) > 0.5 and
            np.linalg.norm(pos_nm[idx] - P2_pos) > 0.5):
            P3 = idx
            break
    if P3 is None:
        P3 = protein_ca_indices[2][0]

    # Select ligand heavy atoms (non-hydrogen, non-terminal)
    # Get heavy atoms only
    ligand_heavy = []
    for idx in ligand_atom_indices:
        atom = list(topology.atoms())[idx]
        if atom.element is None or atom.element.symbol != 'H':
            # Calculate connectivity (number of bonds to other ligand atoms)
            ligand_heavy.append(idx)

    if len(ligand_heavy) < 3:
        return None, "Not enough ligand heavy atoms"

    # Select L1 (closest to P1), L2, L3 to form a triangle
    ligand_heavy_dist = [(idx, np.linalg.norm(pos_nm[idx] - P1_pos)) for idx in ligand_heavy]
    ligand_heavy_dist.sort(key=lambda x: x[1])

    L1 = ligand_heavy_dist[0][0]
    L1_pos = pos_nm[L1]

    # Find L2: at least 0.2 nm from L1
    L2 = None
    for idx, dist in ligand_heavy_dist[1:]:
        if np.linalg.norm(pos_nm[idx] - L1_pos) > 0.2:
            L2 = idx
            break
    if L2 is None:
        L2 = ligand_heavy_dist[1][0]
    L2_pos = pos_nm[L2]

    # Find L3: at least 0.2 nm from both L1 and L2
    L3 = None
    for idx, dist in ligand_heavy_dist[2:]:
        if (np.linalg.norm(pos_nm[idx] - L1_pos) > 0.2 and
            np.linalg.norm(pos_nm[idx] - L2_pos) > 0.2):
            L3 = idx
            break
    if L3 is None:
        L3 = ligand_heavy_dist[2][0]

    # Calculate equilibrium values
    def calc_distance(i, j):
        return np.linalg.norm(pos_nm[i] - pos_nm[j])

    def calc_angle(i, j, k):
        v1 = pos_nm[i] - pos_nm[j]
        v2 = pos_nm[k] - pos_nm[j]
        cos_angle = np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2))
        return np.arccos(np.clip(cos_angle, -1, 1))

    def calc_dihedral(i, j, k, l):
        b1 = pos_nm[j] - pos_nm[i]
        b2 = pos_nm[k] - pos_nm[j]
        b3 = pos_nm[l] - pos_nm[k]
        n1 = np.cross(b1, b2)
        n2 = np.cross(b2, b3)
        m1 = np.cross(n1, b2 / np.linalg.norm(b2))
        x = np.dot(n1, n2)
        y = np.dot(m1, n2)
        return np.arctan2(y, x)

    values = {
        'r_P1_L1': calc_distance(P1, L1),
        'theta_P2_P1_L1': calc_angle(P2, P1, L1),
        'theta_P1_L1_L2': calc_angle(P1, L1, L2),
        'phi_P3_P2_P1_L1': calc_dihedral(P3, P2, P1, L1),
        'phi_P2_P1_L1_L2': calc_dihedral(P2, P1, L1, L2),
        'phi_P1_L1_L2_L3': calc_dihedral(P1, L1, L2, L3),
    }

    return {
        'protein': [P3, P2, P1],  # Order: P3-P2-P1-L1-L2-L3
        'ligand': [L1, L2, L3],
        'values': values,
    }, None


def create_boresch_restraints(anchors, k_distance=10.0, k_angle=10.0, k_dihedral=10.0):
    """
    Create Boresch restraint forces with lambda scaling.

    Force constants in kcal/mol units:
        k_distance: kcal/mol/Å²
        k_angle: kcal/mol/rad²
        k_dihedral: kcal/mol/rad²

    Returns list of (force, particle_indices) tuples.
    """
    from openmm import CustomBondForce, CustomAngleForce, CustomTorsionForce

    P3, P2, P1 = anchors['protein']
    L1, L2, L3 = anchors['ligand']
    values = anchors['values']

    # Convert force constants from kcal/mol to kJ/mol
    # 1 kcal = 4.184 kJ
    k_dist_kjmol = k_distance * 4.184 * 100  # kcal/mol/Å² -> kJ/mol/nm²
    k_ang_kjmol = k_angle * 4.184            # kcal/mol/rad² -> kJ/mol/rad²
    k_dih_kjmol = k_dihedral * 4.184         # kcal/mol/rad² -> kJ/mol/rad²

    restraints = []

    # 1. Distance restraint P1-L1 (standard harmonic: 0.5*k*(r-r0)^2)
    # lambda_restraints: 0 = restraints OFF, 1 = restraints ON
    dist_force = CustomBondForce('lambda_restraints * 0.5 * k * (r - r0)^2')
    dist_force.addGlobalParameter('lambda_restraints', 0.0)
    dist_force.addPerBondParameter('k')
    dist_force.addPerBondParameter('r0')
    dist_force.addBond(P1, L1, [k_dist_kjmol, values['r_P1_L1']])
    dist_force.setForceGroup(10)
    restraints.append(('distance_P1_L1', dist_force))

    # 2. Angle restraint P2-P1-L1 (standard harmonic: 0.5*k*(θ-θ0)^2)
    angle1_force = CustomAngleForce('lambda_restraints * 0.5 * k * (theta - theta0)^2')
    angle1_force.addGlobalParameter('lambda_restraints', 0.0)
    angle1_force.addPerAngleParameter('k')
    angle1_force.addPerAngleParameter('theta0')
    angle1_force.addAngle(P2, P1, L1, [k_ang_kjmol, values['theta_P2_P1_L1']])
    angle1_force.setForceGroup(10)
    restraints.append(('angle_P2_P1_L1', angle1_force))

    # 3. Angle restraint P1-L1-L2 (standard harmonic: 0.5*k*(θ-θ0)^2)
    angle2_force = CustomAngleForce('lambda_restraints * 0.5 * k * (theta - theta0)^2')
    angle2_force.addGlobalParameter('lambda_restraints', 0.0)
    angle2_force.addPerAngleParameter('k')
    angle2_force.addPerAngleParameter('theta0')
    angle2_force.addAngle(P1, L1, L2, [k_ang_kjmol, values['theta_P1_L1_L2']])
    angle2_force.setForceGroup(10)
    restraints.append(('angle_P1_L1_L2', angle2_force))

    # 4. Dihedral restraint P3-P2-P1-L1 (periodic-aware harmonic: 0.5*k*min(Δφ, 2π-Δφ)^2)
    dih1_force = CustomTorsionForce('lambda_restraints * 0.5 * k * min(dtheta, 2*pi - dtheta)^2; dtheta = abs(theta - theta0); pi = 3.14159265359')
    dih1_force.addGlobalParameter('lambda_restraints', 0.0)
    dih1_force.addPerTorsionParameter('k')
    dih1_force.addPerTorsionParameter('theta0')
    dih1_force.addTorsion(P3, P2, P1, L1, [k_dih_kjmol, values['phi_P3_P2_P1_L1']])
    dih1_force.setForceGroup(10)
    restraints.append(('dihedral_P3_P2_P1_L1', dih1_force))

    # 5. Dihedral restraint P2-P1-L1-L2
    dih2_force = CustomTorsionForce('lambda_restraints * 0.5 * k * min(dtheta, 2*pi - dtheta)^2; dtheta = abs(theta - theta0); pi = 3.14159265359')
    dih2_force.addGlobalParameter('lambda_restraints', 0.0)
    dih2_force.addPerTorsionParameter('k')
    dih2_force.addPerTorsionParameter('theta0')
    dih2_force.addTorsion(P2, P1, L1, L2, [k_dih_kjmol, values['phi_P2_P1_L1_L2']])
    dih2_force.setForceGroup(10)
    restraints.append(('dihedral_P2_P1_L1_L2', dih2_force))

    # 6. Dihedral restraint P1-L1-L2-L3
    dih3_force = CustomTorsionForce('lambda_restraints * 0.5 * k * min(dtheta, 2*pi - dtheta)^2; dtheta = abs(theta - theta0); pi = 3.14159265359')
    dih3_force.addGlobalParameter('lambda_restraints', 0.0)
    dih3_force.addPerTorsionParameter('k')
    dih3_force.addPerTorsionParameter('theta0')
    dih3_force.addTorsion(P1, L1, L2, L3, [k_dih_kjmol, values['phi_P1_L1_L2_L3']])
    dih3_force.setForceGroup(10)
    restraints.append(('dihedral_P1_L1_L2_L3', dih3_force))

    return restraints


def calc_boresch_correction(anchors, temperature_K=310.0,
                            k_r=10.0, k_theta=10.0, k_phi=10.0):
    """
    Calculate the standard state correction for Boresch restraints.

    For harmonic restraints with U = 0.5*k*Δx², the partition function gives
    a Gaussian width σ = sqrt(kT/k).

    The Boresch correction accounts for the difference between:
    - The restrained ligand's configurational volume in the binding site
    - The standard state volume V° = 1660 Å³ (1 M concentration)

    ΔG_corr = kT * ln[ (8π² V°) / (r₀² sin(θ_A) sin(θ_B)) *
                       σ_r * σ_θA * σ_θB * σ_φA * σ_φB * σ_φC / (2π)³ ]

    For ABFE (decoupling):
      ΔG_bind = ΔG_decouple_solvent - ΔG_decouple_complex + ΔG_corr

    Args:
        anchors: dict with 'values' containing equilibrium restraint values
        temperature_K: temperature in Kelvin (default 300 K)
        k_r: distance force constant in kcal/mol/Å²
        k_theta: angle force constant in kcal/mol/rad²
        k_phi: dihedral force constant in kcal/mol/rad²

    Returns:
        ΔG_corr in kcal/mol (typically positive, ~2-6 kcal/mol)
    """
    import math

    values = anchors['values']

    # Standard state volume (1 M concentration)
    V_standard = 1660.0  # Å³

    # kT at temperature
    kB = 0.001987204  # kcal/mol/K
    kT = kB * temperature_K

    # Equilibrium values
    r0 = values['r_P1_L1'] * 10  # nm -> Å
    theta_A = values['theta_P2_P1_L1']
    theta_B = values['theta_P1_L1_L2']

    # Check for ill-conditioned angles (close to 0 or π)
    if math.sin(theta_A) < 0.1 or math.sin(theta_B) < 0.1:
        print(f"  [WARNING] Angle close to 0 or π - correction may be inaccurate!")
        print(f"    sin(θ_A) = {math.sin(theta_A):.3f}, sin(θ_B) = {math.sin(theta_B):.3f}")

    # Gaussian widths from partition function
    # For U = 0.5*k*Δx², the Gaussian integral gives: ∫exp(-U/kT)dx = sqrt(2π*kT/k)
    # So the effective width is sqrt(2π) * sqrt(kT/k)
    sigma_r = math.sqrt(2 * math.pi * kT / k_r)           # Å
    sigma_theta = math.sqrt(2 * math.pi * kT / k_theta)   # rad
    sigma_phi = math.sqrt(2 * math.pi * kT / k_phi)       # rad

    # Boresch formula: restrained configurational volume
    # V_site = r₀² * sin(θA) * sin(θB) * σ_r * σ_θA * σ_θB * σ_φA * σ_φB * σ_φC / (8π²)
    # The 8π² normalizes for the orientational degrees of freedom
    V_site = (r0**2 * math.sin(theta_A) * math.sin(theta_B) *
              sigma_r * sigma_theta**2 * sigma_phi**3) / (8 * math.pi**2)

    # Correction: free energy to release restraints to standard state
    # ΔG_corr = -kT * ln(V_site / V_standard) = kT * ln(V_standard / V_site)
    dG_corr = kT * math.log(V_standard / V_site)

    return dG_corr


def calc_pme_charge_correction(ligand_charge, box_length_nm, temperature_K=310.0):
    """
    Calculate PME finite-size correction for charged ligand decoupling.

    When decoupling a charged ligand in periodic PME, there's a finite-size
    artifact from the neutralizing background charge. This correction follows
    Rocklin et al. (J Chem Phys 2013) and is typically small (~0.1-1 kcal/mol).

    For ABFE with charged ligand:
      ΔG_bind = ΔG_solvent - ΔG_complex + ΔG_Boresch + ΔG_PME_correction

    Args:
        ligand_charge: net charge of ligand (in elementary charge units)
        box_length_nm: cubic box edge length in nm
        temperature_K: temperature in Kelvin

    Returns:
        ΔG_PME_correction in kcal/mol (0 if ligand is neutral)
    """
    import math

    if abs(ligand_charge) < 0.01:
        return 0.0  # No correction needed for neutral ligand

    # Physical constants
    kB = 0.001987204  # kcal/mol/K
    kT = kB * temperature_K

    # Convert box to Angstroms
    L = box_length_nm * 10  # nm -> Å

    # Madelung constant for cubic lattice (self-energy term)
    # ξ_EW ≈ -2.837297 for Ewald sum
    xi_ewald = -2.837297

    # Coulomb constant: 332.0637 kcal·Å/(mol·e²)
    coulomb_const = 332.0637

    # Self-interaction correction
    # ΔG_self = ξ_EW * q² * k_c / (2 * L)
    dG_self = xi_ewald * (ligand_charge ** 2) * coulomb_const / (2 * L)

    # Net correction (typically negative for removing charged ligand)
    # For decoupling: we're removing the ligand, so correction is opposite sign
    dG_correction = -dG_self

    return dG_correction


def main():
    # Safe stdout configuration
    try:
        sys.stdout.reconfigure(line_buffering=True)
        sys.stderr.reconfigure(line_buffering=True)
    except (AttributeError, OSError):
        pass
    os.environ["PYTHONUNBUFFERED"] = "1"

    print("="*60)
    print("FEP SETUP WITH PMX (Open Source)")
    print("="*60)
    print()
    print("Working directory:", os.getcwd())
    print()

    # ============================================================
    # SECTION 1: CHECK DEPENDENCIES
    # ============================================================
    print("[1] Checking dependencies...")

    try:
        import pmx
        pmx_version = getattr(pmx, '__version__', 'unknown')
        print(f"  [PASS] pmx version: {pmx_version}")
    except ImportError as e:
        print(f"  [WARNING] pmx not installed: {e}")
        print("  (Optional - mutation done via mutate_r702w.py)")

    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem, rdMolAlign, rdFMCS
        print("  [PASS] RDKit")
    except ImportError:
        print("  [FAIL] RDKit - required for ligand handling")
        sys.exit(1)

    try:
        import mdtraj as md
        print("  [PASS] MDTraj")
    except ImportError:
        print("  [FAIL] MDTraj - required for trajectory handling")
        sys.exit(1)

    try:
        from openmm import app, unit, XmlSerializer
        from openmm.app import PDBFile, Modeller
        print("  [PASS] OpenMM")
    except ImportError:
        print("  [FAIL] OpenMM - required for simulation setup")
        sys.exit(1)

    try:
        from openff.toolkit import Molecule as OFFMolecule
        from openmmforcefields.generators import SystemGenerator
        print("  [PASS] OpenFF/OpenMMForceFields")
    except ImportError as e:
        print(f"  [FAIL] OpenFF: {e}")
        print("  OpenFF is REQUIRED for ligand parameterization")
        sys.exit(1)

    try:
        from openmmtools.alchemy import AbsoluteAlchemicalFactory, AlchemicalRegion
        print("  [PASS] OpenMMTools (alchemy)")
        HAS_ALCHEMY = True
    except ImportError as e:
        print(f"  [WARNING] OpenMMTools alchemy: {e}")
        print("  Install with: pip install openmmtools")
        HAS_ALCHEMY = False

    print()

    # ============================================================
    # SECTION 2: FIND INPUT FILES
    # ============================================================
    print("="*60)
    print("[2] FINDING INPUT FILES")
    print("="*60)

    dcd_patterns = ["**/WT*[Ff]ebuxostat*rep2*.dcd", "**/WT*[Ff]ebuxostat*.dcd",
                    "**/febuxostat*rep2*.dcd", "PHASE_6/**/*.dcd"]
    pdb_patterns = ["**/WT*[Ff]ebuxostat*rep2*solvated*.pdb", "**/WT*[Ff]ebuxostat*rep2*.pdb",
                    "**/WT*[Ff]ebuxostat*solvated*.pdb", "PHASE_6/**/*solvated*.pdb"]
    sdf_patterns = ["**/febuxostat.sdf", "**/[Ff]ebuxostat*.sdf", "**/ligand*.sdf"]

    def find_files(patterns):
        found = []
        for p in patterns:
            found.extend(glob.glob(p, recursive=True))
        return sorted(set(found))

    dcd_files = find_files(dcd_patterns)
    pdb_files = find_files(pdb_patterns)
    sdf_files = find_files(sdf_patterns)

    print(f"Found: {len(dcd_files)} DCD, {len(pdb_files)} PDB, {len(sdf_files)} SDF")

    if not dcd_files:
        print("[FAIL] No DCD files found!")
        sys.exit(1)

    # Select WT febuxostat trajectory
    def is_wt_febuxostat(f):
        fl = f.lower()
        return "febuxostat" in fl and "r702w" not in fl and "g908r" not in fl and "mutant" not in fl

    wt_febuxostat = [f for f in dcd_files if is_wt_febuxostat(f)]
    if not wt_febuxostat:
        print("[FAIL] No WT Febuxostat trajectories found!")
        sys.exit(1)

    rep2_feb = [f for f in wt_febuxostat if "rep2" in f.lower()]
    dcd_path = max(rep2_feb if rep2_feb else wt_febuxostat, key=os.path.getsize)
    print(f"[SELECTED] DCD: {os.path.basename(dcd_path)}")

    # Find matching PDB
    viable_pdbs = []
    for pdb in pdb_files:
        try:
            test_frame = md.load_frame(dcd_path, 0, top=pdb)
            if test_frame.n_atoms > 1000:
                viable_pdbs.append({"path": pdb, "n_atoms": test_frame.n_atoms})
                print(f"  [OK] {os.path.basename(pdb)}: {test_frame.n_atoms} atoms")
        except Exception as e:
            print(f"  [SKIP] {os.path.basename(pdb)}: {str(e)[:50]}")

    if not viable_pdbs:
        print("[FAIL] No PDB matches DCD topology!")
        sys.exit(1)

    def pdb_score(p):
        name = p["path"].lower()
        wt_bonus = 100 if "wt" in name and "febuxostat" in name else 0
        return wt_bonus + p["n_atoms"]

    viable_pdbs.sort(key=pdb_score, reverse=True)
    pdb_path = viable_pdbs[0]["path"]
    print(f"[SELECTED] PDB: {os.path.basename(pdb_path)} ({viable_pdbs[0]['n_atoms']} atoms)")

    # Select SDF and get ligand name explicitly
    if not sdf_files:
        print("[FAIL] No SDF files found!")
        sys.exit(1)

    febuxostat_sdfs = [s for s in sdf_files if "febuxostat" in s.lower()]
    sdf_path = febuxostat_sdfs[0] if febuxostat_sdfs else sdf_files[0]
    print(f"[SELECTED] SDF: {os.path.basename(sdf_path)}")

    # Get ligand info from SDF explicitly
    ligand_name, sdf_mol = get_ligand_name_from_sdf(sdf_path)
    if sdf_mol is None:
        print("[FAIL] Could not read ligand from SDF!")
        sys.exit(1)

    sdf_heavy = sdf_mol.GetNumHeavyAtoms()
    sdf_total = sdf_mol.GetNumAtoms()
    sdf_charge = Chem.GetFormalCharge(sdf_mol)
    print(f"[LIGAND] Name: {ligand_name}")
    print(f"[LIGAND] Atoms: {sdf_total} total, {sdf_heavy} heavy, charge {sdf_charge}")

    print()

    # ============================================================
    # SECTION 3: PREPARE STRUCTURES
    # ============================================================
    print("="*60)
    print("[3] PREPARING STRUCTURES")
    print("="*60)

    # Create output directories
    for d in ["fep_pmx", "fep_pmx/wt_complex", "fep_pmx/mut_complex",
              "fep_pmx/solvent"]:
        os.makedirs(d, exist_ok=True)

    # Load last frame
    print("[INFO] Loading last frame from trajectory...")
    with md.open(dcd_path) as traj_file:
        n_frames = len(traj_file)
    print(f"  Total frames: {n_frames}")
    last_frame = md.load_frame(dcd_path, n_frames - 1, top=pdb_path)
    print(f"  Last frame atoms: {last_frame.n_atoms}")

    try:
        last_frame.image_molecules(inplace=True)
        print("  [PASS] Molecules imaged for PBC")
    except Exception as e:
        print(f"  [WARNING] Imaging failed: {e}")

    # Find ligand residue by matching heavy atom count with SDF
    print(f"[INFO] Finding ligand residue matching SDF ({sdf_heavy} heavy atoms)...")

    EXCLUDE = {"ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE","LEU","LYS","MET","PHE",
               "PRO","SER","THR","TRP","TYR","VAL","HIE","HID","HIP","HSD","HSE","HSP","CYX","CYM",
               "ACE","NME","NH2","HOH","WAT","TIP3","SOL","NA","CL","K","MG","CA","ZN","NA+","CL-"}

    def count_heavy(residue):
        return sum(1 for a in residue.atoms if a.element is None or a.element.symbol != "H")

    # Find residue matching the SDF heavy atom count
    ligand_res = None
    ligand_res_idx = None
    ligand_heavy = None

    for r in last_frame.topology.residues:
        if r.name.upper() not in EXCLUDE:
            heavy = count_heavy(r)
            # Match by heavy atom count (allow +-1 for protonation differences)
            if abs(heavy - sdf_heavy) <= 1:
                ligand_res = r
                ligand_res_idx = r.index
                ligand_heavy = heavy
                break

    if ligand_res is None:
        print(f"[FAIL] No residue found matching ligand heavy atom count ({sdf_heavy})!")
        print("  Non-standard residues found:")
        for r in last_frame.topology.residues:
            if r.name.upper() not in EXCLUDE:
                print(f"    {r.name}: {count_heavy(r)} heavy atoms")
        sys.exit(1)

    print(f"[SELECTED] Ligand residue: {ligand_res.name} (index {ligand_res_idx}, {ligand_heavy} heavy atoms)")

    # Find ARG 702
    PROTEIN_RES = {"ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE",
                   "LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL",
                   "HIE","HID","HIP","HSD","HSE","HSP","CYX","CYM","ARN","AR0"}

    res702_matches = []
    for chain in last_frame.topology.chains:
        for r in chain.residues:
            if r.resSeq == 702 and r.name.upper() in PROTEIN_RES:
                res702_matches.append((chain, r))

    if not res702_matches:
        print("[FAIL] No protein residue 702 found!")
        sys.exit(1)

    if len(res702_matches) > 1:
        print(f"[WARNING] Multiple residue 702 found in {len(res702_matches)} chains")
        ligand_atoms = [a for a in last_frame.topology.atoms if a.residue.index == ligand_res_idx]
        ligand_coords = last_frame.xyz[0, [a.index for a in ligand_atoms], :] * 10
        ligand_centroid = np.mean(ligand_coords, axis=0)

        best_dist = float('inf')
        best_match = None
        for chain, res in res702_matches:
            res_atoms = list(res.atoms)
            res_coords = last_frame.xyz[0, [a.index for a in res_atoms], :] * 10
            dist = np.linalg.norm(np.mean(res_coords, axis=0) - ligand_centroid)
            if dist < best_dist:
                best_dist = dist
                best_match = (chain, res)
        target_chain, target_res = best_match
    else:
        target_chain, target_res = res702_matches[0]

    if target_res.name.upper() not in ["ARG", "ARN", "AR0"]:
        print(f"[FAIL] Residue 702 is {target_res.name}, expected ARG!")
        sys.exit(1)

    chain_id = target_chain.chain_id if target_chain.chain_id else "A"
    print(f"[PASS] Found ARG 702 in chain {chain_id}")
    print(f"[INFO] Mutation: R702W (ARG +1 -> TRP 0, charge-changing)")

    # Extract protein and ligand
    protein_idx = list(last_frame.topology.select("protein"))
    ligand_idx = [a.index for a in last_frame.topology.atoms if a.residue.index == ligand_res_idx]
    ligand_idx_set = set(ligand_idx)
    protein_idx = [i for i in protein_idx if i not in ligand_idx_set]

    print(f"[INFO] Protein atoms: {len(protein_idx)}, Ligand atoms: {len(ligand_idx)}")

    # Save structures
    protein_only = last_frame.atom_slice(protein_idx)
    protein_only.save("fep_pmx/wt_complex/protein.pdb")
    print("[PASS] Saved fep_pmx/wt_complex/protein.pdb")

    ligand_only = last_frame.atom_slice(ligand_idx)
    ligand_only.save("fep_pmx/wt_complex/ligand_md.pdb")
    print("[PASS] Saved fep_pmx/wt_complex/ligand_md.pdb")

    # Align SDF to MD pose using MCS
    print()
    print("[INFO] Aligning SDF to MD pose...")

    if sdf_mol.GetNumConformers() == 0:
        embed_result = AllChem.EmbedMolecule(sdf_mol, AllChem.ETKDG())
        if embed_result == -1:
            print("[FAIL] Could not generate 3D conformer for ligand!")
            sys.exit(1)

    md_ligand_mol = Chem.MolFromPDBFile("fep_pmx/wt_complex/ligand_md.pdb", removeHs=False)
    if md_ligand_mol is None or md_ligand_mol.GetNumConformers() == 0:
        print("[FAIL] Could not load MD ligand for alignment!")
        sys.exit(1)

    # MCS-based alignment (required, no fallback)
    mcs = rdFMCS.FindMCS([sdf_mol, md_ligand_mol],
                          atomCompare=rdFMCS.AtomCompare.CompareAny,
                          bondCompare=rdFMCS.BondCompare.CompareAny,
                          ringMatchesRingOnly=False,
                          completeRingsOnly=False,
                          timeout=10)

    if mcs.numAtoms < sdf_heavy * 0.5:
        print(f"[FAIL] MCS alignment failed - only {mcs.numAtoms} atoms matched!")
        print("  SDF and MD ligand may not be the same molecule")
        sys.exit(1)

    mcs_query = Chem.MolFromSmarts(mcs.smartsString)
    sdf_match = sdf_mol.GetSubstructMatch(mcs_query)
    md_match = md_ligand_mol.GetSubstructMatch(mcs_query)

    if not sdf_match or not md_match:
        print("[FAIL] Could not get substructure matches for alignment!")
        sys.exit(1)

    atom_map = list(zip(sdf_match, md_match))
    rmsd = rdMolAlign.AlignMol(sdf_mol, md_ligand_mol, atomMap=atom_map)
    print(f"[PASS] Aligned SDF to MD pose (RMSD: {rmsd:.2f} A, {len(atom_map)} atoms)")

    writer = Chem.SDWriter("fep_pmx/wt_complex/ligand.sdf")
    writer.write(sdf_mol)
    writer.close()
    print("[PASS] Saved fep_pmx/wt_complex/ligand.sdf")

    print()

    # ============================================================
    # SECTION 4: PMX MUTATION
    # ============================================================
    print("="*60)
    print("[4] PMX MUTATION SETUP (R702W)")
    print("="*60)

    mutation_successful = False
    target_chain_for_check = chain_id if chain_id else ""

    # Check existing mutant
    if os.path.exists("fep_pmx/mut_complex/protein.pdb"):
        print("[INFO] Checking existing mutant file...")
        try:
            existing_mut = PDBFile("fep_pmx/mut_complex/protein.pdb")
            trp702_in_target = False
            arg702_in_target = False

            for res in existing_mut.topology.residues():
                res_id_norm = normalize_resid(res.id)
                res_chain = res.chain.id if hasattr(res, 'chain') and hasattr(res.chain, 'id') else ""
                is_target_chain = (target_chain_for_check == "" or res_chain == "" or
                                   res_chain == target_chain_for_check)

                if res_id_norm == 702 and is_target_chain:
                    if res.name == 'TRP':
                        trp702_in_target = True
                    elif res.name in ['ARG', 'ARN', 'AR0']:
                        arg702_in_target = True

            if trp702_in_target and not arg702_in_target:
                print("  [PASS] Existing mutant has TRP 702 - using it")
                mutation_successful = True
        except Exception as e:
            print(f"  [INFO] Could not check existing mutant: {e}")

    if not mutation_successful:
        # Copy WT protein and instruct user to mutate
        shutil.copy("fep_pmx/wt_complex/protein.pdb", "fep_pmx/mut_complex/protein.pdb")
        print("[WARNING] Automatic mutation not available!")
        print("  Run: python mutate_r702w.py")
        print(f"  Or mutate ARG 702 -> TRP in chain {chain_id} using PyMOL")

    # Copy ligand to mutant complex
    shutil.copy("fep_pmx/wt_complex/ligand.sdf", "fep_pmx/mut_complex/ligand.sdf")

    print()

    # ============================================================
    # SECTION 5: CREATE SOLVATED COMPLEXES
    # ============================================================
    print("="*60)
    print("[5] CREATING SOLVATED SYSTEMS")
    print("="*60)

    # Load ligand with OpenFF
    off_mol = OFFMolecule.from_file("fep_pmx/wt_complex/ligand.sdf", allow_undefined_stereo=True)
    print(f"[PASS] OpenFF molecule: {off_mol.n_atoms} atoms, charge {off_mol.total_charge}")

    # Create system generator with PME settings for periodic systems
    # useDispersionCorrection=True adds analytical LJ tail correction for accuracy
    system_gen = SystemGenerator(
        forcefields=['amber14/protein.ff14SB.xml', 'amber14/tip3p.xml'],
        small_molecule_forcefield='openff-2.1.0',
        molecules=[off_mol],
        periodic_forcefield_kwargs={
            'nonbondedMethod': app.PME,
            'nonbondedCutoff': 1.0 * unit.nanometer,
            'constraints': app.HBonds,
            'useDispersionCorrection': True  # LJ long-range correction
        }
    )
    print("[PASS] SystemGenerator with OpenFF 2.1.0, PME, dispersion correction")

    # Solvation parameters
    PADDING = 1.0 * unit.nanometer
    IONIC_STRENGTH = 150 * unit.millimolar

    systems_to_build = [
        ("wt_complex", "fep_pmx/wt_complex/protein.pdb", True),
        ("mut_complex", "fep_pmx/mut_complex/protein.pdb", True),
        ("solvent", None, False),  # Ligand only in water
    ]

    built_systems = {}

    for sys_name, protein_path, has_protein in systems_to_build:
        print()
        print(f"[INFO] Building {sys_name}...")

        out_dir = f"fep_pmx/{sys_name}"
        os.makedirs(out_dir, exist_ok=True)

        # Load ligand topology and positions
        ligand_top = off_mol.to_topology().to_openmm()
        ligand_pos = off_mol.conformers[0].to_openmm()

        if has_protein:
            # Load protein
            protein_pdb = PDBFile(protein_path)
            modeller = Modeller(protein_pdb.topology, protein_pdb.positions)

            # Add ligand to protein
            modeller.add(ligand_top, ligand_pos)
            print(f"  Protein+Ligand: {modeller.topology.getNumAtoms()} atoms")
        else:
            # Ligand only
            modeller = Modeller(ligand_top, ligand_pos)
            print(f"  Ligand only: {modeller.topology.getNumAtoms()} atoms")

        # Add solvent (TIP3P with ions)
        modeller.addSolvent(
            system_gen.forcefield,
            model='tip3p',
            padding=PADDING,
            ionicStrength=IONIC_STRENGTH,
            positiveIon='Na+',
            negativeIon='Cl-'
        )
        print(f"  Solvated: {modeller.topology.getNumAtoms()} atoms")

        # Count ions and waters
        n_na = sum(1 for r in modeller.topology.residues() if r.name in ['NA', 'Na+'])
        n_cl = sum(1 for r in modeller.topology.residues() if r.name in ['CL', 'Cl-'])
        n_water = sum(1 for r in modeller.topology.residues() if r.name in ['HOH', 'WAT'])
        print(f"  Ions: {n_na} Na+, {n_cl} Cl-")
        print(f"  Waters: {n_water} molecules ({n_water * 3} atoms)")

        # Create system (PME settings from SystemGenerator constructor)
        system = system_gen.create_system(modeller.topology)
        print(f"  [PASS] System created: {system.getNumParticles()} particles")

        # Add MonteCarloBarostat for NPT equilibration
        from openmm import MonteCarloBarostat
        barostat = MonteCarloBarostat(1*unit.bar, TEMPERATURE_K*unit.kelvin, 25)
        system.addForce(barostat)
        print(f"  [PASS] Added MonteCarloBarostat (1 bar, {TEMPERATURE_K:.0f} K)")

        # Save
        with open(f"{out_dir}/topology.pdb", "w") as f:
            PDBFile.writeFile(modeller.topology, modeller.positions, f)

        with open(f"{out_dir}/system.xml", "w") as f:
            f.write(XmlSerializer.serialize(system))

        # Save positions (handle both Vec3 list and numpy array formats)
        pos_unit = modeller.positions.value_in_unit(unit.nanometer)
        try:
            # Try Vec3 format first
            pos_array = np.array([[p.x, p.y, p.z] for p in pos_unit])
        except (AttributeError, TypeError):
            # Already numpy array or array-like
            pos_array = np.asarray(pos_unit)
        np.save(f"{out_dir}/positions.npy", pos_array)

        print(f"  [PASS] Saved {sys_name} topology, system, positions")

        built_systems[sys_name] = {
            'topology': modeller.topology,
            'positions': modeller.positions,
            'system': system,
            'n_atoms': modeller.topology.getNumAtoms()
        }

    print()

    # ============================================================
    # SECTION 6: CREATE ALCHEMICAL SYSTEMS WITH BORESCH RESTRAINTS
    # ============================================================
    print("="*60)
    print("[6] CREATING ALCHEMICAL SYSTEMS WITH BORESCH RESTRAINTS")
    print("="*60)

    if not HAS_ALCHEMY:
        print("[WARNING] OpenMMTools not available - skipping alchemical setup")
        print("  Install with: pip install openmmtools")
    else:
        # Create staged lambda schedule with restraints:
        # - Electrostatics: 1.0 -> 0.0 (turn off)
        # - Sterics: 1.0 -> 0.0 (turn off after electrostatics)
        # - Restraints: 0.0 -> 1.0 (turn ON as ligand decouples)
        N_ELEC = 10   # Windows for electrostatics
        N_VDW = 10    # Windows for sterics
        N_TOTAL = N_ELEC + N_VDW

        # Build lambda schedule arrays with 20 unique windows
        # Stage 1: Turn off electrostatics (1->0), turn on restraints (0->1)
        # Stage 2: Turn off sterics (0.9->0), restraints stay ON
        # Note: Stage 2 starts sterics at 0.9 (not 1.0) to avoid duplicate boundary state
        lambda_elec_arr = []
        lambda_sterics_arr = []
        lambda_restraints_arr = []

        # Stage 1: electrostatics 1.0 -> 0.0, sterics = 1.0, restraints 0.0 -> 1.0
        for i in range(N_ELEC):
            elec = 1.0 - (i / (N_ELEC - 1))      # 1.0 -> 0.0
            ster = 1.0                            # stays ON
            rest = i / (N_ELEC - 1)               # 0.0 -> 1.0
            lambda_elec_arr.append(elec)
            lambda_sterics_arr.append(ster)
            lambda_restraints_arr.append(rest)

        # Stage 2: electrostatics = 0.0, sterics 0.9 -> 0.0, restraints = 1.0
        # Sterics starts at 0.9 (not 1.0) to give 20 unique windows
        for i in range(N_VDW):
            elec = 0.0                            # stays OFF
            ster = 0.9 - (i * 0.9 / (N_VDW - 1)) # 0.9 -> 0.0
            rest = 1.0                            # stays ON
            lambda_elec_arr.append(elec)
            lambda_sterics_arr.append(ster)
            lambda_restraints_arr.append(rest)

        lambda_elec = np.array(lambda_elec_arr)
        lambda_sterics = np.array(lambda_sterics_arr)
        lambda_restraints = np.array(lambda_restraints_arr)
        N_TOTAL = len(lambda_elec)  # Should be 20

        print(f"[INFO] Creating {N_TOTAL} lambda windows (staged protocol with restraints)")
        print(f"  Stage 1: {N_ELEC} windows (0-{N_ELEC-1}) - electrostatics 1→0, restraints 0→1")
        print(f"  Stage 2: {N_TOTAL - N_ELEC} windows ({N_ELEC}-{N_TOTAL-1}) - sterics 1→0, restraints=1")
        print()
        print("  Lambda schedule:")
        print("  Window  Elec   Sterics  Restraints")
        for i in range(N_TOTAL):
            print(f"  {i:5d}  {lambda_elec[i]:.3f}   {lambda_sterics[i]:.3f}    {lambda_restraints[i]:.3f}")

        # Store anchor info for all complex systems
        all_anchors = {}
        # Reference anchor info from WT (residue IDs and atom names) for consistency
        wt_anchor_ref = None

        # Process systems in order: wt_complex first, then mut_complex, then solvent
        system_order = ['wt_complex', 'mut_complex', 'solvent']
        for sys_name in system_order:
            if sys_name not in built_systems:
                continue
            sys_data = built_systems[sys_name]

            print()
            print(f"[INFO] Alchemical setup for {sys_name}...")

            system = sys_data['system']
            topology = sys_data['topology']
            positions = sys_data['positions']
            has_protein = sys_name != 'solvent'

            # Find ligand atoms by matching residue name from trajectory (NO FALLBACK)
            ligand_atom_indices = []
            target_resname = ligand_res.name.upper()

            for res in topology.residues():
                if res.name.upper() == target_resname:
                    for atom in res.atoms():
                        ligand_atom_indices.append(atom.index)

            # Sort for deterministic ordering
            ligand_atom_indices = sorted(ligand_atom_indices)

            # Hard fail if ligand not found by exact residue name match
            if not ligand_atom_indices:
                raise RuntimeError(
                    f"No ligand atoms found in {sys_name}! "
                    f"Expected residue name '{target_resname}' from trajectory. "
                    f"Ensure ligand residue name matches exactly."
                )

            print(f"  Ligand atoms: {len(ligand_atom_indices)} (indices {ligand_atom_indices[0]}-{ligand_atom_indices[-1]})")

            # Add Boresch restraints for complex systems (not solvent-only)
            anchors = None
            if has_protein:
                print()
                print("  [INFO] Setting up Boresch restraints...")

                atom_list = list(topology.atoms())

                if sys_name == 'wt_complex' or wt_anchor_ref is None:
                    # Select anchors fresh for WT (or if no reference exists)
                    anchors, error = select_boresch_anchors(
                        topology, positions, ligand_atom_indices, cutoff_nm=1.0
                    )

                    if anchors is None:
                        print(f"  [FAIL] Could not select anchor atoms: {error}")
                        print("  [FAIL] Boresch restraints REQUIRED for complex FEP!")
                        sys.exit(1)

                    # Store reference info (residue IDs and atom names) for mutant
                    P3, P2, P1 = anchors['protein']
                    L1, L2, L3 = anchors['ligand']
                    wt_anchor_ref = {
                        'protein_resids': [atom_list[P3].residue.id, atom_list[P2].residue.id, atom_list[P1].residue.id],
                        'protein_names': [atom_list[P3].name, atom_list[P2].name, atom_list[P1].name],
                        'ligand_names': [atom_list[L1].name, atom_list[L2].name, atom_list[L3].name],
                    }
                    print(f"  [INFO] Selected anchors from {sys_name} (reference for mutant)")

                else:
                    # Use WT reference to find matching atoms in mutant
                    print(f"  [INFO] Using WT anchor reference for consistency...")

                    # Find protein anchors by residue ID and atom name
                    P3, P2, P1 = None, None, None
                    for atom in topology.atoms():
                        if atom.residue.id == wt_anchor_ref['protein_resids'][0] and atom.name == wt_anchor_ref['protein_names'][0]:
                            P3 = atom.index
                        elif atom.residue.id == wt_anchor_ref['protein_resids'][1] and atom.name == wt_anchor_ref['protein_names'][1]:
                            P2 = atom.index
                        elif atom.residue.id == wt_anchor_ref['protein_resids'][2] and atom.name == wt_anchor_ref['protein_names'][2]:
                            P1 = atom.index

                    # Find ligand anchors by atom name
                    L1, L2, L3 = None, None, None
                    for idx in ligand_atom_indices:
                        atom = atom_list[idx]
                        if atom.name == wt_anchor_ref['ligand_names'][0]:
                            L1 = idx
                        elif atom.name == wt_anchor_ref['ligand_names'][1]:
                            L2 = idx
                        elif atom.name == wt_anchor_ref['ligand_names'][2]:
                            L3 = idx

                    if None in [P3, P2, P1, L1, L2, L3]:
                        print(f"  [FAIL] Could not find matching anchor atoms in {sys_name}")
                        print(f"    P3={P3}, P2={P2}, P1={P1}, L1={L1}, L2={L2}, L3={L3}")
                        sys.exit(1)

                    # Recalculate equilibrium values from mutant positions
                    from openmm import unit as omm_unit
                    pos_unit = positions.value_in_unit(omm_unit.nanometer)
                    try:
                        pos_nm = np.array([[p.x, p.y, p.z] for p in pos_unit])
                    except (AttributeError, TypeError):
                        pos_nm = np.asarray(pos_unit)

                    def calc_dist(i, j):
                        return np.linalg.norm(pos_nm[i] - pos_nm[j])

                    def calc_angle(i, j, k):
                        v1 = pos_nm[i] - pos_nm[j]
                        v2 = pos_nm[k] - pos_nm[j]
                        cos_angle = np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2))
                        return np.arccos(np.clip(cos_angle, -1, 1))

                    def calc_dihedral(i, j, k, l):
                        b1 = pos_nm[j] - pos_nm[i]
                        b2 = pos_nm[k] - pos_nm[j]
                        b3 = pos_nm[l] - pos_nm[k]
                        n1 = np.cross(b1, b2)
                        n2 = np.cross(b2, b3)
                        m1 = np.cross(n1, b2 / np.linalg.norm(b2))
                        x = np.dot(n1, n2)
                        y = np.dot(m1, n2)
                        return np.arctan2(y, x)

                    anchors = {
                        'protein': [P3, P2, P1],
                        'ligand': [L1, L2, L3],
                        'values': {
                            'r_P1_L1': calc_dist(P1, L1),
                            'theta_P2_P1_L1': calc_angle(P2, P1, L1),
                            'theta_P1_L1_L2': calc_angle(P1, L1, L2),
                            'phi_P3_P2_P1_L1': calc_dihedral(P3, P2, P1, L1),
                            'phi_P2_P1_L1_L2': calc_dihedral(P2, P1, L1, L2),
                            'phi_P1_L1_L2_L3': calc_dihedral(P1, L1, L2, L3),
                        }
                    }
                    print(f"  [PASS] Matched anchors from WT reference")

                P3, P2, P1 = anchors['protein']
                L1, L2, L3 = anchors['ligand']

                print(f"  Protein anchors:")
                print(f"    P3: atom {P3} ({atom_list[P3].name} in {atom_list[P3].residue.name}{atom_list[P3].residue.id})")
                print(f"    P2: atom {P2} ({atom_list[P2].name} in {atom_list[P2].residue.name}{atom_list[P2].residue.id})")
                print(f"    P1: atom {P1} ({atom_list[P1].name} in {atom_list[P1].residue.name}{atom_list[P1].residue.id})")
                print(f"  Ligand anchors:")
                print(f"    L1: atom {L1} ({atom_list[L1].name})")
                print(f"    L2: atom {L2} ({atom_list[L2].name})")
                print(f"    L3: atom {L3} ({atom_list[L3].name})")

                print(f"  Equilibrium values:")
                vals = anchors['values']
                print(f"    r(P1-L1):      {vals['r_P1_L1']*10:.2f} Å")
                print(f"    θ(P2-P1-L1):   {np.degrees(vals['theta_P2_P1_L1']):.1f}°")
                print(f"    θ(P1-L1-L2):   {np.degrees(vals['theta_P1_L1_L2']):.1f}°")
                print(f"    φ(P3-P2-P1-L1): {np.degrees(vals['phi_P3_P2_P1_L1']):.1f}°")
                print(f"    φ(P2-P1-L1-L2): {np.degrees(vals['phi_P2_P1_L1_L2']):.1f}°")
                print(f"    φ(P1-L1-L2-L3): {np.degrees(vals['phi_P1_L1_L2_L3']):.1f}°")

                # Create and add restraint forces to system BEFORE alchemical factory
                restraint_forces = create_boresch_restraints(
                    anchors, k_distance=K_DISTANCE, k_angle=K_ANGLE, k_dihedral=K_DIHEDRAL
                )
                for name, force in restraint_forces:
                    system.addForce(force)
                print(f"  [PASS] Added {len(restraint_forces)} Boresch restraint forces")

                # Calculate standard state correction
                dG_corr = calc_boresch_correction(
                    anchors, temperature_K=TEMPERATURE_K,
                    k_r=K_DISTANCE, k_theta=K_ANGLE, k_phi=K_DIHEDRAL
                )
                print(f"  [INFO] Standard state correction: {dG_corr:.3f} kcal/mol")

                all_anchors[sys_name] = anchors

            # Create alchemical region (ANNIHILATION mode for both)
            # annihilate_electrostatics=True: required for proper PME handling
            # annihilate_sterics=True: consistent thermodynamic path (removes all ligand interactions)
            # Note: Annihilation removes ligand-ligand intramolecular interactions at lambda=0
            alchemical_region = AlchemicalRegion(
                alchemical_atoms=ligand_atom_indices,
                annihilate_electrostatics=True,   # Annihilate for proper PME
                annihilate_sterics=True,          # Annihilate for consistency
                softcore_alpha=0.5,
                softcore_a=1,
                softcore_b=1,
                softcore_c=6
            )

            # Create alchemical factory with direct-space PME treatment
            # Direct-space is required for alchemical electrostatics with PME
            factory = AbsoluteAlchemicalFactory(
                alchemical_pme_treatment='direct-space'
            )
            alchemical_system = factory.create_alchemical_system(system, alchemical_region)

            # Verify barostat is preserved in alchemical system (REQUIRED)
            from openmm import MonteCarloBarostat
            has_barostat = False
            for i in range(alchemical_system.getNumForces()):
                force = alchemical_system.getForce(i)
                if isinstance(force, MonteCarloBarostat):
                    has_barostat = True
                    break
            if has_barostat:
                print(f"  [PASS] Barostat preserved in alchemical system")
            else:
                raise RuntimeError(f"Barostat missing from alchemical system for {sys_name}! "
                                   "NPT equilibration requires barostat.")

            # Save alchemical system
            out_dir = f"fep_pmx/{sys_name}"
            with open(f"{out_dir}/alchemical_system.xml", "w") as f:
                f.write(XmlSerializer.serialize(alchemical_system))

            # Save lambda schedule (3 columns: elec, sterics, restraints)
            lambda_schedule = np.column_stack([lambda_elec, lambda_sterics, lambda_restraints])
            np.save(f"{out_dir}/lambda_schedule.npy", lambda_schedule)

            # Also save individual arrays for clarity
            np.save(f"{out_dir}/lambda_electrostatics.npy", lambda_elec)
            np.save(f"{out_dir}/lambda_sterics.npy", lambda_sterics)
            np.save(f"{out_dir}/lambda_restraints.npy", lambda_restraints)

            # Save anchor information for complex systems
            if anchors is not None:
                anchor_data = {
                    'protein_atoms': anchors['protein'],
                    'ligand_atoms': anchors['ligand'],
                    'equilibrium_values': anchors['values'],
                    'standard_state_correction_kcal_mol': calc_boresch_correction(
                        anchors, temperature_K=TEMPERATURE_K,
                        k_r=K_DISTANCE, k_theta=K_ANGLE, k_phi=K_DIHEDRAL
                    ),
                }
                np.save(f"{out_dir}/boresch_anchors.npy", anchor_data, allow_pickle=True)

            print(f"  [PASS] Alchemical system saved with {len(ligand_atom_indices)} alchemical atoms")
            print(f"  [PASS] Lambda schedule: {N_TOTAL} windows (elec+sterics+restraints)")

        # Verify Boresch anchors match between WT and mutant
        print()
        print("[INFO] Verifying Boresch anchor consistency...")
        if 'wt_complex' in all_anchors and 'mut_complex' in all_anchors:
            wt_anchors = all_anchors['wt_complex']
            mut_anchors = all_anchors['mut_complex']

            # Get atom names for comparison
            wt_atoms = list(built_systems['wt_complex']['topology'].atoms())
            mut_atoms = list(built_systems['mut_complex']['topology'].atoms())

            print("  WT anchor residues:")
            for i, idx in enumerate(wt_anchors['protein']):
                a = wt_atoms[idx]
                print(f"    P{3-i}: {a.residue.name}{a.residue.id} {a.name}")

            print("  Mutant anchor residues:")
            for i, idx in enumerate(mut_anchors['protein']):
                a = mut_atoms[idx]
                print(f"    P{3-i}: {a.residue.name}{a.residue.id} {a.name}")

            # Check if anchor residue IDs match
            wt_res_ids = [wt_atoms[idx].residue.id for idx in wt_anchors['protein']]
            mut_res_ids = [mut_atoms[idx].residue.id for idx in mut_anchors['protein']]

            if wt_res_ids == mut_res_ids:
                print("  [PASS] Anchor residue IDs match between WT and mutant")
            else:
                print(f"  [WARNING] Anchor residue IDs differ: WT={wt_res_ids}, Mut={mut_res_ids}")

        # Verify lambda direction
        print()
        print("[INFO] Lambda direction verification:")
        print(f"  Window 0:  elec={lambda_elec[0]:.3f}, sterics={lambda_sterics[0]:.3f}, restraints={lambda_restraints[0]:.3f}")
        print(f"  Window {N_TOTAL-1}: elec={lambda_elec[-1]:.3f}, sterics={lambda_sterics[-1]:.3f}, restraints={lambda_restraints[-1]:.3f}")
        if lambda_elec[0] == 1.0 and lambda_sterics[0] == 1.0 and lambda_restraints[0] == 0.0:
            print("  [PASS] Window 0 = fully coupled (correct)")
        else:
            print("  [WARNING] Window 0 values unexpected!")
        # Use approximate comparison for floating point
        if (np.isclose(lambda_elec[-1], 0.0) and np.isclose(lambda_sterics[-1], 0.0) and
            np.isclose(lambda_restraints[-1], 1.0)):
            print(f"  [PASS] Window {N_TOTAL-1} = fully decoupled with restraints ON (correct)")
        else:
            print(f"  [WARNING] Window {N_TOTAL-1} values unexpected!")

        # ============================================================
        # SECTION 6b: CREATE PER-WINDOW FOLDERS AND RUN SCRIPTS
        # ============================================================
        print()
        print("="*60)
        print("[6b] CREATING PER-WINDOW FOLDERS AND RUN SCRIPTS")
        print("="*60)

        # Template for run_window.py
        run_window_template = '''#!/usr/bin/env python
"""
FEP Window Runner - Window {window_idx}
System: {sys_name}
Lambda: elec={lambda_elec:.4f}, sterics={lambda_sterics:.4f}, restraints={lambda_restraints:.4f}

This script:
1. Loads the alchemical system
2. Runs NPT equilibration (with barostat)
3. Removes barostat for NVT production
4. Runs production MD and computes u_nk at all lambda states
5. Saves u_nk.npy for MBAR analysis
"""
import os
import sys
import numpy as np

# Add parent directory to path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from openmm import XmlSerializer, LangevinMiddleIntegrator, Platform, MonteCarloBarostat
from openmm.app import PDBFile, Simulation, StateDataReporter
from openmm import unit

# Configuration
WINDOW_IDX = {window_idx}
LAMBDA_ELEC = {lambda_elec}
LAMBDA_STERICS = {lambda_sterics}
LAMBDA_RESTRAINTS = {lambda_restraints}
SYS_NAME = "{sys_name}"
HAS_RESTRAINTS = {has_restraints}

# Simulation parameters
TEMPERATURE = {temperature} * unit.kelvin
FRICTION = 1.0 / unit.picosecond
TIMESTEP = 2.0 * unit.femtoseconds
EQUIL_STEPS = 50000    # 100 ps NPT equilibration
PROD_STEPS = 500000    # 1 ns production
ENERGY_INTERVAL = 500  # Compute energies every 1 ps

def main():
    print(f"="*60)
    print(f"FEP Window {{WINDOW_IDX}} - {{SYS_NAME}}")
    print(f"="*60)
    print(f"Lambda: elec={{LAMBDA_ELEC:.4f}}, sterics={{LAMBDA_STERICS:.4f}}, restraints={{LAMBDA_RESTRAINTS:.4f}}")
    print()

    # Load system and positions
    parent_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

    with open(os.path.join(parent_dir, "alchemical_system.xml"), "r") as f:
        system = XmlSerializer.deserialize(f.read())

    positions = np.load(os.path.join(parent_dir, "positions.npy"))
    positions = positions * unit.nanometer

    pdb = PDBFile(os.path.join(parent_dir, "topology.pdb"))
    topology = pdb.topology

    # Load lambda schedule
    lambda_schedule = np.load(os.path.join(parent_dir, "lambda_schedule.npy"))
    n_windows = len(lambda_schedule)
    print(f"Loaded lambda schedule: {{n_windows}} windows")

    # Create integrator
    integrator = LangevinMiddleIntegrator(TEMPERATURE, FRICTION, TIMESTEP)

    # Select platform (prefer CUDA)
    try:
        platform = Platform.getPlatformByName('CUDA')
        properties = {{'Precision': 'mixed'}}
        print("Using CUDA platform")
    except:
        try:
            platform = Platform.getPlatformByName('OpenCL')
            properties = {{}}
            print("Using OpenCL platform")
        except:
            platform = Platform.getPlatformByName('CPU')
            properties = {{}}
            print("Using CPU platform")

    # Create simulation
    simulation = Simulation(topology, system, integrator, platform, properties)
    simulation.context.setPositions(positions)

    # Set lambda parameters
    context = simulation.context
    context.setParameter('lambda_electrostatics', LAMBDA_ELEC)
    context.setParameter('lambda_sterics', LAMBDA_STERICS)
    if HAS_RESTRAINTS:
        context.setParameter('lambda_restraints', LAMBDA_RESTRAINTS)

    # Minimize
    print("Minimizing energy...")
    simulation.minimizeEnergy(maxIterations=1000)

    # NPT Equilibration (barostat already in system)
    print(f"Running NPT equilibration ({{EQUIL_STEPS}} steps)...")
    simulation.step(EQUIL_STEPS)

    # Remove barostat for NVT production
    print("Removing barostat for NVT production...")

    # Get state BEFORE modifying system (positions + box vectors)
    state = context.getState(getPositions=True)
    positions_after_equil = state.getPositions()
    box_vectors = state.getPeriodicBoxVectors()

    # Remove barostat from system
    for i in range(system.getNumForces()):
        force = system.getForce(i)
        if isinstance(force, MonteCarloBarostat):
            system.removeForce(i)
            break

    # Create fresh simulation with new integrator
    integrator2 = LangevinMiddleIntegrator(TEMPERATURE, FRICTION, TIMESTEP)
    simulation = Simulation(topology, system, integrator2, platform, properties)
    context = simulation.context

    # Restore box vectors FIRST, then positions
    context.setPeriodicBoxVectors(*box_vectors)
    context.setPositions(positions_after_equil)

    # Reinitialize velocities from Maxwell-Boltzmann (safer than copying)
    context.setVelocitiesToTemperature(TEMPERATURE)

    # Set lambda parameters
    context.setParameter('lambda_electrostatics', LAMBDA_ELEC)
    context.setParameter('lambda_sterics', LAMBDA_STERICS)
    if HAS_RESTRAINTS:
        context.setParameter('lambda_restraints', LAMBDA_RESTRAINTS)

    # Verify energy is sane before production
    state = context.getState(getEnergy=True)
    E_nvt = state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
    print(f"  NVT initial energy: {{E_nvt:.1f}} kJ/mol")
    if E_nvt > 0 or abs(E_nvt) > 1e10:
        raise RuntimeError(f"Energy looks wrong after NVT setup: {{E_nvt}}")

    # Production run - collect u_nk
    print(f"Running NVT production ({{PROD_STEPS}} steps)...")
    n_samples = PROD_STEPS // ENERGY_INTERVAL
    u_nk = np.zeros((n_samples, n_windows))

    for sample_idx in range(n_samples):
        # Run dynamics
        simulation.step(ENERGY_INTERVAL)

        # Compute reduced potential at all lambda states
        for k in range(n_windows):
            lam_e, lam_s, lam_r = lambda_schedule[k]
            context.setParameter('lambda_electrostatics', lam_e)
            context.setParameter('lambda_sterics', lam_s)
            if HAS_RESTRAINTS:
                context.setParameter('lambda_restraints', lam_r)

            state = context.getState(getEnergy=True)
            # Reduced potential u = U / kT
            kT = unit.MOLAR_GAS_CONSTANT_R * TEMPERATURE
            u_nk[sample_idx, k] = state.getPotentialEnergy() / kT

        # Reset to current window's lambda
        context.setParameter('lambda_electrostatics', LAMBDA_ELEC)
        context.setParameter('lambda_sterics', LAMBDA_STERICS)
        if HAS_RESTRAINTS:
            context.setParameter('lambda_restraints', LAMBDA_RESTRAINTS)

        if (sample_idx + 1) % 100 == 0:
            print(f"  Sample {{sample_idx + 1}}/{{n_samples}}")

    # Save u_nk
    np.save("u_nk.npy", u_nk)
    print(f"[PASS] Saved u_nk.npy: shape {{u_nk.shape}}")

    # Save final positions
    state = context.getState(getPositions=True)
    pos = state.getPositions(asNumpy=True).value_in_unit(unit.nanometer)
    np.save("final_positions.npy", pos)
    print("[PASS] Saved final_positions.npy")

    print()
    print("="*60)
    print("Window complete!")
    print("="*60)

if __name__ == "__main__":
    main()
'''

        # Create per-window folders for each system
        for sys_name in system_order:
            if sys_name not in built_systems:
                continue

            sys_dir = f"fep_pmx/{sys_name}"
            has_restraints = sys_name != 'solvent'

            print(f"[INFO] Creating window folders for {sys_name}...")

            for window_idx in range(N_TOTAL):
                window_dir = f"{sys_dir}/window_{window_idx:02d}"
                os.makedirs(window_dir, exist_ok=True)

                # Get lambda values for this window
                lam_e = lambda_elec[window_idx]
                lam_s = lambda_sterics[window_idx]
                lam_r = lambda_restraints[window_idx]

                # Write run_window.py
                script_content = run_window_template.format(
                    window_idx=window_idx,
                    sys_name=sys_name,
                    lambda_elec=lam_e,
                    lambda_sterics=lam_s,
                    lambda_restraints=lam_r,
                    has_restraints=has_restraints,
                    temperature=TEMPERATURE_K
                )

                with open(f"{window_dir}/run_window.py", "w") as f:
                    f.write(script_content)

            print(f"  [PASS] Created {N_TOTAL} window folders with run_window.py")

    print()

    # ============================================================
    # SECTION 7: SUMMARY
    # ============================================================
    print("="*60)
    print("SETUP SUMMARY")
    print("="*60)
    print()
    print("Systems created:")
    for sys_name, sys_data in built_systems.items():
        print(f"  {sys_name}: {sys_data['n_atoms']} atoms (solvated)")
    print()
    print("Files per system:")
    print("  topology.pdb           - Solvated structure")
    print("  system.xml             - OpenMM system (PME, solvated)")
    print("  positions.npy          - Coordinates")
    if HAS_ALCHEMY:
        print("  alchemical_system.xml  - Alchemical system for FEP")
        print("  lambda_schedule.npy    - Lambda schedule (20 windows, 3 cols: elec, sterics, restraints)")
        print("  lambda_*.npy           - Individual lambda arrays")
        print("  boresch_anchors.npy    - Anchor atoms and equilibrium values (complex systems)")
        print("  window_XX/             - Per-window folders (20 total)")
        print("    run_window.py        - Script to run this window and save u_nk.npy")
    print()
    print("Boresch restraints (complex systems only):")
    print("  - 6 restraints: 1 distance + 2 angles + 3 dihedrals")
    print("  - Force constants: k_r=10 kcal/mol/Å², k_θ=10 kcal/mol/rad², k_φ=10 kcal/mol/rad²")
    print("  - lambda_restraints: 0=OFF (coupled), 1=ON (decoupled)")
    print("  - Standard state correction MUST be applied after MBAR!")
    print()
    print("Mutation: R702W (ARG -> TRP)")
    print("  Charge change: +1 -> 0")
    print()
    if not mutation_successful:
        print("ACTION REQUIRED:")
        print("  Run: python mutate_r702w.py")
        print("  Then re-run this script")
        print()
    print("Next steps:")
    print("  1. Verify mutation in mut_complex")
    print("  2. Run NPT equilibration for each system:")
    print("     - lambda_elec=1, lambda_sterics=1, lambda_restraints=0 (fully coupled)")
    print("     - System has MonteCarloBarostat for NPT")
    print("  3. Run NVT production FEP at each lambda window:")
    print("     - REMOVE barostat for production (cleaner MBAR, no pV correction needed)")
    print("     - Set context.setParameter('lambda_electrostatics', value)")
    print("     - Set context.setParameter('lambda_sterics', value)")
    print("     - Set context.setParameter('lambda_restraints', value) [COMPLEX ONLY]")
    print("     NOTE: Solvent leg has no restraints - skip lambda_restraints for solvent")
    print("  4. Collect u_nk energies and analyze with MBAR")
    print("  5. Apply corrections to ΔG_bind (ABFE sign convention):")
    print("     ΔG_bind = ΔG_decouple_solvent - ΔG_decouple_complex + ΔG_Boresch + ΔG_PME")
    print("     - ΔG_Boresch: Standard state correction (~7.6 kcal/mol)")
    print("     - ΔG_PME: Charge correction (0 for neutral ligand, see calc_pme_charge_correction)")
    print("     (Negative ΔG_bind = favorable binding)")
    print()
    print("NPT vs NVT protocol:")
    print("  - Equilibration: NPT (with barostat) to relax box size")
    print("  - Production FEP: NVT (remove barostat) for clean MBAR analysis")
    print()
    print("Window count: 20 (K=20 for MBAR stacking)")
    print()
    print("="*60)


if __name__ == "__main__":
    main()
