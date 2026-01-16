#!/usr/bin/env python
"""
Deep verification of Natural Product CID_10120 FEP results.
Outputs to: fep_pmx_natural/verification_results.txt
"""
import numpy as np
import json
import sys
from pathlib import Path
from io import StringIO

# Constants
R = 8.314462618e-3  # kJ/(mol·K)
TEMPERATURE = 310.0  # K
RT = R * TEMPERATURE  # ~2.577 kJ/mol at 310K
KCAL_PER_KJ = 1 / 4.184

BASE_DIR = Path(__file__).parent.resolve()
OUTPUT_FILE = BASE_DIR / "verification_results.txt"


class OutputCapture:
    """Capture output to both file and stdout."""
    def __init__(self, filepath):
        self.filepath = filepath
        self.buffer = StringIO()

    def write(self, text):
        print(text, end='')
        self.buffer.write(text)

    def writeln(self, text=""):
        self.write(text + "\n")

    def save(self):
        with open(self.filepath, 'w') as f:
            f.write(self.buffer.getvalue())


def load_pdb_coordinates(pdb_file):
    """Load atom coordinates from PDB file."""
    atoms = []
    with open(pdb_file, 'r') as f:
        for line in f:
            if line.startswith(('ATOM', 'HETATM')):
                atom_name = line[12:16].strip()
                res_name = line[17:20].strip()
                chain = line[21:22].strip()
                res_id_str = line[22:26].strip()
                # Handle non-numeric residue IDs (like water "A000")
                try:
                    res_id = int(res_id_str)
                except ValueError:
                    res_id = -1  # Mark as non-standard
                try:
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                except ValueError:
                    continue  # Skip malformed lines
                atoms.append({
                    'name': atom_name,
                    'resname': res_name,
                    'chain': chain,
                    'resid': res_id,
                    'coords': np.array([x, y, z])
                })
    return atoms


def check_a_distances(out):
    """Check A: R702-to-ligand distance and protein dimensions."""
    out.writeln("\nA. DISTANCE CHECK")
    out.writeln("-" * 50)

    for sys_name in ['wt_complex', 'mut_complex']:
        pdb_file = BASE_DIR / sys_name / "topology.pdb"
        if not pdb_file.exists():
            out.writeln(f"  ERROR: {pdb_file} not found")
            continue

        atoms = load_pdb_coordinates(pdb_file)

        # Find R702 (or W702 for mutant) CA atom
        r702_ca = None
        mutation_residue = None
        for atom in atoms:
            if atom['resid'] == 702 and atom['name'] == 'CA':
                r702_ca = atom
                mutation_residue = atom['resname']
                break

        # Find ligand atoms (various possible residue names)
        ligand_resnames = ['UNK', 'NAT', 'LIG', 'MOL']
        ligand_atoms = [a for a in atoms if a['resname'] in ligand_resnames]

        # Also check for chain B residue 1 (common ligand location)
        if not ligand_atoms:
            ligand_atoms = [a for a in atoms if a['chain'] == 'B']

        if r702_ca is None:
            out.writeln(f"  {sys_name}: ERROR - R702 CA not found")
            continue

        if not ligand_atoms:
            out.writeln(f"  {sys_name}: ERROR - No ligand atoms found")
            continue

        # Calculate minimum distance to ligand
        min_dist = float('inf')
        closest_ligand_atom = None
        for lig_atom in ligand_atoms:
            dist = np.linalg.norm(r702_ca['coords'] - lig_atom['coords'])
            if dist < min_dist:
                min_dist = dist
                closest_ligand_atom = lig_atom

        # Calculate protein dimensions
        protein_atoms = [a for a in atoms if a['resname'] not in ligand_resnames + ['HOH', 'WAT', 'NA', 'CL', 'SOL']]
        if protein_atoms:
            coords = np.array([a['coords'] for a in protein_atoms])
            dims = coords.max(axis=0) - coords.min(axis=0)
        else:
            dims = np.array([0, 0, 0])

        out.writeln(f"  {sys_name}:")
        out.writeln(f"    Residue 702: {mutation_residue}")
        out.writeln(f"    R702 CA to ligand: {min_dist:.1f} A")
        out.writeln(f"    Closest ligand atom: {closest_ligand_atom['name']} ({closest_ligand_atom['resname']})")
        out.writeln(f"    Protein dimensions: {dims[0]:.1f} x {dims[1]:.1f} x {dims[2]:.1f} A")

    # Return average distance from both systems
    return min_dist if min_dist != float('inf') else None


def check_b_movement(out):
    """Check B: Protein movement during simulation."""
    out.writeln("\nB. PROTEIN MOVEMENT")
    out.writeln("-" * 50)

    results = {}

    for sys_name in ['wt_complex', 'mut_complex']:
        sys_dir = BASE_DIR / sys_name

        # Load initial positions
        init_pos_file = sys_dir / "positions.npy"
        if not init_pos_file.exists():
            out.writeln(f"  {sys_name}: ERROR - positions.npy not found")
            continue

        init_pos = np.load(init_pos_file)

        # Load topology to get CA indices
        pdb_file = sys_dir / "topology.pdb"
        atoms = load_pdb_coordinates(pdb_file)
        ca_indices = [i for i, a in enumerate(atoms) if a['name'] == 'CA' and a['resname'] not in ['UNK', 'NAT', 'LIG', 'HOH', 'WAT']]

        out.writeln(f"  {sys_name} (using {len(ca_indices)} CA atoms):")

        for window_idx in [0, 10, 19]:
            final_pos_file = sys_dir / f"window_{window_idx:02d}" / "final_positions.npy"
            if not final_pos_file.exists():
                out.writeln(f"    Window {window_idx}: MISSING")
                continue

            final_pos = np.load(final_pos_file)

            # Calculate RMSD for CA atoms
            if len(init_pos) != len(final_pos):
                out.writeln(f"    Window {window_idx}: Size mismatch ({len(init_pos)} vs {len(final_pos)})")
                continue

            # Extract CA coordinates (convert nm to Angstrom)
            init_ca = init_pos[ca_indices] * 10  # nm to A
            final_ca = final_pos[ca_indices] * 10

            # Simple RMSD (no alignment)
            diff = final_ca - init_ca
            rmsd = np.sqrt(np.mean(np.sum(diff**2, axis=1)))

            status = "OK" if rmsd > 1.0 else "LOW" if rmsd > 0.5 else "FROZEN"
            out.writeln(f"    Window {window_idx}: RMSD = {rmsd:.2f} A [{status}]")

            results[f"{sys_name}_w{window_idx}"] = rmsd

    return results


def load_u_nk_with_cap(sys_dir, n_windows=20, max_u=1e6):
    """Load u_nk data with value capping."""
    u_nk_list = []
    n_samples_list = []
    cap_counts = []

    for win_idx in range(n_windows):
        u_nk_file = sys_dir / f"window_{win_idx:02d}" / "u_nk.npy"
        if not u_nk_file.exists():
            raise FileNotFoundError(f"Missing: {u_nk_file}")

        u_nk = np.load(u_nk_file)

        # Validate shape
        if u_nk.ndim != 2 or u_nk.shape[1] != n_windows:
            raise ValueError(f"Window {win_idx}: bad shape {u_nk.shape}, expected (N, {n_windows})")

        # Check for NaN/Inf
        if np.isnan(u_nk).any() or np.isinf(u_nk).any():
            raise ValueError(f"Window {win_idx}: contains NaN or Inf values")
        u_min = u_nk.min()
        n_capped = np.sum(u_nk > u_min + max_u)

        if n_capped > 0:
            u_nk = np.clip(u_nk, None, u_min + max_u)

        u_nk_list.append(u_nk)
        n_samples_list.append(u_nk.shape[0])
        cap_counts.append(n_capped)

    return u_nk_list, np.array(n_samples_list), cap_counts


def run_mbar_subset(u_nk_list, N_k, start_win=0, end_win=None):
    """Run MBAR on a subset of windows."""
    from pymbar import MBAR

    if end_win is None:
        end_win = len(u_nk_list)

    # Extract subset
    subset_u_nk = u_nk_list[start_win:end_win]
    subset_N_k = N_k[start_win:end_win]
    n_subset = len(subset_u_nk)

    # Stack and transpose
    u_kn_stacked = np.vstack(subset_u_nk)
    # Only use columns for the subset windows
    u_kn_subset = u_kn_stacked[:, start_win:end_win]
    u_kn = u_kn_subset.T

    mbar = MBAR(u_kn, subset_N_k, verbose=False)
    results = mbar.compute_free_energy_differences()

    dG_kT = results['Delta_f'][0, -1]
    dG_err_kT = results['dDelta_f'][0, -1]

    dG_kcal = dG_kT * RT * KCAL_PER_KJ
    dG_err_kcal = dG_err_kT * RT * KCAL_PER_KJ

    return dG_kcal, dG_err_kcal, mbar


def check_c_phase_breakdown(out):
    """Check C: Per-phase dG breakdown."""
    out.writeln("\nC. PER-PHASE dG BREAKDOWN (kcal/mol)")
    out.writeln("-" * 50)

    results = {}

    for sys_name in ['wt_complex', 'mut_complex']:
        sys_dir = BASE_DIR / sys_name

        try:
            u_nk_list, N_k, _ = load_u_nk_with_cap(sys_dir)

            # Full analysis
            dG_total, dG_total_err, _ = run_mbar_subset(u_nk_list, N_k, 0, 20)

            # Phase 1: Electrostatics (windows 0-9)
            dG_elec, dG_elec_err, _ = run_mbar_subset(u_nk_list, N_k, 0, 10)

            # Phase 2: Sterics (windows 10-19)
            dG_sterics, dG_sterics_err, _ = run_mbar_subset(u_nk_list, N_k, 10, 20)

            results[sys_name] = {
                'elec': (dG_elec, dG_elec_err),
                'sterics': (dG_sterics, dG_sterics_err),
                'total': (dG_total, dG_total_err)
            }

            out.writeln(f"  {sys_name}:")
            out.writeln(f"    Electrostatics (0-9):  {dG_elec:+.2f} +/- {dG_elec_err:.2f}")
            out.writeln(f"    Sterics (10-19):       {dG_sterics:+.2f} +/- {dG_sterics_err:.2f}")
            out.writeln(f"    Total (0-19):          {dG_total:+.2f} +/- {dG_total_err:.2f}")

        except Exception as e:
            out.writeln(f"  {sys_name}: ERROR - {e}")

    # Calculate differences
    if 'wt_complex' in results and 'mut_complex' in results:
        out.writeln(f"\n  Phase contributions to ddG:")
        ddG_elec = results['mut_complex']['elec'][0] - results['wt_complex']['elec'][0]
        ddG_sterics = results['mut_complex']['sterics'][0] - results['wt_complex']['sterics'][0]
        ddG_total = results['mut_complex']['total'][0] - results['wt_complex']['total'][0]
        out.writeln(f"    Electrostatics: {ddG_elec:+.2f} kcal/mol")
        out.writeln(f"    Sterics:        {ddG_sterics:+.2f} kcal/mol")
        out.writeln(f"    Total ddG:      {ddG_total:+.2f} kcal/mol")

    return results


def check_d_window_consistency(out):
    """Check D: Windows 0-14 vs 15-19 consistency."""
    out.writeln("\nD. WINDOW CONSISTENCY")
    out.writeln("-" * 50)

    results = {}

    for sys_name in ['wt_complex', 'mut_complex']:
        sys_dir = BASE_DIR / sys_name

        try:
            u_nk_list, N_k, _ = load_u_nk_with_cap(sys_dir)

            # Windows 0-14
            dG_0_14, err_0_14, _ = run_mbar_subset(u_nk_list, N_k, 0, 15)

            # Windows 15-19
            dG_15_19, err_15_19, _ = run_mbar_subset(u_nk_list, N_k, 15, 20)

            # Full
            dG_full, err_full, _ = run_mbar_subset(u_nk_list, N_k, 0, 20)

            results[sys_name] = {
                '0-14': dG_0_14,
                '15-19': dG_15_19,
                'full': dG_full
            }

            out.writeln(f"  {sys_name}:")
            out.writeln(f"    Windows 0-14:  {dG_0_14:+.2f} +/- {err_0_14:.2f} kcal/mol")
            out.writeln(f"    Windows 15-19: {dG_15_19:+.2f} +/- {err_15_19:.2f} kcal/mol")
            out.writeln(f"    Full (0-19):   {dG_full:+.2f} +/- {err_full:.2f} kcal/mol")

        except Exception as e:
            out.writeln(f"  {sys_name}: ERROR - {e}")

    return results


def check_e_bootstrap(out, n_bootstrap=50):
    """Check E: Bootstrap analysis."""
    out.writeln("\nE. BOOTSTRAP ANALYSIS")
    out.writeln("-" * 50)
    out.writeln(f"  Running {n_bootstrap} bootstrap iterations...")

    from pymbar import MBAR

    results = {}

    for sys_name in ['wt_complex', 'mut_complex']:
        sys_dir = BASE_DIR / sys_name

        try:
            u_nk_list, N_k, _ = load_u_nk_with_cap(sys_dir)

            # Original MBAR
            dG_orig, dG_err_orig, _ = run_mbar_subset(u_nk_list, N_k)

            # Bootstrap
            bootstrap_dGs = []
            for i in range(n_bootstrap):
                # Resample within each window
                resampled = []
                for u_nk in u_nk_list:
                    n_samples = len(u_nk)
                    idx = np.random.choice(n_samples, size=n_samples, replace=True)
                    resampled.append(u_nk[idx])

                try:
                    dG_boot, _, _ = run_mbar_subset(resampled, N_k)
                    bootstrap_dGs.append(dG_boot)
                except:
                    pass  # Skip failed bootstraps

            bootstrap_dGs = np.array(bootstrap_dGs)
            boot_mean = np.mean(bootstrap_dGs)
            boot_std = np.std(bootstrap_dGs)
            boot_ci = np.percentile(bootstrap_dGs, [2.5, 97.5])

            ratio = boot_std / dG_err_orig if dG_err_orig > 0 else float('inf')

            results[sys_name] = {
                'mbar_dG': dG_orig,
                'mbar_err': dG_err_orig,
                'boot_mean': boot_mean,
                'boot_std': boot_std,
                'boot_ci': boot_ci,
                'ratio': ratio
            }

            out.writeln(f"  {sys_name}:")
            out.writeln(f"    MBAR dG:       {dG_orig:+.2f} +/- {dG_err_orig:.2f} kcal/mol")
            out.writeln(f"    Bootstrap:     {boot_mean:+.2f} +/- {boot_std:.2f} kcal/mol")
            out.writeln(f"    95% CI:        [{boot_ci[0]:.2f}, {boot_ci[1]:.2f}]")
            out.writeln(f"    Ratio (boot/MBAR): {ratio:.2f} {'[OK]' if ratio < 1.5 else '[HIGH]'}")

        except Exception as e:
            out.writeln(f"  {sys_name}: ERROR - {e}")

    # DDG bootstrap
    if 'wt_complex' in results and 'mut_complex' in results:
        ddG = results['mut_complex']['mbar_dG'] - results['wt_complex']['mbar_dG']
        ddG_err = np.sqrt(results['wt_complex']['mbar_err']**2 + results['mut_complex']['mbar_err']**2)
        ddG_boot_std = np.sqrt(results['wt_complex']['boot_std']**2 + results['mut_complex']['boot_std']**2)

        out.writeln(f"\n  ddG Analysis:")
        out.writeln(f"    MBAR ddG:      {ddG:+.2f} +/- {ddG_err:.2f} kcal/mol")
        out.writeln(f"    Bootstrap std: +/- {ddG_boot_std:.2f} kcal/mol")

    return results


def check_f_capping(out):
    """Check F: Value capping impact analysis."""
    out.writeln("\nF. VALUE CAPPING IMPACT")
    out.writeln("-" * 50)

    results = {}

    for sys_name in ['wt_complex', 'mut_complex']:
        sys_dir = BASE_DIR / sys_name

        out.writeln(f"  {sys_name}:")

        # Count overflow values per window
        out.writeln(f"    Overflow counts (at cap=1e6):")
        u_nk_list, N_k, cap_counts = load_u_nk_with_cap(sys_dir, max_u=1e6)
        total_capped = sum(cap_counts)
        for i, count in enumerate(cap_counts):
            if count > 0:
                out.writeln(f"      Window {i}: {count} values capped")
        out.writeln(f"      Total: {total_capped}")

        # Sensitivity test
        out.writeln(f"    Sensitivity test:")
        dG_results = {}
        for cap_val in [1e5, 1e6, 1e7]:
            try:
                u_nk_list, N_k, _ = load_u_nk_with_cap(sys_dir, max_u=cap_val)
                dG, _, _ = run_mbar_subset(u_nk_list, N_k)
                dG_results[cap_val] = dG
                out.writeln(f"      Cap {cap_val:.0e}: dG = {dG:+.2f} kcal/mol")
            except Exception as e:
                out.writeln(f"      Cap {cap_val:.0e}: FAILED - {e}")

        # Test excluding window 19
        try:
            u_nk_list, N_k, _ = load_u_nk_with_cap(sys_dir)
            dG_no19, _, _ = run_mbar_subset(u_nk_list, N_k, 0, 19)
            out.writeln(f"      Excluding win19: dG = {dG_no19:+.2f} kcal/mol")
            dG_results['no_win19'] = dG_no19
        except Exception as e:
            out.writeln(f"      Excluding win19: FAILED - {e}")

        results[sys_name] = dG_results

    # DDG sensitivity
    if 'wt_complex' in results and 'mut_complex' in results:
        out.writeln(f"\n  ddG sensitivity:")
        for key in [1e5, 1e6, 1e7, 'no_win19']:
            if key in results['wt_complex'] and key in results['mut_complex']:
                ddG = results['mut_complex'][key] - results['wt_complex'][key]
                label = f"Cap {key:.0e}" if isinstance(key, float) else "No win19"
                out.writeln(f"    {label}: ddG = {ddG:+.2f} kcal/mol")

    return results


def check_g_overlap(out):
    """Check G: Overlap matrix analysis."""
    out.writeln("\nG. OVERLAP MATRIX")
    out.writeln("-" * 50)

    from pymbar import MBAR

    results = {}

    for sys_name in ['wt_complex', 'mut_complex']:
        sys_dir = BASE_DIR / sys_name

        try:
            u_nk_list, N_k, _ = load_u_nk_with_cap(sys_dir)
            _, _, mbar = run_mbar_subset(u_nk_list, N_k)

            overlap_results = mbar.compute_overlap()
            overlap_matrix = overlap_results['matrix']

            out.writeln(f"  {sys_name}:")

            # Find min adjacent overlap
            min_overlap = 1.0
            min_pair = (0, 1)
            low_overlaps = []

            for i in range(len(overlap_matrix) - 1):
                ov = overlap_matrix[i, i+1]
                if ov < min_overlap:
                    min_overlap = ov
                    min_pair = (i, i+1)
                if ov < 0.03:
                    low_overlaps.append((i, i+1, ov))

            out.writeln(f"    Min adjacent overlap: {min_overlap:.4f} at windows {min_pair[0]}-{min_pair[1]}")

            if low_overlaps:
                out.writeln(f"    LOW OVERLAP PAIRS (<3%):")
                for i, j, ov in low_overlaps:
                    out.writeln(f"      {i}-{j}: {ov:.4f}")
            else:
                out.writeln(f"    All adjacent overlaps > 3%: OK")

            # Print condensed matrix (diagonal neighbors only)
            out.writeln(f"    Adjacent overlaps:")
            row1 = "      "
            for i in range(19):
                row1 += f"{overlap_matrix[i, i+1]:.2f} "
            out.writeln(row1)

            results[sys_name] = {
                'matrix': overlap_matrix,
                'min_overlap': min_overlap,
                'min_pair': min_pair,
                'low_pairs': low_overlaps
            }

        except Exception as e:
            out.writeln(f"  {sys_name}: ERROR - {e}")

    return results


def main():
    out = OutputCapture(OUTPUT_FILE)

    out.writeln("=" * 80)
    out.writeln("NATURAL PRODUCT CID_10120 (BUFADIENOLIDE) FEP VERIFICATION")
    out.writeln("=" * 80)
    out.writeln(f"Base directory: {BASE_DIR}")
    out.writeln(f"Output file: {OUTPUT_FILE}")

    # Run all checks
    try:
        dist = check_a_distances(out)
    except Exception as e:
        out.writeln(f"\nA. DISTANCE CHECK: ERROR - {e}")
        dist = None

    try:
        movement = check_b_movement(out)
    except Exception as e:
        out.writeln(f"\nB. PROTEIN MOVEMENT: ERROR - {e}")
        movement = None

    try:
        phases = check_c_phase_breakdown(out)
    except Exception as e:
        out.writeln(f"\nC. PER-PHASE BREAKDOWN: ERROR - {e}")
        phases = None

    try:
        consistency = check_d_window_consistency(out)
    except Exception as e:
        out.writeln(f"\nD. WINDOW CONSISTENCY: ERROR - {e}")
        consistency = None

    try:
        bootstrap = check_e_bootstrap(out, n_bootstrap=50)
    except Exception as e:
        out.writeln(f"\nE. BOOTSTRAP ANALYSIS: ERROR - {e}")
        bootstrap = None

    try:
        capping = check_f_capping(out)
    except Exception as e:
        out.writeln(f"\nF. VALUE CAPPING: ERROR - {e}")
        capping = None

    try:
        overlap = check_g_overlap(out)
    except Exception as e:
        out.writeln(f"\nG. OVERLAP MATRIX: ERROR - {e}")
        overlap = None

    # Summary
    out.writeln("\n" + "=" * 80)
    out.writeln("VERIFICATION SUMMARY")
    out.writeln("=" * 80)

    checks = []

    # A: Distance
    if dist and dist > 50:
        checks.append(("[X] A. Distance confirmed >50 A (allosteric)", True))
    else:
        checks.append(("[ ] A. Distance NOT confirmed allosteric", False))

    # B: Movement
    if movement:
        all_moved = all(v > 1.0 for v in movement.values())
        if all_moved:
            checks.append(("[X] B. Protein RMSD > 1 A (not frozen)", True))
        else:
            checks.append(("[ ] B. Some windows show low RMSD", False))
    else:
        checks.append(("[ ] B. Movement check failed", False))

    # C: Phase breakdown
    if phases:
        checks.append(("[X] C. Per-phase breakdown computed", True))
    else:
        checks.append(("[ ] C. Phase breakdown failed", False))

    # D: Consistency
    if consistency:
        checks.append(("[X] D. Window consistency computed", True))
    else:
        checks.append(("[ ] D. Consistency check failed", False))

    # E: Bootstrap
    if bootstrap:
        wt_ratio = bootstrap.get('wt_complex', {}).get('ratio', 999)
        mut_ratio = bootstrap.get('mut_complex', {}).get('ratio', 999)
        if wt_ratio < 1.5 and mut_ratio < 1.5:
            checks.append(("[X] E. Bootstrap agrees with MBAR error", True))
        else:
            checks.append(("[ ] E. Bootstrap shows higher variance than MBAR", False))
    else:
        checks.append(("[ ] E. Bootstrap analysis failed", False))

    # F: Capping
    if capping:
        checks.append(("[X] F. Capping sensitivity analyzed", True))
    else:
        checks.append(("[ ] F. Capping analysis failed", False))

    # G: Overlap
    if overlap:
        wt_min = overlap.get('wt_complex', {}).get('min_overlap', 0)
        mut_min = overlap.get('mut_complex', {}).get('min_overlap', 0)
        if wt_min > 0.03 and mut_min > 0.03:
            checks.append(("[X] G. All overlaps > 3%", True))
        else:
            checks.append(("[ ] G. Some overlaps < 3%", False))
    else:
        checks.append(("[ ] G. Overlap analysis failed", False))

    for check_str, _ in checks:
        out.writeln(check_str)

    passed = sum(1 for _, ok in checks if ok)
    total = len(checks)

    out.writeln(f"\nPassed: {passed}/{total}")

    if passed == total:
        out.writeln("\nOVERALL ASSESSMENT: RELIABLE")
    elif passed >= total - 2:
        out.writeln("\nOVERALL ASSESSMENT: MARGINAL")
    else:
        out.writeln("\nOVERALL ASSESSMENT: SUSPECT")

    out.writeln("\n" + "=" * 80)

    # Save output
    out.save()
    print(f"\nResults saved to: {OUTPUT_FILE}")


if __name__ == "__main__":
    main()
