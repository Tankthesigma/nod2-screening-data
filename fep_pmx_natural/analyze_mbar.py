#!/usr/bin/env python
"""
MBAR Analysis for Natural Product CID_10120 (Bufadienolide) FEP
Computes ΔG for WT and MUT complexes, then calculates ΔΔG for R702W mutation effect.
"""
import numpy as np
import json
from pathlib import Path
from pymbar import MBAR

# Constants
R = 8.314462618e-3  # kJ/(mol·K)
TEMPERATURE = 310.0  # K
RT = R * TEMPERATURE  # ~2.577 kJ/mol at 310K

BASE_DIR = Path(__file__).parent.resolve()


def load_u_nk_data(sys_dir, n_windows=20, cap_extreme=True, max_u=1e6):
    """Load u_nk data from all windows for a system.

    Args:
        sys_dir: Directory containing window subdirectories
        n_windows: Number of windows to load
        cap_extreme: If True, cap extreme values to prevent overflow
        max_u: Maximum allowed value (in kT units)
    """
    u_nk_list = []
    n_samples_list = []
    n_capped_total = 0

    for win_idx in range(n_windows):
        win_dir = sys_dir / f"window_{win_idx:02d}"
        u_nk_file = win_dir / "u_nk.npy"

        if not u_nk_file.exists():
            raise FileNotFoundError(f"Missing u_nk file: {u_nk_file}")

        u_nk = np.load(u_nk_file)

        # Validate shape
        if u_nk.ndim != 2:
            raise ValueError(f"u_nk has wrong dimensions: {u_nk.ndim}, expected 2")
        if u_nk.shape[1] != n_windows:
            raise ValueError(f"u_nk has {u_nk.shape[1]} columns, expected {n_windows}")

        # Check for NaN/Inf
        if np.isnan(u_nk).any():
            nan_count = np.isnan(u_nk).sum()
            raise ValueError(f"Window {win_idx} has {nan_count} NaN values")
        if np.isinf(u_nk).any():
            inf_count = np.isinf(u_nk).sum()
            raise ValueError(f"Window {win_idx} has {inf_count} Inf values")

        # Cap extreme values to prevent MBAR overflow
        if cap_extreme:
            # Find the minimum (most favorable) energy in this window
            u_min = u_nk.min()
            # Cap values that are more than max_u above the minimum
            n_capped = np.sum(u_nk > u_min + max_u)
            if n_capped > 0:
                u_nk = np.clip(u_nk, None, u_min + max_u)
                n_capped_total += n_capped
                print(f"  Window {win_idx}: capped {n_capped} extreme values")

        u_nk_list.append(u_nk)
        n_samples_list.append(u_nk.shape[0])

    if n_capped_total > 0:
        print(f"  Total capped values: {n_capped_total}")

    return u_nk_list, n_samples_list


def run_mbar_analysis(sys_name, sys_dir, n_windows=20):
    """Run MBAR analysis for a single system."""
    print(f"\n{'='*60}")
    print(f"MBAR Analysis: {sys_name}")
    print(f"{'='*60}")

    # Load data
    print("Loading u_nk data...")
    u_nk_list, n_samples_list = load_u_nk_data(sys_dir, n_windows)

    N_k = np.array(n_samples_list)
    print(f"  Windows: {n_windows}")
    print(f"  Samples per window: {N_k}")
    print(f"  Total samples: {N_k.sum()}")

    # Stack u_nk matrices: shape (total_samples, n_windows)
    u_kn_stacked = np.vstack(u_nk_list)
    print(f"  Combined u_kn shape: {u_kn_stacked.shape}")

    # MBAR expects u_kn with shape (K, N) where K=states, N=total samples
    # Our u_nk has shape (N, K), so transpose
    u_kn = u_kn_stacked.T
    print(f"  Transposed for MBAR: {u_kn.shape}")

    # Run MBAR
    print("\nRunning MBAR...")
    try:
        mbar = MBAR(u_kn, N_k, verbose=False)
    except Exception as e:
        print(f"  MBAR failed: {e}")
        raise

    # Get free energy differences
    results = mbar.compute_free_energy_differences()
    f_ij = results['Delta_f']
    df_ij = results['dDelta_f']

    # Total ΔG from state 0 (fully coupled) to state K-1 (fully decoupled)
    dG = f_ij[0, -1] * RT  # Convert from kT to kJ/mol
    dG_err = df_ij[0, -1] * RT

    print(f"\n  ΔG (0 → {n_windows-1}): {dG:.2f} ± {dG_err:.2f} kJ/mol")
    print(f"  ΔG (0 → {n_windows-1}): {dG/4.184:.2f} ± {dG_err/4.184:.2f} kcal/mol")

    # Compute overlap matrix
    print("\nComputing overlap matrix...")
    overlap_results = mbar.compute_overlap()
    overlap_matrix = overlap_results['matrix']

    # Check adjacent window overlaps
    print("\nAdjacent window overlaps:")
    min_overlap = 1.0
    for i in range(n_windows - 1):
        ov = overlap_matrix[i, i+1]
        status = "OK" if ov > 0.03 else "LOW" if ov > 0.01 else "POOR"
        print(f"  {i} ↔ {i+1}: {ov:.4f} [{status}]")
        min_overlap = min(min_overlap, ov)

    print(f"\n  Minimum overlap: {min_overlap:.4f}")
    if min_overlap < 0.03:
        print("  WARNING: Low overlap may affect accuracy")

    return {
        'dG_kJ': dG,
        'dG_err_kJ': dG_err,
        'dG_kcal': dG / 4.184,
        'dG_err_kcal': dG_err / 4.184,
        'min_overlap': min_overlap,
        'f_ij': f_ij,
        'df_ij': df_ij,
        'overlap_matrix': overlap_matrix
    }


def load_boresch_correction(sys_dir):
    """Load Boresch restraint correction from JSON file."""
    json_file = sys_dir / "boresch_anchors.json"
    if not json_file.exists():
        print(f"  Warning: No Boresch anchors file at {json_file}")
        return 0.0

    with open(json_file, 'r') as f:
        data = json.load(f)

    # Extract equilibrium values
    eq = data.get('equilibrium', {})
    r0 = eq.get('r0_nm', 0.4)  # nm

    # Standard state correction for Boresch restraints
    # ΔG_restraint = -RT * ln(V_site / V_standard)
    # V_standard = 1660 Å³ (1 M standard state)
    # V_site ≈ 4π * r0³ / 3 for spherical approximation

    V_standard = 1.660  # nm³ (1660 Å³)
    V_site = 4 * np.pi * r0**3 / 3  # nm³

    # This is approximate - full Boresch correction is more complex
    dG_restraint = -RT * np.log(V_site / V_standard)

    print(f"  Boresch r0 = {r0:.3f} nm")
    print(f"  V_site ≈ {V_site:.4f} nm³")
    print(f"  ΔG_restraint ≈ {dG_restraint:.2f} kJ/mol")

    return dG_restraint


def main():
    print("="*70)
    print("MBAR ANALYSIS - Natural Product CID_10120 (Bufadienolide)")
    print("NOD2 R702W Mutation Effect on Binding")
    print("="*70)

    results = {}

    # Analyze WT complex
    wt_dir = BASE_DIR / "wt_complex"
    if not wt_dir.exists():
        print(f"ERROR: WT directory not found: {wt_dir}")
        return

    results['wt'] = run_mbar_analysis("WT Complex", wt_dir)

    # Analyze MUT complex
    mut_dir = BASE_DIR / "mut_complex"
    if not mut_dir.exists():
        print(f"ERROR: MUT directory not found: {mut_dir}")
        return

    results['mut'] = run_mbar_analysis("MUT Complex (R702W)", mut_dir)

    # Load Boresch corrections
    print(f"\n{'='*60}")
    print("Boresch Restraint Corrections")
    print(f"{'='*60}")

    print("\nWT Complex:")
    dG_restraint_wt = load_boresch_correction(wt_dir)

    print("\nMUT Complex:")
    dG_restraint_mut = load_boresch_correction(mut_dir)

    # Calculate ΔΔG
    print(f"\n{'='*70}")
    print("FINAL RESULTS")
    print(f"{'='*70}")

    dG_wt = results['wt']['dG_kcal']
    dG_wt_err = results['wt']['dG_err_kcal']
    dG_mut = results['mut']['dG_kcal']
    dG_mut_err = results['mut']['dG_err_kcal']

    # ΔΔG = ΔG_mut - ΔG_wt
    # Positive ΔΔG means mutation weakens binding
    ddG = dG_mut - dG_wt
    ddG_err = np.sqrt(dG_wt_err**2 + dG_mut_err**2)

    print(f"\nΔG (WT complex):  {dG_wt:+.2f} ± {dG_wt_err:.2f} kcal/mol")
    print(f"ΔG (MUT complex): {dG_mut:+.2f} ± {dG_mut_err:.2f} kcal/mol")
    print(f"\n{'─'*50}")
    print(f"ΔΔG (R702W effect): {ddG:+.2f} ± {ddG_err:.2f} kcal/mol")
    print(f"{'─'*50}")

    if ddG > 0:
        print(f"\n→ R702W mutation WEAKENS binding by {abs(ddG):.2f} kcal/mol")
    else:
        print(f"\n→ R702W mutation STRENGTHENS binding by {abs(ddG):.2f} kcal/mol")

    # Quality assessment
    print(f"\n{'='*60}")
    print("Quality Assessment")
    print(f"{'='*60}")
    wt_overlap = results['wt']['min_overlap']
    mut_overlap = results['mut']['min_overlap']
    print(f"  WT min overlap:  {wt_overlap:.4f} {'✓' if wt_overlap > 0.03 else '⚠'}")
    print(f"  MUT min overlap: {mut_overlap:.4f} {'✓' if mut_overlap > 0.03 else '⚠'}")

    if wt_overlap > 0.03 and mut_overlap > 0.03:
        print("\n  ✓ Good overlap - results are reliable")
    else:
        print("\n  ⚠ Low overlap detected - consider more sampling or windows")

    # Save results
    results_file = BASE_DIR / "mbar_results.json"
    save_data = {
        'ligand': 'CID_10120 (Bufadienolide)',
        'temperature_K': TEMPERATURE,
        'wt_dG_kcal': float(dG_wt),
        'wt_dG_err_kcal': float(dG_wt_err),
        'mut_dG_kcal': float(dG_mut),
        'mut_dG_err_kcal': float(dG_mut_err),
        'ddG_kcal': float(ddG),
        'ddG_err_kcal': float(ddG_err),
        'wt_min_overlap': float(wt_overlap),
        'mut_min_overlap': float(mut_overlap)
    }
    with open(results_file, 'w') as f:
        json.dump(save_data, f, indent=2)
    print(f"\nResults saved to: {results_file}")

    return results


if __name__ == "__main__":
    main()
