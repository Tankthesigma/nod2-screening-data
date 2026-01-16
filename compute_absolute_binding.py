#!/usr/bin/env python
"""
Compute absolute binding free energies from FEP data.

ΔG_bind = ΔG_complex - ΔG_solvent

For each compound:
- Load complex WT and MUT results
- Run MBAR on solvent leg
- Compute absolute ΔG_bind
"""
import numpy as np
import json
from pathlib import Path
from pymbar import MBAR

# Constants
R = 8.314462618e-3  # kJ/(mol·K)
TEMPERATURE = 310.0  # K
RT = R * TEMPERATURE
KCAL_PER_KJ = 1 / 4.184


def load_u_nk_with_cap(sys_dir, n_windows=20, max_u=1e6):
    """Load u_nk data with value capping to handle overflow."""
    u_nk_list = []
    n_samples_list = []
    n_capped_total = 0

    for win_idx in range(n_windows):
        u_nk_file = sys_dir / f"window_{win_idx:02d}" / "u_nk.npy"
        if not u_nk_file.exists():
            raise FileNotFoundError(f"Missing: {u_nk_file}")

        u_nk = np.load(u_nk_file)

        # Validate
        if u_nk.ndim != 2 or u_nk.shape[1] != n_windows:
            raise ValueError(f"Window {win_idx}: bad shape {u_nk.shape}")

        if np.isnan(u_nk).any() or np.isinf(u_nk).any():
            raise ValueError(f"Window {win_idx}: NaN/Inf values")

        # Cap extreme values
        u_min = u_nk.min()
        n_capped = np.sum(u_nk > u_min + max_u)
        if n_capped > 0:
            u_nk = np.clip(u_nk, None, u_min + max_u)
            n_capped_total += n_capped

        u_nk_list.append(u_nk)
        n_samples_list.append(u_nk.shape[0])

    if n_capped_total > 0:
        print(f"    Capped {n_capped_total} extreme values")

    return u_nk_list, np.array(n_samples_list)


def run_mbar(sys_dir, name):
    """Run MBAR analysis on a system directory."""
    print(f"\n  Analyzing {name}...")
    print(f"    Directory: {sys_dir}")

    u_nk_list, N_k = load_u_nk_with_cap(sys_dir)

    # Stack and transpose for MBAR
    u_kn_stacked = np.vstack(u_nk_list)
    u_kn = u_kn_stacked.T

    print(f"    Windows: {len(N_k)}, Samples: {N_k.sum()}")

    # Run MBAR
    mbar = MBAR(u_kn, N_k, verbose=False)
    results = mbar.compute_free_energy_differences()

    dG_kT = results['Delta_f'][0, -1]
    dG_err_kT = results['dDelta_f'][0, -1]

    dG_kcal = dG_kT * RT * KCAL_PER_KJ
    dG_err_kcal = dG_err_kT * RT * KCAL_PER_KJ

    print(f"    ΔG = {dG_kcal:+.2f} ± {dG_err_kcal:.2f} kcal/mol")

    # Get overlap
    overlap = mbar.compute_overlap()
    min_overlap = min(overlap['matrix'][i, i+1] for i in range(len(N_k)-1))
    print(f"    Min overlap: {min_overlap:.3f}")

    return dG_kcal, dG_err_kcal, min_overlap


def main():
    base_dir = Path(__file__).parent.resolve()

    print("=" * 70)
    print("ABSOLUTE BINDING FREE ENERGY CALCULATION")
    print("=" * 70)

    results = {}

    # =========================================================================
    # FEBUXOSTAT
    # =========================================================================
    print("\n" + "=" * 70)
    print("FEBUXOSTAT (Drug)")
    print("=" * 70)

    fep_pmx = base_dir / "fep_pmx"

    # Complex WT
    try:
        dG_complex_wt, err_complex_wt, _ = run_mbar(fep_pmx / "wt_complex", "WT Complex")
        results['feb_complex_wt'] = (dG_complex_wt, err_complex_wt)
    except Exception as e:
        print(f"  WT Complex ERROR: {e}")
        results['feb_complex_wt'] = None

    # Complex MUT
    try:
        dG_complex_mut, err_complex_mut, _ = run_mbar(fep_pmx / "mut_complex", "MUT Complex (R702W)")
        results['feb_complex_mut'] = (dG_complex_mut, err_complex_mut)
    except Exception as e:
        print(f"  MUT Complex ERROR: {e}")
        results['feb_complex_mut'] = None

    # Solvent
    try:
        dG_solvent_feb, err_solvent_feb, _ = run_mbar(fep_pmx / "solvent", "Solvent")
        results['feb_solvent'] = (dG_solvent_feb, err_solvent_feb)
    except Exception as e:
        print(f"  Solvent ERROR: {e}")
        results['feb_solvent'] = None

    # =========================================================================
    # NATURAL PRODUCT (CID_10120 Bufadienolide)
    # NOTE: Previously mislabeled as CID_10592/Dihydrocortisol - correct ID is CID_10120
    # =========================================================================
    print("\n" + "=" * 70)
    print("CID_10120 (Bufadienolide - Natural Product)")
    print("=" * 70)

    fep_natural = base_dir / "fep_pmx_natural"

    # Complex WT
    try:
        dG_complex_wt_nat, err_complex_wt_nat, _ = run_mbar(fep_natural / "wt_complex", "WT Complex")
        results['nat_complex_wt'] = (dG_complex_wt_nat, err_complex_wt_nat)
    except Exception as e:
        print(f"  WT Complex ERROR: {e}")
        results['nat_complex_wt'] = None

    # Complex MUT
    try:
        dG_complex_mut_nat, err_complex_mut_nat, _ = run_mbar(fep_natural / "mut_complex", "MUT Complex (R702W)")
        results['nat_complex_mut'] = (dG_complex_mut_nat, err_complex_mut_nat)
    except Exception as e:
        print(f"  MUT Complex ERROR: {e}")
        results['nat_complex_mut'] = None

    # Solvent
    try:
        dG_solvent_nat, err_solvent_nat, _ = run_mbar(fep_natural / "solvent", "Solvent")
        results['nat_solvent'] = (dG_solvent_nat, err_solvent_nat)
    except Exception as e:
        print(f"  Solvent ERROR: {e}")
        results['nat_solvent'] = None

    # =========================================================================
    # COMPUTE ABSOLUTE BINDING FREE ENERGIES
    # =========================================================================
    print("\n" + "=" * 70)
    print("ABSOLUTE BINDING FREE ENERGIES")
    print("=" * 70)
    print("\nΔG_bind = ΔG_complex - ΔG_solvent")
    print("(More negative = stronger binding)\n")

    def compute_binding(complex_result, solvent_result):
        if complex_result is None or solvent_result is None:
            return None, None
        dG_complex, err_complex = complex_result
        dG_solvent, err_solvent = solvent_result
        dG_bind = dG_complex - dG_solvent
        err_bind = np.sqrt(err_complex**2 + err_solvent**2)
        return dG_bind, err_bind

    # Febuxostat
    print("FEBUXOSTAT:")
    if results.get('feb_solvent'):
        print(f"  ΔG_solvent = {results['feb_solvent'][0]:+.2f} ± {results['feb_solvent'][1]:.2f} kcal/mol")

    feb_bind_wt, feb_bind_wt_err = compute_binding(results.get('feb_complex_wt'), results.get('feb_solvent'))
    feb_bind_mut, feb_bind_mut_err = compute_binding(results.get('feb_complex_mut'), results.get('feb_solvent'))

    if feb_bind_wt is not None:
        print(f"  ΔG_bind(WT)  = {feb_bind_wt:+.2f} ± {feb_bind_wt_err:.2f} kcal/mol")
    if feb_bind_mut is not None:
        print(f"  ΔG_bind(MUT) = {feb_bind_mut:+.2f} ± {feb_bind_mut_err:.2f} kcal/mol")
    if feb_bind_wt is not None and feb_bind_mut is not None:
        ddG_feb = feb_bind_mut - feb_bind_wt
        ddG_feb_err = np.sqrt(feb_bind_wt_err**2 + feb_bind_mut_err**2)
        print(f"  ΔΔG (MUT-WT) = {ddG_feb:+.2f} ± {ddG_feb_err:.2f} kcal/mol")

    # Natural Product
    print("\nCID_10120 (Bufadienolide):")
    if results.get('nat_solvent'):
        print(f"  ΔG_solvent = {results['nat_solvent'][0]:+.2f} ± {results['nat_solvent'][1]:.2f} kcal/mol")

    nat_bind_wt, nat_bind_wt_err = compute_binding(results.get('nat_complex_wt'), results.get('nat_solvent'))
    nat_bind_mut, nat_bind_mut_err = compute_binding(results.get('nat_complex_mut'), results.get('nat_solvent'))

    if nat_bind_wt is not None:
        print(f"  ΔG_bind(WT)  = {nat_bind_wt:+.2f} ± {nat_bind_wt_err:.2f} kcal/mol")
    if nat_bind_mut is not None:
        print(f"  ΔG_bind(MUT) = {nat_bind_mut:+.2f} ± {nat_bind_mut_err:.2f} kcal/mol")
    if nat_bind_wt is not None and nat_bind_mut is not None:
        ddG_nat = nat_bind_mut - nat_bind_wt
        ddG_nat_err = np.sqrt(nat_bind_wt_err**2 + nat_bind_mut_err**2)
        print(f"  ΔΔG (MUT-WT) = {ddG_nat:+.2f} ± {ddG_nat_err:.2f} kcal/mol")

    # =========================================================================
    # COMPARISON TABLE
    # =========================================================================
    print("\n" + "=" * 70)
    print("COMPARISON TABLE")
    print("=" * 70)
    print("\n{:<25} {:>15} {:>15} {:>15}".format(
        "Compound / System", "ΔG_complex", "ΔG_solvent", "ΔG_bind"))
    print("-" * 70)

    if results.get('feb_complex_wt') and results.get('feb_solvent'):
        print("{:<25} {:>+12.2f} {:>+12.2f} {:>+12.2f}".format(
            "Febuxostat (WT)",
            results['feb_complex_wt'][0],
            results['feb_solvent'][0],
            feb_bind_wt))

    if results.get('feb_complex_mut') and results.get('feb_solvent'):
        print("{:<25} {:>+12.2f} {:>+12.2f} {:>+12.2f}".format(
            "Febuxostat (R702W)",
            results['feb_complex_mut'][0],
            results['feb_solvent'][0],
            feb_bind_mut))

    if results.get('nat_complex_wt') and results.get('nat_solvent'):
        print("{:<25} {:>+12.2f} {:>+12.2f} {:>+12.2f}".format(
            "CID_10120 (WT)",
            results['nat_complex_wt'][0],
            results['nat_solvent'][0],
            nat_bind_wt))

    if results.get('nat_complex_mut') and results.get('nat_solvent'):
        print("{:<25} {:>+12.2f} {:>+12.2f} {:>+12.2f}".format(
            "CID_10120 (R702W)",
            results['nat_complex_mut'][0],
            results['nat_solvent'][0],
            nat_bind_mut))

    print("-" * 70)
    print("(All values in kcal/mol; more negative ΔG_bind = stronger binding)")

    # =========================================================================
    # SAVE RESULTS
    # =========================================================================
    output = {
        'febuxostat': {
            'dG_complex_wt': results.get('feb_complex_wt'),
            'dG_complex_mut': results.get('feb_complex_mut'),
            'dG_solvent': results.get('feb_solvent'),
            'dG_bind_wt': (feb_bind_wt, feb_bind_wt_err) if feb_bind_wt else None,
            'dG_bind_mut': (feb_bind_mut, feb_bind_mut_err) if feb_bind_mut else None,
        },
        'cid_10120': {  # NOTE: Previously mislabeled as cid_10592
            'dG_complex_wt': results.get('nat_complex_wt'),
            'dG_complex_mut': results.get('nat_complex_mut'),
            'dG_solvent': results.get('nat_solvent'),
            'dG_bind_wt': (nat_bind_wt, nat_bind_wt_err) if nat_bind_wt else None,
            'dG_bind_mut': (nat_bind_mut, nat_bind_mut_err) if nat_bind_mut else None,
        }
    }

    # Convert tuples to lists for JSON
    def convert_for_json(obj):
        if isinstance(obj, dict):
            return {k: convert_for_json(v) for k, v in obj.items()}
        elif isinstance(obj, tuple):
            return list(obj)
        elif isinstance(obj, np.ndarray):
            return obj.tolist()
        return obj

    output_file = base_dir / "absolute_binding_results.json"
    with open(output_file, 'w') as f:
        json.dump(convert_for_json(output), f, indent=2)
    print(f"\nResults saved to: {output_file}")


if __name__ == "__main__":
    main()
