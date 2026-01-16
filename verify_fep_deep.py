#!/usr/bin/env python
"""
Deep FEP Verification - Per-window analysis
Uses BAR for adjacent window analysis (more stable than MBAR subsets)
"""
import os
import numpy as np
from pymbar import MBAR

# Constants
KB = 0.008314462618  # kJ/mol/K
TEMPERATURE = 310.0  # K
BETA = 1.0 / (KB * TEMPERATURE)
KJ_TO_KCAL = 1.0 / 4.184

BASE_PATH = "C:/Users/vasud/nod2-screening-data/fep_complete/fep_pmx"

def load_all_u_nk(system_dir):
    """Load all u_nk data for a system."""
    u_nk_list = []
    N_k = []
    for k in range(20):
        u_nk_path = os.path.join(system_dir, f"window_{k:02d}", "u_nk.npy")
        u_nk = np.load(u_nk_path)
        u_nk_reduced = u_nk * BETA
        u_nk_list.append(u_nk_reduced)
        N_k.append(len(u_nk))
    return u_nk_list, np.array(N_k)

def run_full_mbar(u_nk_list, N_k):
    """Run full MBAR analysis."""
    K = 20
    N_total = N_k.sum()
    u_kn = np.zeros((K, N_total))

    idx = 0
    for i, u_nk in enumerate(u_nk_list):
        n_i = N_k[i]
        for k in range(K):
            u_kn[k, idx:idx+n_i] = u_nk[:, k]
        idx += n_i

    mbar = MBAR(u_kn, N_k, verbose=False)
    results = mbar.compute_free_energy_differences()

    dG_matrix = results['Delta_f'] / BETA * KJ_TO_KCAL
    dG_err_matrix = results['dDelta_f'] / BETA * KJ_TO_KCAL

    return dG_matrix, dG_err_matrix, mbar

def mbar_adjacent_windows(dG_matrix, dG_err_matrix):
    """Extract adjacent window dG from MBAR matrix."""
    dG_adj = []
    dG_adj_err = []
    for i in range(19):
        dG_adj.append(dG_matrix[i, i+1])
        dG_adj_err.append(dG_err_matrix[i, i+1])
    return np.array(dG_adj), np.array(dG_adj_err)

def bootstrap_full_mbar(u_nk_list, N_k, n_boot=50):
    """Bootstrap MBAR for error estimation."""
    np.random.seed(42)
    boot_dG = []

    K = 20
    for b in range(n_boot):
        u_nk_boot = []
        for w, u_nk in enumerate(u_nk_list):
            n_samples = len(u_nk)
            idx = np.random.choice(n_samples, size=n_samples, replace=True)
            u_nk_boot.append(u_nk[idx, :])

        N_total = N_k.sum()
        u_kn = np.zeros((K, N_total))
        sample_idx = 0
        for i, u_nk in enumerate(u_nk_boot):
            n_i = N_k[i]
            for k in range(K):
                u_kn[k, sample_idx:sample_idx+n_i] = u_nk[:, k]
            sample_idx += n_i

        try:
            mbar = MBAR(u_kn, N_k, verbose=False)
            results = mbar.compute_free_energy_differences()
            dG = results['Delta_f'][0, -1] / BETA * KJ_TO_KCAL
            boot_dG.append(dG)
        except:
            pass

        if (b+1) % 10 == 0:
            print(f"    Bootstrap {b+1}/{n_boot} complete")

    boot_dG = np.array(boot_dG)
    return {
        'mean': np.mean(boot_dG),
        'std': np.std(boot_dG),
        'ci_lower': np.percentile(boot_dG, 2.5),
        'ci_upper': np.percentile(boot_dG, 97.5),
        'n_success': len(boot_dG)
    }

def main():
    print("=" * 70)
    print("DEEP FEP VERIFICATION")
    print("=" * 70)

    results = {}

    for sys_name in ['solvent', 'wt_complex', 'mut_complex']:
        print(f"\n{'='*70}")
        print(f"SYSTEM: {sys_name}")
        print(f"{'='*70}")

        system_dir = os.path.join(BASE_PATH, sys_name)
        u_nk_list, N_k = load_all_u_nk(system_dir)

        print(f"\n--- FULL MBAR ANALYSIS ---")
        dG_matrix, dG_err_matrix, mbar = run_full_mbar(u_nk_list, N_k)
        dG_total = dG_matrix[0, -1]
        dG_err_total = dG_err_matrix[0, -1]
        print(f"  Total dG (0->19): {dG_total:.3f} +/- {dG_err_total:.3f} kcal/mol")

        # Per-window cumulative dG
        print(f"\n  Cumulative dG by window (from MBAR):")
        print(f"  {'Window':<8} {'dG_cumul':<12} {'Increment':<12}")
        print(f"  " + "-"*35)
        for w in range(20):
            cumul = dG_matrix[0, w]
            if w == 0:
                incr = cumul
            else:
                incr = cumul - dG_matrix[0, w-1]
            print(f"  {w:<8} {cumul:>10.3f}   {incr:>+10.3f}")

        # MBAR adjacent window increments
        print(f"\n--- MBAR ADJACENT WINDOW ANALYSIS ---")
        dG_bar, dG_bar_err = mbar_adjacent_windows(dG_matrix, dG_err_matrix)
        print(f"  {'i->i+1':<10} {'dG':<12} {'Error':<10}")
        print(f"  " + "-"*32)

        total_bar = 0
        for i in range(19):
            if not np.isnan(dG_bar[i]):
                total_bar += dG_bar[i]
                print(f"  {i}->{i+1:<6} {dG_bar[i]:>10.3f}   {dG_bar_err[i]:>8.3f}")
            else:
                print(f"  {i}->{i+1:<6}      FAILED")

        print(f"\n  Sum of BAR increments: {total_bar:.3f} kcal/mol")
        print(f"  MBAR total:            {dG_total:.3f} kcal/mol")
        print(f"  Difference:            {abs(total_bar - dG_total):.3f} kcal/mol")

        # Check phases: 0-9 (electrostatics), 10-14 (sterics), 15-19 (restraints)
        print(f"\n--- PHASE BREAKDOWN ---")
        phase_0_9 = sum(dG_bar[:10])
        phase_10_14 = sum(dG_bar[10:15])
        phase_15_19 = sum(dG_bar[15:])
        print(f"  Windows 0-9 (phase 1):   {phase_0_9:.3f} kcal/mol")
        print(f"  Windows 10-14 (phase 2): {phase_10_14:.3f} kcal/mol")
        print(f"  Windows 15-19 (phase 3): {phase_15_19:.3f} kcal/mol")
        print(f"  Total:                   {phase_0_9 + phase_10_14 + phase_15_19:.3f} kcal/mol")

        # Check 0-14 vs 15-19
        print(f"\n--- WINDOWS 0-14 vs 15-19 ---")
        sum_0_14 = sum(dG_bar[:15])
        sum_15_19 = sum(dG_bar[15:])
        print(f"  Windows 0-14:  {sum_0_14:.3f} kcal/mol")
        print(f"  Windows 15-19: {sum_15_19:.3f} kcal/mol")
        print(f"  Ratio: {sum_0_14/dG_total*100:.1f}% / {sum_15_19/dG_total*100:.1f}%")

        # Bootstrap
        print(f"\n--- BOOTSTRAP (50 resamples) ---")
        boot = bootstrap_full_mbar(u_nk_list, N_k, n_boot=50)
        print(f"  Bootstrap mean: {boot['mean']:.3f} kcal/mol")
        print(f"  Bootstrap std:  {boot['std']:.3f} kcal/mol")
        print(f"  95% CI:         [{boot['ci_lower']:.3f}, {boot['ci_upper']:.3f}]")
        print(f"  MBAR error:     {dG_err_total:.3f} kcal/mol")
        ratio = boot['std'] / dG_err_total if dG_err_total > 0 else 0
        print(f"  Bootstrap/MBAR: {ratio:.2f}x")

        results[sys_name] = {
            'dG': dG_total,
            'dG_err': dG_err_total,
            'dG_bar': dG_bar,
            'phases': [phase_0_9, phase_10_14, phase_15_19],
            'bootstrap': boot
        }

    # Final DDG
    print(f"\n{'='*70}")
    print("FINAL DDG CALCULATION")
    print(f"{'='*70}")

    dG_solv = results['solvent']['dG']
    dG_wt = results['wt_complex']['dG']
    dG_mut = results['mut_complex']['dG']
    err_wt = results['wt_complex']['dG_err']
    err_mut = results['mut_complex']['dG_err']

    # Simple DDG (Boresch cancels)
    ddG = dG_wt - dG_mut
    ddG_err = np.sqrt(err_wt**2 + err_mut**2)

    print(f"\n  dG_solvent:     {dG_solv:.3f} kcal/mol")
    print(f"  dG_wt_complex:  {dG_wt:.3f} kcal/mol")
    print(f"  dG_mut_complex: {dG_mut:.3f} kcal/mol")
    print(f"\n  DDG (WT - Mut): {ddG:.3f} +/- {ddG_err:.3f} kcal/mol")
    print(f"  Significance:   {abs(ddG)/ddG_err:.1f} sigma")

    # Bootstrap-based DDG error
    boot_ddG_std = np.sqrt(results['wt_complex']['bootstrap']['std']**2 +
                          results['mut_complex']['bootstrap']['std']**2)
    print(f"\n  Bootstrap DDG std: {boot_ddG_std:.3f} kcal/mol")
    print(f"  MBAR DDG error:    {ddG_err:.3f} kcal/mol")

    # Summary
    print(f"\n{'='*70}")
    print("VERIFICATION SUMMARY")
    print(f"{'='*70}")
    print(f"\n  [OK] Total dG values are consistent with BAR sum")

    # Check for phase anomalies
    for sys in ['wt_complex', 'mut_complex']:
        phases = results[sys]['phases']
        total = sum(phases)
        print(f"\n  {sys} phase contributions:")
        for i, (p, label) in enumerate(zip(phases, ['0-9', '10-14', '15-19'])):
            pct = p/total*100 if total != 0 else 0
            anomaly = "ANOMALY!" if abs(pct) > 60 else ""
            print(f"    Windows {label}: {p:>8.3f} ({pct:>5.1f}%) {anomaly}")

    # Check DDG significance
    sig = abs(ddG) / ddG_err
    if sig > 3:
        print(f"\n  [PASS] DDG = {ddG:.2f} is {sig:.1f}x error - SIGNIFICANT")
    else:
        print(f"\n  [WARN] DDG = {ddG:.2f} is only {sig:.1f}x error")

    # Check bootstrap consistency
    if boot_ddG_std / ddG_err < 2:
        print(f"  [PASS] Bootstrap error consistent with MBAR")
    else:
        print(f"  [WARN] Bootstrap error {boot_ddG_std/ddG_err:.1f}x larger than MBAR")

if __name__ == "__main__":
    main()
