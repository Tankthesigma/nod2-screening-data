#!/usr/bin/env python
"""
Comprehensive MBAR analysis for FEP verification
Checks D, E, F: Per-phase breakdown, window consistency, bootstrap convergence
"""
import numpy as np
import os
from pymbar import MBAR

BASE = "C:/Users/vasud/nod2-screening-data/fep_complete/fep_pmx"

# Constants
TEMPERATURE = 310.0  # K
KB = 0.001987204  # kcal/(mol·K)
BETA = 1.0 / (KB * TEMPERATURE)
KJ_TO_KCAL = 1.0 / 4.184

def load_lambda_schedule(system_dir):
    """Load lambda schedule and determine phases."""
    sched = np.load(os.path.join(system_dir, "lambda_schedule.npy"))
    # Columns: electrostatics, sterics, restraints
    return sched

def load_u_nk_data(system_dir, window_range=None):
    """Load u_nk data from specified windows."""
    if window_range is None:
        window_range = range(20)

    u_nk_list = []
    N_k = []

    for w in window_range:
        u_nk_file = os.path.join(system_dir, f"window_{w:02d}", "u_nk.npy")
        if os.path.exists(u_nk_file):
            u_nk = np.load(u_nk_file)
            # Convert kJ/mol to reduced units (dimensionless)
            u_nk_reduced = u_nk * KJ_TO_KCAL * BETA
            u_nk_list.append(u_nk_reduced)
            N_k.append(len(u_nk))
        else:
            print(f"  Warning: {u_nk_file} not found")

    return u_nk_list, np.array(N_k)

def build_u_kn_matrix(u_nk_list, N_k, state_range=None):
    """Build the u_kn matrix for MBAR."""
    if state_range is None:
        K = u_nk_list[0].shape[1]
        state_range = range(K)

    K = len(state_range)
    N_total = sum(N_k)

    u_kn = np.zeros((K, N_total))

    n_start = 0
    for w, u_nk in enumerate(u_nk_list):
        n_samples = N_k[w]
        for k_idx, k in enumerate(state_range):
            u_kn[k_idx, n_start:n_start+n_samples] = u_nk[:, k]
        n_start += n_samples

    return u_kn

def run_mbar_subset(system_dir, window_range, label):
    """Run MBAR on a subset of windows."""
    u_nk_list, N_k = load_u_nk_data(system_dir, window_range)

    if len(u_nk_list) == 0:
        return None, None

    # Build matrix with only the relevant states
    state_range = list(window_range)
    u_kn = build_u_kn_matrix(u_nk_list, N_k, state_range)

    # Adjust N_k for the subset
    N_k_subset = np.zeros(len(state_range), dtype=int)
    for i, w in enumerate(window_range):
        idx = list(window_range).index(w)
        if idx < len(N_k):
            N_k_subset[i] = N_k[idx]

    try:
        mbar = MBAR(u_kn, N_k_subset, verbose=False)
        results = mbar.compute_free_energy_differences()
        dG = results['Delta_f'][0, -1] / BETA  # kcal/mol
        dG_err = results['dDelta_f'][0, -1] / BETA
        return dG, dG_err
    except Exception as e:
        print(f"  MBAR error for {label}: {e}")
        return None, None

def run_full_mbar(system_dir):
    """Run MBAR on all 20 windows."""
    u_nk_list, N_k = load_u_nk_data(system_dir)

    K = 20
    N_total = sum(N_k)
    u_kn = np.zeros((K, N_total))

    n_start = 0
    for w, u_nk in enumerate(u_nk_list):
        n_samples = N_k[w]
        for k in range(K):
            u_kn[k, n_start:n_start+n_samples] = u_nk[:, k]
        n_start += n_samples

    mbar = MBAR(u_kn, N_k, verbose=False)
    results = mbar.compute_free_energy_differences()

    # Get cumulative dG for each state
    dG_cumulative = []
    for k in range(K):
        dG_cumulative.append(results['Delta_f'][0, k] / BETA)

    dG_total = results['Delta_f'][0, -1] / BETA
    dG_err = results['dDelta_f'][0, -1] / BETA

    return dG_total, dG_err, dG_cumulative, u_nk_list, N_k

def bootstrap_mbar(u_nk_list, N_k, n_bootstrap=100):
    """Bootstrap resampling of MBAR results."""
    np.random.seed(42)
    bootstrap_dG = []

    K = 20

    for b in range(n_bootstrap):
        # Resample frames within each window
        u_nk_boot = []
        for w, u_nk in enumerate(u_nk_list):
            n_samples = len(u_nk)
            indices = np.random.choice(n_samples, size=n_samples, replace=True)
            u_nk_boot.append(u_nk[indices, :])

        # Build u_kn matrix
        N_total = sum(N_k)
        u_kn = np.zeros((K, N_total))
        n_start = 0
        for w, u_nk in enumerate(u_nk_boot):
            n_samples = N_k[w]
            for k in range(K):
                u_kn[k, n_start:n_start+n_samples] = u_nk[:, k]
            n_start += n_samples

        try:
            mbar = MBAR(u_kn, N_k, verbose=False)
            results = mbar.compute_free_energy_differences()
            dG = results['Delta_f'][0, -1] / BETA
            bootstrap_dG.append(dG)
        except:
            pass

    bootstrap_dG = np.array(bootstrap_dG)
    return {
        'mean': np.mean(bootstrap_dG),
        'std': np.std(bootstrap_dG),
        'ci_lower': np.percentile(bootstrap_dG, 2.5),
        'ci_upper': np.percentile(bootstrap_dG, 97.5),
        'n_success': len(bootstrap_dG)
    }

def main():
    print("=" * 70)
    print("COMPREHENSIVE MBAR ANALYSIS - FEP VERIFICATION")
    print("=" * 70)

    results = {}

    for sys_name in ['solvent', 'wt_complex', 'mut_complex']:
        print(f"\n{'='*70}")
        print(f"SYSTEM: {sys_name}")
        print(f"{'='*70}")

        system_dir = f"{BASE}/{sys_name}"

        # Load lambda schedule
        sched = load_lambda_schedule(system_dir)
        print(f"\nLambda schedule:")
        print(f"  Window 0:  elec={sched[0,0]:.2f}, sterics={sched[0,1]:.2f}, restraints={sched[0,2]:.2f}")
        print(f"  Window 19: elec={sched[19,0]:.2f}, sterics={sched[19,1]:.2f}, restraints={sched[19,2]:.2f}")

        # Check D: Full MBAR and per-state breakdown
        print(f"\n--- CHECK D: FULL MBAR & PER-STATE BREAKDOWN ---")
        dG, dG_err, dG_cumulative, u_nk_list, N_k = run_full_mbar(system_dir)
        print(f"  Total dG (0→19): {dG:.3f} ± {dG_err:.3f} kcal/mol")
        print(f"\n  Cumulative dG by window:")
        print(f"  {'Window':<8} {'Cumulative dG':<15} {'Increment'}")
        print(f"  " + "-" * 40)
        prev_dG = 0
        for w in range(20):
            increment = dG_cumulative[w] - prev_dG
            print(f"  {w:<8} {dG_cumulative[w]:>12.3f}    {increment:>+8.3f}")
            prev_dG = dG_cumulative[w]

        results[sys_name] = {
            'dG': dG,
            'dG_err': dG_err,
            'dG_cumulative': dG_cumulative
        }

        # Check E: Windows 0-14 vs 15-19
        print(f"\n--- CHECK E: WINDOW SUBSET ANALYSIS ---")
        dG_0_14, err_0_14 = run_mbar_subset(system_dir, range(15), "windows 0-14")
        dG_15_19, err_15_19 = run_mbar_subset(system_dir, range(15, 20), "windows 15-19")

        if dG_0_14 is not None:
            print(f"  Windows 0-14:  dG = {dG_0_14:.3f} ± {err_0_14:.3f} kcal/mol")
        if dG_15_19 is not None:
            print(f"  Windows 15-19: dG = {dG_15_19:.3f} ± {err_15_19:.3f} kcal/mol")

        # Also check 0-9 and 10-19
        dG_0_9, err_0_9 = run_mbar_subset(system_dir, range(10), "windows 0-9")
        dG_10_19, err_10_19 = run_mbar_subset(system_dir, range(10, 20), "windows 10-19")

        if dG_0_9 is not None:
            print(f"  Windows 0-9:   dG = {dG_0_9:.3f} ± {err_0_9:.3f} kcal/mol")
        if dG_10_19 is not None:
            print(f"  Windows 10-19: dG = {dG_10_19:.3f} ± {err_10_19:.3f} kcal/mol")

        # Check F: Bootstrap
        print(f"\n--- CHECK F: BOOTSTRAP ANALYSIS (100 resamples) ---")
        boot = bootstrap_mbar(u_nk_list, N_k, n_bootstrap=100)
        print(f"  Bootstrap mean:   {boot['mean']:.3f} kcal/mol")
        print(f"  Bootstrap std:    {boot['std']:.3f} kcal/mol")
        print(f"  95% CI:           [{boot['ci_lower']:.3f}, {boot['ci_upper']:.3f}]")
        print(f"  MBAR error:       {dG_err:.3f} kcal/mol")
        print(f"  Bootstrap/MBAR ratio: {boot['std']/dG_err:.2f}x")

        results[sys_name]['bootstrap'] = boot

    # Final summary: DDG calculation
    print(f"\n{'='*70}")
    print("FINAL DDG CALCULATION")
    print(f"{'='*70}")

    dG_solv = results['solvent']['dG']
    dG_wt = results['wt_complex']['dG']
    dG_mut = results['mut_complex']['dG']

    err_solv = results['solvent']['dG_err']
    err_wt = results['wt_complex']['dG_err']
    err_mut = results['mut_complex']['dG_err']

    # Boresch corrections (from stored data)
    boresch_wt = 7.588
    boresch_mut = 7.588

    dG_bind_wt = dG_solv - dG_wt + boresch_wt
    dG_bind_mut = dG_solv - dG_mut + boresch_mut
    ddG = dG_bind_mut - dG_bind_wt

    # Error propagation
    err_bind_wt = np.sqrt(err_solv**2 + err_wt**2)
    err_bind_mut = np.sqrt(err_solv**2 + err_mut**2)
    err_ddG = np.sqrt(err_wt**2 + err_mut**2)  # Solvent cancels

    print(f"\nMBAR Results:")
    print(f"  dG_solvent:     {dG_solv:.3f} ± {err_solv:.3f} kcal/mol")
    print(f"  dG_wt_complex:  {dG_wt:.3f} ± {err_wt:.3f} kcal/mol")
    print(f"  dG_mut_complex: {dG_mut:.3f} ± {err_mut:.3f} kcal/mol")

    print(f"\nBinding Free Energies:")
    print(f"  dG_bind(WT):  {dG_bind_wt:.3f} ± {err_bind_wt:.3f} kcal/mol")
    print(f"  dG_bind(Mut): {dG_bind_mut:.3f} ± {err_bind_mut:.3f} kcal/mol")

    print(f"\n  DDG_bind = {ddG:.3f} ± {err_ddG:.3f} kcal/mol")
    print(f"  Significance: {abs(ddG)/err_ddG:.1f} sigma from zero")

    # Bootstrap-based DDG
    boot_wt = results['wt_complex']['bootstrap']
    boot_mut = results['mut_complex']['bootstrap']
    boot_ddG_std = np.sqrt(boot_wt['std']**2 + boot_mut['std']**2)

    print(f"\nBootstrap-based error:")
    print(f"  DDG std (bootstrap): {boot_ddG_std:.3f} kcal/mol")
    print(f"  DDG std (MBAR):      {err_ddG:.3f} kcal/mol")
    print(f"  Ratio:               {boot_ddG_std/err_ddG:.2f}x")

    print(f"\n{'='*70}")
    print("VERIFICATION SUMMARY")
    print(f"{'='*70}")
    print(f"  [{'PASS' if abs(ddG)/err_ddG > 3 else 'WARN'}] DDG significantly different from zero: {abs(ddG)/err_ddG:.1f}σ")
    print(f"  [{'PASS' if boot_ddG_std/err_ddG < 2 else 'WARN'}] Bootstrap error consistent with MBAR: {boot_ddG_std/err_ddG:.2f}x")

    # Check for outlier windows
    print(f"\n  Checking for outlier windows...")
    for sys_name in ['wt_complex', 'mut_complex']:
        cumul = results[sys_name]['dG_cumulative']
        increments = [cumul[0]] + [cumul[i] - cumul[i-1] for i in range(1, 20)]
        mean_incr = np.mean(np.abs(increments))
        std_incr = np.std(increments)
        outliers = [i for i, inc in enumerate(increments) if abs(inc) > mean_incr + 3*std_incr]
        if outliers:
            print(f"  [{sys_name}] Potential outlier windows: {outliers}")
        else:
            print(f"  [{sys_name}] No outlier windows detected")

if __name__ == "__main__":
    main()
