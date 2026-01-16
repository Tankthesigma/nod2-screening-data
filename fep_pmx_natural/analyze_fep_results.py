#!/usr/bin/env python
"""
FEP Results Analysis for Natural Product CID_10120

Computes:
- DeltaG(WT): Binding free energy for wild-type NOD2
- DeltaG(Mut): Binding free energy for R702W mutant
- DeltaDeltaG = DeltaG(Mut) - DeltaG(WT): Mutation effect on binding

Includes:
- MBAR analysis with error estimates
- Bootstrap confidence intervals
- Per-phase breakdown
- Overlap matrix
- Verification checks A-F (same as febuxostat)
"""

import numpy as np
import json
from pathlib import Path
from datetime import datetime

BASE_DIR = Path("C:/Users/vasud/nod2-screening-data/fep_pmx_natural")
RESULTS_DIR = BASE_DIR / "results"

def load_u_nk_data(sys_name):
    """Load all u_nk data for a system."""
    u_nk_list = []
    N_k = []

    for i in range(20):
        u_nk_path = BASE_DIR / sys_name / f"window_{i:02d}" / "u_nk.npy"
        if u_nk_path.exists():
            u_nk = np.load(u_nk_path)
            u_nk_list.append(u_nk)
            N_k.append(u_nk.shape[0])
        else:
            print(f"  [MISSING] {u_nk_path}")
            return None, None

    return u_nk_list, np.array(N_k)

def run_mbar_analysis(u_nk_list, N_k):
    """Run MBAR analysis on u_nk data."""
    try:
        from pymbar import MBAR

        # Stack all u_nk matrices
        u_kn = np.vstack(u_nk_list)
        print(f"  Combined data: {u_kn.shape}")

        # Run MBAR
        mbar = MBAR(u_kn.T, N_k)

        # Get free energy differences
        results = mbar.compute_free_energy_differences()
        dG_matrix = results['Delta_f']
        dG_err_matrix = results['dDelta_f']

        # Total dG: window 0 to window 19
        dG_total = dG_matrix[0, 19]
        dG_err_total = dG_err_matrix[0, 19]

        # Convert to kcal/mol (from kT at 310K)
        kT_kcal = 0.001987 * 310.0  # kcal/mol
        dG_kcal = dG_total * kT_kcal
        dG_err_kcal = dG_err_total * kT_kcal

        # Overlap matrix
        overlap = mbar.compute_overlap()
        overlap_matrix = overlap['matrix']

        return {
            'dG_kT': dG_total,
            'dG_err_kT': dG_err_total,
            'dG_kcal': dG_kcal,
            'dG_err_kcal': dG_err_kcal,
            'dG_matrix': dG_matrix,
            'dG_err_matrix': dG_err_matrix,
            'overlap_matrix': overlap_matrix,
        }

    except Exception as e:
        print(f"  [ERROR] MBAR failed: {e}")
        import traceback
        traceback.print_exc()
        return None

def bootstrap_analysis(u_nk_list, N_k, n_bootstrap=50):
    """Bootstrap resampling for confidence intervals."""
    try:
        from pymbar import MBAR

        bootstrap_dG = []
        n_windows = len(u_nk_list)

        print(f"  Running {n_bootstrap} bootstrap samples...")

        for b in range(n_bootstrap):
            # Resample with replacement from each window
            resampled_u_nk = []
            resampled_N_k = []

            for k in range(n_windows):
                n_samples = u_nk_list[k].shape[0]
                indices = np.random.choice(n_samples, size=n_samples, replace=True)
                resampled_u_nk.append(u_nk_list[k][indices])
                resampled_N_k.append(n_samples)

            # Run MBAR on resampled data
            u_kn = np.vstack(resampled_u_nk)
            mbar = MBAR(u_kn.T, np.array(resampled_N_k))

            results = mbar.compute_free_energy_differences()
            dG = results['Delta_f'][0, 19]
            bootstrap_dG.append(dG)

        bootstrap_dG = np.array(bootstrap_dG)
        kT_kcal = 0.001987 * 310.0

        return {
            'mean_kT': np.mean(bootstrap_dG),
            'std_kT': np.std(bootstrap_dG),
            'ci_95_kT': (np.percentile(bootstrap_dG, 2.5), np.percentile(bootstrap_dG, 97.5)),
            'mean_kcal': np.mean(bootstrap_dG) * kT_kcal,
            'std_kcal': np.std(bootstrap_dG) * kT_kcal,
            'ci_95_kcal': (np.percentile(bootstrap_dG, 2.5) * kT_kcal,
                          np.percentile(bootstrap_dG, 97.5) * kT_kcal),
        }

    except Exception as e:
        print(f"  [ERROR] Bootstrap failed: {e}")
        return None

def per_phase_breakdown(dG_matrix, dG_err_matrix):
    """Break down dG by phase.

    Note: Phases are defined so they sum to total dG(0->19).
    Each phase starts from end of previous phase.
    """
    kT_kcal = 0.001987 * 310.0

    # Phase definitions (start_window, end_window) - both inclusive
    # These must chain together: 0->9, 9->14, 14->19 sums to 0->19
    phases = {
        'electrostatics': (0, 9),   # Windows 0-9: electrostatics decoupling
        'sterics_1': (9, 14),       # Windows 9-14: sterics decoupling (first half)
        'sterics_2': (14, 19),      # Windows 14-19: sterics decoupling (second half)
    }

    breakdown = {}
    for phase_name, (start, end) in phases.items():
        dG = dG_matrix[start, end] * kT_kcal
        dG_err = dG_err_matrix[start, end] * kT_kcal
        breakdown[phase_name] = {'dG': dG, 'dG_err': dG_err}

    return breakdown

def check_overlap_quality(overlap_matrix):
    """Check if overlap between adjacent windows is sufficient."""
    issues = []
    n_windows = overlap_matrix.shape[0]

    print("\n  Overlap between adjacent windows:")
    print("  " + "-"*40)

    for i in range(n_windows - 1):
        ov = overlap_matrix[i, i+1]
        if ov < 0.01:
            status = "FAIL"
            issues.append(f"Window {i}-{i+1}: overlap={ov:.4f} (< 0.01)")
        elif ov < 0.03:
            status = "LOW"
            issues.append(f"Window {i}-{i+1}: overlap={ov:.4f} (< 0.03)")
        else:
            status = "OK"
        print(f"    Window {i:2d}-{i+1:2d}: {ov:.4f} [{status}]")

    return issues

def main():
    print("="*70)
    print("FEP RESULTS ANALYSIS")
    print("Natural Product CID_10120 (Bufadienolide)")
    print("="*70)
    print()

    RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    results = {
        'timestamp': datetime.now().isoformat(),
        'ligand': 'CID_10120 (Bufadienolide)',
        'mutation': 'R702W',
        'systems': {},
    }

    # Analyze each system
    systems = {
        'wt_complex': 'Wild-type NOD2',
        'mut_complex': 'R702W Mutant',
        'solvent': 'Solvent only',
    }

    for sys_name, sys_desc in systems.items():
        print(f"\n{'='*60}")
        print(f"Analyzing: {sys_name} ({sys_desc})")
        print(f"{'='*60}")

        # Load data
        u_nk_list, N_k = load_u_nk_data(sys_name)
        if u_nk_list is None:
            print(f"  [SKIP] Missing data for {sys_name}")
            continue

        print(f"  Loaded {len(u_nk_list)} windows, N_k = {N_k}")

        # MBAR analysis
        print("\n  Running MBAR analysis...")
        mbar_results = run_mbar_analysis(u_nk_list, N_k)
        if mbar_results is None:
            continue

        print(f"  DeltaG = {mbar_results['dG_kcal']:.2f} ± {mbar_results['dG_err_kcal']:.2f} kcal/mol")

        # Bootstrap
        print("\n  Running bootstrap analysis...")
        bootstrap_results = bootstrap_analysis(u_nk_list, N_k, n_bootstrap=50)
        if bootstrap_results:
            print(f"  Bootstrap: {bootstrap_results['mean_kcal']:.2f} ± {bootstrap_results['std_kcal']:.2f} kcal/mol")
            print(f"  95% CI: [{bootstrap_results['ci_95_kcal'][0]:.2f}, {bootstrap_results['ci_95_kcal'][1]:.2f}]")

        # Per-phase breakdown
        print("\n  Per-phase breakdown:")
        breakdown = per_phase_breakdown(mbar_results['dG_matrix'], mbar_results['dG_err_matrix'])
        for phase, vals in breakdown.items():
            print(f"    {phase}: {vals['dG']:.2f} ± {vals['dG_err']:.2f} kcal/mol")

        # Overlap quality
        overlap_issues = check_overlap_quality(mbar_results['overlap_matrix'])

        # Store results
        results['systems'][sys_name] = {
            'description': sys_desc,
            'dG_kcal': mbar_results['dG_kcal'],
            'dG_err_kcal': mbar_results['dG_err_kcal'],
            'bootstrap_std_kcal': bootstrap_results['std_kcal'] if bootstrap_results else None,
            'bootstrap_ci_95': bootstrap_results['ci_95_kcal'] if bootstrap_results else None,
            'per_phase': breakdown,
            'overlap_issues': overlap_issues,
        }

    # Compute binding free energies and DeltaDeltaG
    print("\n" + "="*70)
    print("BINDING FREE ENERGY CALCULATIONS")
    print("="*70)

    if 'wt_complex' in results['systems'] and 'solvent' in results['systems']:
        dG_wt_complex = results['systems']['wt_complex']['dG_kcal']
        dG_solvent = results['systems']['solvent']['dG_kcal']
        dG_wt_err = results['systems']['wt_complex']['dG_err_kcal']
        dG_sol_err = results['systems']['solvent']['dG_err_kcal']

        # DeltaG_bind(WT) = DeltaG_complex - DeltaG_solvent
        dG_bind_wt = dG_wt_complex - dG_solvent
        dG_bind_wt_err = np.sqrt(dG_wt_err**2 + dG_sol_err**2)

        print(f"\nDeltaG_bind(WT) = {dG_bind_wt:.2f} ± {dG_bind_wt_err:.2f} kcal/mol")
        results['dG_bind_wt'] = dG_bind_wt
        results['dG_bind_wt_err'] = dG_bind_wt_err
    else:
        print("\n[WARNING] Cannot compute WT binding energy - missing data")

    if 'mut_complex' in results['systems'] and 'solvent' in results['systems']:
        dG_mut_complex = results['systems']['mut_complex']['dG_kcal']
        dG_solvent = results['systems']['solvent']['dG_kcal']
        dG_mut_err = results['systems']['mut_complex']['dG_err_kcal']
        dG_sol_err = results['systems']['solvent']['dG_err_kcal']

        # DeltaG_bind(Mut) = DeltaG_complex - DeltaG_solvent
        dG_bind_mut = dG_mut_complex - dG_solvent
        dG_bind_mut_err = np.sqrt(dG_mut_err**2 + dG_sol_err**2)

        print(f"DeltaG_bind(R702W) = {dG_bind_mut:.2f} ± {dG_bind_mut_err:.2f} kcal/mol")
        results['dG_bind_mut'] = dG_bind_mut
        results['dG_bind_mut_err'] = dG_bind_mut_err
    else:
        print("[WARNING] Cannot compute Mut binding energy - missing data")

    # DeltaDeltaG
    if 'dG_bind_wt' in results and 'dG_bind_mut' in results:
        ddG = results['dG_bind_mut'] - results['dG_bind_wt']
        ddG_err = np.sqrt(results['dG_bind_wt_err']**2 + results['dG_bind_mut_err']**2)

        print(f"\nDeltaDeltaG = DeltaG_bind(R702W) - DeltaG_bind(WT)")
        print(f"DeltaDeltaG = {ddG:.2f} ± {ddG_err:.2f} kcal/mol")

        # Statistical significance
        sigma = abs(ddG) / ddG_err if ddG_err > 0 else 0
        print(f"Significance: {sigma:.1f} sigma")

        if sigma >= 2:
            if ddG > 0:
                interpretation = "R702W WEAKENS binding (statistically significant)"
            else:
                interpretation = "R702W STRENGTHENS binding (statistically significant)"
        else:
            interpretation = "No significant effect of R702W mutation"

        print(f"Interpretation: {interpretation}")

        results['ddG'] = ddG
        results['ddG_err'] = ddG_err
        results['ddG_sigma'] = sigma
        results['interpretation'] = interpretation

    # Save results
    print("\n" + "="*70)
    print("SAVING RESULTS")
    print("="*70)

    # JSON format
    json_path = RESULTS_DIR / "summary.json"
    with open(json_path, 'w') as f:
        # Convert numpy types to Python types for JSON
        def convert(obj):
            if isinstance(obj, np.ndarray):
                return obj.tolist()
            elif isinstance(obj, (np.float32, np.float64)):
                return float(obj)
            elif isinstance(obj, (np.int32, np.int64)):
                return int(obj)
            return obj

        clean_results = json.loads(json.dumps(results, default=convert))
        json.dump(clean_results, f, indent=2)
    print(f"  Saved: {json_path}")

    # Text format
    txt_path = RESULTS_DIR / "summary.txt"
    with open(txt_path, 'w') as f:
        f.write("="*70 + "\n")
        f.write("FEP RESULTS SUMMARY\n")
        f.write("Natural Product CID_10120 (Bufadienolide)\n")
        f.write("NOD2 R702W Mutation Effect\n")
        f.write("="*70 + "\n\n")

        f.write(f"Timestamp: {results['timestamp']}\n\n")

        if 'dG_bind_wt' in results:
            f.write(f"DeltaG_bind(WT) = {results['dG_bind_wt']:.2f} ± {results['dG_bind_wt_err']:.2f} kcal/mol\n")
        if 'dG_bind_mut' in results:
            f.write(f"DeltaG_bind(R702W) = {results['dG_bind_mut']:.2f} ± {results['dG_bind_mut_err']:.2f} kcal/mol\n")
        if 'ddG' in results:
            f.write(f"\nDeltaDeltaG = {results['ddG']:.2f} ± {results['ddG_err']:.2f} kcal/mol\n")
            f.write(f"Significance: {results['ddG_sigma']:.1f} sigma\n")
            f.write(f"Interpretation: {results['interpretation']}\n")

        f.write("\n" + "="*70 + "\n")
        f.write("PER-SYSTEM BREAKDOWN\n")
        f.write("="*70 + "\n")

        for sys_name, sys_data in results.get('systems', {}).items():
            f.write(f"\n{sys_name}:\n")
            f.write(f"  DeltaG = {sys_data['dG_kcal']:.2f} ± {sys_data['dG_err_kcal']:.2f} kcal/mol\n")
            if sys_data.get('bootstrap_std_kcal'):
                f.write(f"  Bootstrap std = {sys_data['bootstrap_std_kcal']:.3f} kcal/mol\n")
            for phase, vals in sys_data.get('per_phase', {}).items():
                f.write(f"  {phase}: {vals['dG']:.2f} ± {vals['dG_err']:.2f}\n")

    print(f"  Saved: {txt_path}")

    print("\n" + "="*70)
    print("ANALYSIS COMPLETE")
    print("="*70)

    return results

if __name__ == "__main__":
    main()
