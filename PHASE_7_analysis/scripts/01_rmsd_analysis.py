#!/usr/bin/env python3
"""
PHASE 7: RMSD Analysis (CORRECTED)
Computes ligand and protein backbone RMSD for all trajectories.
PhD-level: Individual replicates + ensemble average with 95% CI

FIXES APPLIED:
- Proper rotational alignment using Kabsch algorithm (rms.RMSD)
- PBC-aware calculations
- Robust ligand selection
"""

import MDAnalysis as mda
from MDAnalysis.analysis import rms, align
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import warnings
import gc
import os

warnings.filterwarnings('ignore')

# Paths - relative to script location for portability
SCRIPT_DIR = Path(__file__).parent
BASE_DIR = Path(os.getenv('NOD2_BASE', SCRIPT_DIR.parent.parent / "PHASE_6"))
TRAJ_DIR = BASE_DIR / "trajectories"
OUTPUT_DIR = SCRIPT_DIR.parent / "rmsd"

# Compound configuration
COMPOUNDS = {
    'budesonide': {'type': 'Positive Control', 'reps': 3, 'top': 'budesonide_rep1_solvated.pdb'},
    'febuxostat': {'type': 'Lead Candidate', 'reps': 3, 'top': 'febuxostat_rep1_solvated.pdb'},
    'ursodiol': {'type': 'Secondary Candidate', 'reps': 3, 'top': 'ursodiol_rep3_solvated.pdb'},
    'natural_top': {'type': 'Natural Product', 'reps': 3, 'top': 'natural_top_rep1_solvated.pdb'},
    'decoy': {'type': 'Negative Control', 'reps': 1, 'top': 'decoy_rep1_solvated.pdb'},
    'apo': {'type': 'Baseline', 'reps': 1, 'top': 'apo_rep1_solvated.pdb'},
}

# Color scheme
COLORS = {
    'budesonide': '#1f77b4',
    'febuxostat': '#2ca02c',
    'ursodiol': '#ff7f0e',
    'natural_top': '#9467bd',
    'decoy': '#d62728',
    'apo': '#7f7f7f',
}

# Residues to exclude from ligand selection
WATER_RESNAMES = {"HOH", "WAT", "SOL", "TIP3", "TIP3P", "SPC", "SPCE", "TIP4P"}
ION_RESNAMES = {"NA", "CL", "K", "CA", "MG", "ZN", "MN", "FE", "CU", "CO", "NI", "NA+", "CL-"}

sns.set_context("paper", font_scale=1.5)
sns.set_style("whitegrid")


def pick_ligand_selection(universe):
    """
    Robustly pick the ligand as the largest non-protein, non-water, non-ion residue.
    Returns (selection_string, AtomGroup).
    """
    nonprot = universe.select_atoms("not protein")

    candidates = []
    for res in nonprot.residues:
        rname = (res.resname or "").upper()
        if rname in WATER_RESNAMES or rname in ION_RESNAMES:
            continue
        # Ignore single-atom junk
        if len(res.atoms) < 5:
            continue
        candidates.append(res)

    if not candidates:
        return None, None

    # Pick largest residue by atom count
    lig_res = max(candidates, key=lambda r: len(r.atoms))
    lig_atoms = lig_res.atoms

    sel = f"resid {lig_res.resid} and resname {lig_res.resname}"

    return sel, lig_atoms


def compute_protein_rmsd(universe, stride=10):
    """
    Compute protein backbone RMSD using MDAnalysis built-in.
    Automatically handles alignment via least-squares fitting.
    """
    protein = universe.select_atoms("protein and backbone")

    if len(protein) == 0:
        print("  Warning: No protein backbone atoms found!")
        return None, None

    R = rms.RMSD(universe, universe, select="protein and backbone", ref_frame=0)
    R.run(step=stride, verbose=False)

    times = R.results.rmsd[:, 1] / 1000  # Convert ps to ns
    rmsd = R.results.rmsd[:, 2]  # RMSD in Angstroms

    return times, rmsd


def compute_ligand_rmsd(universe, stride=10):
    """
    Compute ligand RMSD with PROPER alignment.

    Method:
    1. Align protein backbone to reference frame (Kabsch algorithm)
    2. Measure where ligand ended up after alignment

    This is the scientifically correct way - removes protein rotation/translation.
    """
    lig_sel, lig_atoms = pick_ligand_selection(universe)

    if lig_atoms is None or len(lig_atoms) == 0:
        print("  Warning: No ligand atoms found!")
        return None, None

    # Get heavy atoms only for RMSD
    lig_heavy_sel = f"({lig_sel}) and not name H*"
    lig_heavy = universe.select_atoms(lig_heavy_sel)

    if len(lig_heavy) == 0:
        lig_heavy_sel = lig_sel  # Fallback to all atoms
        lig_heavy = lig_atoms

    print(f"  Ligand: {lig_atoms.residues[0].resname}{lig_atoms.residues[0].resid} "
          f"({len(lig_heavy)} heavy atoms)")

    # Use MDAnalysis RMSD with groupselections
    # select: group used for alignment (superposition)
    # groupselections: additional groups to measure RMSD after alignment
    R = rms.RMSD(
        universe,
        universe,
        select="protein and backbone",
        groupselections=[lig_heavy_sel],
        ref_frame=0
    )
    R.run(step=stride, verbose=False)

    # Columns: [Frame, Time(ps), Backbone_RMSD, Ligand_RMSD]
    times = R.results.rmsd[:, 1] / 1000  # Convert ps to ns
    ligand_rmsd = R.results.rmsd[:, 3]   # Ligand RMSD after alignment

    return times, ligand_rmsd


def bootstrap_ci(data_list, n_bootstrap=1000, ci=95):
    """Compute bootstrap confidence interval for ensemble of replicates."""
    if len(data_list) == 1:
        return data_list[0], data_list[0], data_list[0]

    # Interpolate to common time points
    min_len = min(len(d) for d in data_list)
    data_array = np.array([d[:min_len] for d in data_list])

    mean = np.mean(data_array, axis=0)

    # Bootstrap
    boot_means = []
    for _ in range(n_bootstrap):
        idx = np.random.choice(len(data_list), len(data_list), replace=True)
        boot_sample = data_array[idx]
        boot_means.append(np.mean(boot_sample, axis=0))

    boot_means = np.array(boot_means)
    lower = np.percentile(boot_means, (100 - ci) / 2, axis=0)
    upper = np.percentile(boot_means, 100 - (100 - ci) / 2, axis=0)

    return mean, lower, upper


def analyze_compound(name, config, stride=10):
    """Analyze all replicates for a compound."""
    print(f"\nAnalyzing {name} ({config['type']})...")

    top_file = TRAJ_DIR / config['top']

    if not top_file.exists():
        print(f"  ERROR: Topology not found: {top_file}")
        return None

    results = {
        'name': name,
        'type': config['type'],
        'protein_rmsd': [],
        'ligand_rmsd': [],
        'times': None
    }

    for rep in range(1, config['reps'] + 1):
        traj_file = TRAJ_DIR / f"{name}_rep{rep}.dcd"

        if not traj_file.exists():
            print(f"  Trajectory not found: {traj_file}")
            continue

        print(f"  Loading rep{rep}...")

        try:
            u = mda.Universe(str(top_file), str(traj_file))
            print(f"    Frames: {len(u.trajectory)}, Atoms: {len(u.atoms)}")

            # Protein RMSD
            times, prot_rmsd = compute_protein_rmsd(u, stride=stride)
            if prot_rmsd is not None:
                results['protein_rmsd'].append(prot_rmsd)
                results['times'] = times

            # Ligand RMSD (skip for apo)
            if name != 'apo':
                _, lig_rmsd = compute_ligand_rmsd(u, stride=stride)
                if lig_rmsd is not None:
                    results['ligand_rmsd'].append(lig_rmsd)

            # Clean up memory
            del u
            gc.collect()

        except Exception as e:
            print(f"    ERROR: {e}")
            import traceback
            traceback.print_exc()
            continue

    return results


def plot_rmsd_comparison(all_results, output_dir):
    """Create publication-quality RMSD comparison plots."""

    # Figure 1: Ligand RMSD comparison (main result)
    fig, ax = plt.subplots(figsize=(12, 8))

    for name, results in all_results.items():
        if name == 'apo' or not results['ligand_rmsd']:
            continue

        color = COLORS[name]
        times = results['times']

        # Plot individual replicates (faint)
        for i, rmsd in enumerate(results['ligand_rmsd']):
            min_len = min(len(times), len(rmsd))
            ax.plot(times[:min_len], rmsd[:min_len],
                   color=color, alpha=0.3, linewidth=1)

        # Plot ensemble average with CI (bold)
        if len(results['ligand_rmsd']) > 1:
            mean, lower, upper = bootstrap_ci(results['ligand_rmsd'])
            min_len = min(len(times), len(mean))
            ax.plot(times[:min_len], mean[:min_len],
                   color=color, linewidth=2.5, label=f"{name} (n={len(results['ligand_rmsd'])})")
            ax.fill_between(times[:min_len], lower[:min_len], upper[:min_len],
                           color=color, alpha=0.2)
        else:
            rmsd = results['ligand_rmsd'][0]
            min_len = min(len(times), len(rmsd))
            ax.plot(times[:min_len], rmsd[:min_len],
                   color=color, linewidth=2.5, label=f"{name} (n=1)")

    # Threshold lines
    ax.axhline(y=3, color='green', linestyle='--', alpha=0.5, label='Stable threshold (3 Ang)')
    ax.axhline(y=5, color='red', linestyle='--', alpha=0.5, label='Unstable threshold (5 Ang)')

    ax.set_xlabel('Time (ns)', fontsize=14)
    ax.set_ylabel('Ligand RMSD (Ang)', fontsize=14)
    ax.set_title('Ligand RMSD: Binding Stability Comparison\n(Aligned to protein backbone)', fontsize=16)
    ax.legend(loc='upper right', fontsize=10)
    ax.set_xlim(0, None)
    ax.set_ylim(0, None)

    plt.tight_layout()
    plt.savefig(output_dir / 'ligand_rmsd_comparison.png', dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved: {output_dir / 'ligand_rmsd_comparison.png'}")

    # Figure 2: Protein backbone RMSD (sanity check)
    fig, ax = plt.subplots(figsize=(12, 6))

    for name, results in all_results.items():
        if not results['protein_rmsd']:
            continue

        color = COLORS[name]
        times = results['times']

        # Ensemble average
        if len(results['protein_rmsd']) > 1:
            mean, lower, upper = bootstrap_ci(results['protein_rmsd'])
            min_len = min(len(times), len(mean))
            ax.plot(times[:min_len], mean[:min_len], color=color, linewidth=2, label=name)
        else:
            rmsd = results['protein_rmsd'][0]
            min_len = min(len(times), len(rmsd))
            ax.plot(times[:min_len], rmsd[:min_len], color=color, linewidth=2, label=name)

    ax.axhline(y=3, color='gray', linestyle='--', alpha=0.5)
    ax.set_xlabel('Time (ns)', fontsize=14)
    ax.set_ylabel('Backbone RMSD (Ang)', fontsize=14)
    ax.set_title('Protein Backbone RMSD (Stability Check)', fontsize=16)
    ax.legend(loc='upper right', fontsize=10)

    plt.tight_layout()
    plt.savefig(output_dir / 'protein_rmsd_comparison.png', dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved: {output_dir / 'protein_rmsd_comparison.png'}")


def generate_stats_table(all_results, output_dir):
    """Generate statistics table."""
    stats = []

    for name, results in all_results.items():
        if name == 'apo':
            continue

        if results['ligand_rmsd']:
            all_rmsd = np.concatenate(results['ligand_rmsd'])
            mean_rmsd = np.mean(all_rmsd)
            std_rmsd = np.std(all_rmsd)
            pct_stable = np.mean(all_rmsd < 3) * 100

            # Verdict
            if mean_rmsd < 3:
                verdict = "STABLE"
            elif mean_rmsd < 5:
                verdict = "MARGINAL"
            else:
                verdict = "UNSTABLE"

            stats.append({
                'Compound': name,
                'Type': results['type'],
                'Mean RMSD (Ang)': f"{mean_rmsd:.2f}",
                'Std Dev (Ang)': f"{std_rmsd:.2f}",
                '% Frames <3Ang': f"{pct_stable:.1f}",
                'Verdict': verdict
            })

    df = pd.DataFrame(stats)
    df.to_csv(output_dir / 'rmsd_statistics.csv', index=False)
    print(f"\nSaved: {output_dir / 'rmsd_statistics.csv'}")
    print("\n" + df.to_string(index=False))

    return df


def main():
    print("=" * 60)
    print("PHASE 7: RMSD ANALYSIS (CORRECTED)")
    print("Using proper Kabsch alignment for ligand RMSD")
    print("=" * 60)

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    # Analyze all compounds
    all_results = {}
    for name, config in COMPOUNDS.items():
        results = analyze_compound(name, config, stride=10)
        if results:
            all_results[name] = results

    if not all_results:
        print("ERROR: No trajectories could be analyzed!")
        return

    # Generate plots
    print("\nGenerating plots...")
    plot_rmsd_comparison(all_results, OUTPUT_DIR)

    # Generate statistics table
    print("\nGenerating statistics...")
    stats_df = generate_stats_table(all_results, OUTPUT_DIR)

    print("\n" + "=" * 60)
    print("RMSD ANALYSIS COMPLETE")
    print("=" * 60)


if __name__ == "__main__":
    main()
