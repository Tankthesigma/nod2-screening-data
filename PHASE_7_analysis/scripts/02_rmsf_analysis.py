#!/usr/bin/env python3
"""
PHASE 7: RMSF Analysis
Computes per-residue root mean square fluctuation.
Shows which protein regions are flexible vs rigid.
"""

import MDAnalysis as mda
from MDAnalysis.analysis import rms
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import warnings
import gc

warnings.filterwarnings('ignore')

# Paths
BASE_DIR = Path(r"C:\Users\vasud\nod2-screening-data\PHASE_6")
TRAJ_DIR = BASE_DIR / "trajectories"
OUTPUT_DIR = Path(r"C:\Users\vasud\nod2-screening-data\PHASE_7_analysis\rmsf")

# Compound configuration
COMPOUNDS = {
    'budesonide': {'type': 'Positive Control', 'reps': 3, 'top': 'budesonide_rep1_solvated.pdb'},
    'febuxostat': {'type': 'Lead Candidate', 'reps': 3, 'top': 'febuxostat_rep1_solvated.pdb'},
    'ursodiol': {'type': 'Secondary Candidate', 'reps': 3, 'top': 'ursodiol_rep3_solvated.pdb'},
    'natural_top': {'type': 'Natural Product', 'reps': 3, 'top': 'natural_top_rep1_solvated.pdb'},
    'decoy': {'type': 'Negative Control', 'reps': 1, 'top': 'decoy_rep1_solvated.pdb'},
    'apo': {'type': 'Baseline', 'reps': 1, 'top': 'apo_rep1_solvated.pdb'},
}

COLORS = {
    'budesonide': '#1f77b4',
    'febuxostat': '#2ca02c',
    'ursodiol': '#ff7f0e',
    'natural_top': '#9467bd',
    'decoy': '#d62728',
    'apo': '#7f7f7f',
}

sns.set_context("paper", font_scale=1.5)
sns.set_style("whitegrid")


def compute_rmsf(universe, stride=10):
    """Compute per-residue RMSF for CA atoms."""
    protein_ca = universe.select_atoms("protein and name CA")

    if len(protein_ca) == 0:
        print("  Warning: No CA atoms found!")
        return None, None

    # Collect positions
    positions = []
    for ts in universe.trajectory[::stride]:
        positions.append(protein_ca.positions.copy())

    positions = np.array(positions)

    # Compute mean position
    mean_pos = positions.mean(axis=0)

    # Compute RMSF
    diff = positions - mean_pos
    rmsf = np.sqrt((diff ** 2).sum(axis=2).mean(axis=0))

    # Get residue IDs
    resids = protein_ca.resids

    return resids, rmsf


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
        'rmsf_values': [],
        'resids': None
    }

    for rep in range(1, config['reps'] + 1):
        traj_file = TRAJ_DIR / f"{name}_rep{rep}.dcd"

        if not traj_file.exists():
            print(f"  Trajectory not found: {traj_file}")
            continue

        print(f"  Loading rep{rep}...")

        try:
            u = mda.Universe(str(top_file), str(traj_file))

            resids, rmsf = compute_rmsf(u, stride=stride)

            if rmsf is not None:
                results['rmsf_values'].append(rmsf)
                results['resids'] = resids

            del u
            gc.collect()

        except Exception as e:
            print(f"    ERROR: {e}")
            continue

    return results


def plot_rmsf_comparison(all_results, output_dir):
    """Create RMSF comparison plot."""
    fig, ax = plt.subplots(figsize=(14, 6))

    for name, results in all_results.items():
        if not results['rmsf_values']:
            continue

        color = COLORS[name]
        resids = results['resids']

        # Ensemble average
        if len(results['rmsf_values']) > 1:
            rmsf_array = np.array(results['rmsf_values'])
            mean_rmsf = np.mean(rmsf_array, axis=0)
            std_rmsf = np.std(rmsf_array, axis=0)

            ax.plot(resids, mean_rmsf, color=color, linewidth=2, label=name)
            ax.fill_between(resids, mean_rmsf - std_rmsf, mean_rmsf + std_rmsf,
                           color=color, alpha=0.2)
        else:
            ax.plot(resids, results['rmsf_values'][0], color=color, linewidth=2, label=name)

    ax.set_xlabel('Residue Number', fontsize=14)
    ax.set_ylabel('RMSF (Å)', fontsize=14)
    ax.set_title('Per-Residue Flexibility (RMSF)', fontsize=16)
    ax.legend(loc='upper right', fontsize=10)

    plt.tight_layout()
    plt.savefig(output_dir / 'rmsf_comparison.png', dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved: {output_dir / 'rmsf_comparison.png'}")


def plot_rmsf_heatmap(all_results, output_dir):
    """Create RMSF heatmap across compounds."""
    # Collect data for heatmap
    data = {}
    resids = None

    for name, results in all_results.items():
        if results['rmsf_values']:
            mean_rmsf = np.mean(results['rmsf_values'], axis=0)
            data[name] = mean_rmsf
            if resids is None:
                resids = results['resids']

    if not data:
        return

    # Create DataFrame
    df = pd.DataFrame(data, index=resids)

    # Sample every 10 residues for cleaner heatmap
    df_sampled = df.iloc[::10]

    fig, ax = plt.subplots(figsize=(10, 8))
    sns.heatmap(df_sampled.T, cmap='YlOrRd', ax=ax,
                xticklabels=20, cbar_kws={'label': 'RMSF (Å)'})

    ax.set_xlabel('Residue Number', fontsize=14)
    ax.set_ylabel('Compound', fontsize=14)
    ax.set_title('RMSF Heatmap: Flexibility by Compound', fontsize=16)

    plt.tight_layout()
    plt.savefig(output_dir / 'rmsf_heatmap.png', dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved: {output_dir / 'rmsf_heatmap.png'}")


def save_rmsf_data(all_results, output_dir):
    """Save RMSF data to CSV."""
    data = {'Residue': None}

    for name, results in all_results.items():
        if results['rmsf_values']:
            mean_rmsf = np.mean(results['rmsf_values'], axis=0)
            data[name] = mean_rmsf
            if data['Residue'] is None:
                data['Residue'] = results['resids']

    if data['Residue'] is not None:
        df = pd.DataFrame(data)
        df.to_csv(output_dir / 'rmsf_data.csv', index=False)
        print(f"Saved: {output_dir / 'rmsf_data.csv'}")


def main():
    print("=" * 60)
    print("PHASE 7: RMSF ANALYSIS")
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
    plot_rmsf_comparison(all_results, OUTPUT_DIR)
    plot_rmsf_heatmap(all_results, OUTPUT_DIR)

    # Save data
    save_rmsf_data(all_results, OUTPUT_DIR)

    print("\n" + "=" * 60)
    print("RMSF ANALYSIS COMPLETE")
    print("=" * 60)


if __name__ == "__main__":
    main()
