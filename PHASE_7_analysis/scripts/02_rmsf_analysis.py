#!/usr/bin/env python3
"""
PHASE 7: RMSF Analysis (CORRECTED)
Computes per-residue root mean square fluctuation.
Shows which protein regions are flexible vs rigid.

FIXES APPLIED:
- Trajectory alignment to reference frame before RMSF calculation
- Streaming computation to avoid memory issues
- Proper use of MDAnalysis AlignTraj
"""

import MDAnalysis as mda
from MDAnalysis.analysis import align, rms
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
OUTPUT_DIR = SCRIPT_DIR.parent / "rmsf"

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
    """
    Compute per-residue RMSF for CA atoms with PROPER alignment.

    Method:
    1. Create reference structure (first frame)
    2. Align each frame to reference using backbone CA
    3. Compute RMSF on aligned trajectory

    Uses streaming (Welford's algorithm) to avoid memory issues.
    """
    protein_ca = universe.select_atoms("protein and name CA")
    protein_bb = universe.select_atoms("protein and backbone")

    if len(protein_ca) == 0:
        print("  Warning: No CA atoms found!")
        return None, None

    # Create reference - use first frame
    ref = universe.copy()
    ref.trajectory[0]
    ref_ca = ref.select_atoms("protein and name CA")

    # Streaming RMSF calculation (Welford's online algorithm)
    # This avoids storing all frames in memory
    n = 0
    mean = np.zeros((len(protein_ca), 3), dtype=np.float64)
    M2 = np.zeros((len(protein_ca), 3), dtype=np.float64)

    print(f"  Computing RMSF with alignment...")

    for ts in universe.trajectory[::stride]:
        # Align current frame to reference (Kabsch algorithm)
        # This modifies positions in place
        align.alignto(universe, ref, select="protein and name CA")

        # Get aligned CA positions
        x = protein_ca.positions.astype(np.float64)

        # Welford's online algorithm for variance
        n += 1
        delta = x - mean
        mean += delta / n
        delta2 = x - mean
        M2 += delta * delta2

    if n == 0:
        return None, None

    # Population variance over frames, then RMSF
    var = M2 / n
    rmsf = np.sqrt(var.sum(axis=1))  # Sum over x,y,z then sqrt

    return protein_ca.resids, rmsf


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
            print(f"    Frames: {len(u.trajectory)}, CA atoms: {len(u.select_atoms('protein and name CA'))}")

            resids, rmsf = compute_rmsf(u, stride=stride)

            if rmsf is not None:
                results['rmsf_values'].append(rmsf)
                results['resids'] = resids
                print(f"    Mean RMSF: {np.mean(rmsf):.2f} Ang")

            del u
            gc.collect()

        except Exception as e:
            print(f"    ERROR: {e}")
            import traceback
            traceback.print_exc()
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
    ax.set_ylabel('RMSF (Ang)', fontsize=14)
    ax.set_title('Per-Residue Flexibility (RMSF)\n(Aligned to first frame)', fontsize=16)
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
                xticklabels=20, cbar_kws={'label': 'RMSF (Ang)'})

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
    print("PHASE 7: RMSF ANALYSIS (CORRECTED)")
    print("Using trajectory alignment for accurate flexibility")
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
