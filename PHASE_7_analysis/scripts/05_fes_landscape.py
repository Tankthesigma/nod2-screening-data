#!/usr/bin/env python3
"""
PHASE 7: Free Energy Surface (FES) Proxy
Creates 2D density plot of Ligand RMSD vs Radius of Gyration.
PhD-level: Shows conformational stability landscape.
"""

import MDAnalysis as mda
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
OUTPUT_DIR = Path(r"C:\Users\vasud\nod2-screening-data\PHASE_7_analysis\fes")

# Compound configuration (skip apo)
COMPOUNDS = {
    'budesonide': {'type': 'Positive Control', 'reps': 3, 'top': 'budesonide_rep1_solvated.pdb'},
    'febuxostat': {'type': 'Lead Candidate', 'reps': 3, 'top': 'febuxostat_rep1_solvated.pdb'},
    'ursodiol': {'type': 'Secondary Candidate', 'reps': 3, 'top': 'ursodiol_rep3_solvated.pdb'},
    'natural_top': {'type': 'Natural Product', 'reps': 3, 'top': 'natural_top_rep1_solvated.pdb'},
    'decoy': {'type': 'Negative Control', 'reps': 1, 'top': 'decoy_rep1_solvated.pdb'},
}

COLORS = {
    'budesonide': '#1f77b4',
    'febuxostat': '#2ca02c',
    'ursodiol': '#ff7f0e',
    'natural_top': '#9467bd',
    'decoy': '#d62728',
}

sns.set_context("paper", font_scale=1.5)
sns.set_style("whitegrid")


def compute_ligand_rmsd_rg(universe, stride=10):
    """
    Compute ligand RMSD and radius of gyration over trajectory.
    """
    protein_bb = universe.select_atoms("protein and backbone")
    ligand = universe.select_atoms("not protein and not resname HOH WAT SOL NA CL")

    if len(ligand) == 0:
        print("  Warning: No ligand atoms found!")
        return None, None

    print(f"  Ligand atoms: {len(ligand)}")

    # Reference
    universe.trajectory[0]
    ref_ligand = ligand.positions.copy()
    ref_ligand_com = ref_ligand.mean(axis=0)
    ref_ligand_centered = ref_ligand - ref_ligand_com

    rmsd_values = []
    rg_values = []

    for ts in universe.trajectory[::stride]:
        # Ligand RMSD
        lig_pos = ligand.positions.copy()
        lig_com = lig_pos.mean(axis=0)
        lig_centered = lig_pos - lig_com

        # Simple RMSD
        diff = lig_centered - ref_ligand_centered
        rmsd = np.sqrt((diff ** 2).sum() / len(diff))
        rmsd_values.append(rmsd)

        # Radius of gyration
        rg = np.sqrt(np.mean(np.sum((lig_pos - lig_com) ** 2, axis=1)))
        rg_values.append(rg)

    return np.array(rmsd_values), np.array(rg_values)


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
        'rmsd': [],
        'rg': []
    }

    for rep in range(1, config['reps'] + 1):
        traj_file = TRAJ_DIR / f"{name}_rep{rep}.dcd"

        if not traj_file.exists():
            print(f"  Trajectory not found: {traj_file}")
            continue

        print(f"  Loading rep{rep}...")

        try:
            u = mda.Universe(str(top_file), str(traj_file))

            rmsd, rg = compute_ligand_rmsd_rg(u, stride=stride)

            if rmsd is not None:
                results['rmsd'].extend(rmsd)
                results['rg'].extend(rg)

            del u
            gc.collect()

        except Exception as e:
            print(f"    ERROR: {e}")
            continue

    return results


def plot_individual_fes(all_results, output_dir):
    """Create individual FES plots for each compound."""

    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    axes = axes.flatten()

    for idx, (name, results) in enumerate(all_results.items()):
        if idx >= 6:
            break

        ax = axes[idx]

        if not results['rmsd']:
            ax.set_visible(False)
            continue

        rmsd = np.array(results['rmsd'])
        rg = np.array(results['rg'])

        # Create 2D histogram / KDE
        try:
            sns.kdeplot(x=rmsd, y=rg, ax=ax, cmap='Blues', fill=True, levels=20, alpha=0.8)
        except:
            ax.hist2d(rmsd, rg, bins=30, cmap='Blues')

        ax.scatter(rmsd[0], rg[0], color='red', s=100, marker='*',
                  label='Start', zorder=5)
        ax.scatter(rmsd[-1], rg[-1], color='green', s=100, marker='s',
                  label='End', zorder=5)

        ax.axvline(x=3, color='orange', linestyle='--', alpha=0.7, label='Stable (3Å)')
        ax.axvline(x=5, color='red', linestyle='--', alpha=0.7, label='Unstable (5Å)')

        ax.set_xlabel('Ligand RMSD (Å)')
        ax.set_ylabel('Radius of Gyration (Å)')
        ax.set_title(f'{name}', fontsize=14, fontweight='bold')
        ax.legend(loc='upper right', fontsize=8)

    # Hide unused axes
    for idx in range(len(all_results), 6):
        axes[idx].set_visible(False)

    plt.suptitle('Free Energy Surface Proxy: RMSD vs Rg', fontsize=16, fontweight='bold')
    plt.tight_layout()
    plt.savefig(output_dir / 'fes_individual.png', dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved: {output_dir / 'fes_individual.png'}")


def plot_combined_scatter(all_results, output_dir):
    """Create combined scatter plot showing all compounds."""

    fig, ax = plt.subplots(figsize=(12, 8))

    for name, results in all_results.items():
        if not results['rmsd']:
            continue

        rmsd = np.array(results['rmsd'])
        rg = np.array(results['rg'])
        color = COLORS[name]

        # Subsample for clarity
        step = max(1, len(rmsd) // 500)
        ax.scatter(rmsd[::step], rg[::step], c=color, alpha=0.3, s=10, label=name)

        # Plot centroid
        ax.scatter(np.mean(rmsd), np.mean(rg), c=color, s=200, marker='X',
                  edgecolors='black', linewidths=2)

    ax.axvline(x=3, color='green', linestyle='--', alpha=0.5, label='Stable threshold')
    ax.axvline(x=5, color='red', linestyle='--', alpha=0.5, label='Unstable threshold')

    ax.set_xlabel('Ligand RMSD (Å)', fontsize=14)
    ax.set_ylabel('Radius of Gyration (Å)', fontsize=14)
    ax.set_title('Conformational Landscape: All Compounds', fontsize=16)
    ax.legend(loc='upper right', fontsize=10)

    plt.tight_layout()
    plt.savefig(output_dir / 'fes_combined.png', dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved: {output_dir / 'fes_combined.png'}")


def generate_fes_stats(all_results, output_dir):
    """Generate FES statistics."""

    stats = []
    for name, results in all_results.items():
        if results['rmsd']:
            rmsd = np.array(results['rmsd'])
            rg = np.array(results['rg'])

            stats.append({
                'Compound': name,
                'Mean RMSD (Å)': f"{np.mean(rmsd):.2f}",
                'Std RMSD (Å)': f"{np.std(rmsd):.2f}",
                'Mean Rg (Å)': f"{np.mean(rg):.2f}",
                'Std Rg (Å)': f"{np.std(rg):.2f}",
                '% RMSD < 3Å': f"{np.mean(rmsd < 3) * 100:.1f}",
            })

    df = pd.DataFrame(stats)
    df.to_csv(output_dir / 'fes_statistics.csv', index=False)
    print(f"\nSaved: {output_dir / 'fes_statistics.csv'}")
    print("\n" + df.to_string(index=False))


def main():
    print("=" * 60)
    print("PHASE 7: FREE ENERGY SURFACE ANALYSIS")
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

    # Generate outputs
    print("\nGenerating plots...")
    plot_individual_fes(all_results, OUTPUT_DIR)
    plot_combined_scatter(all_results, OUTPUT_DIR)
    generate_fes_stats(all_results, OUTPUT_DIR)

    print("\n" + "=" * 60)
    print("FES ANALYSIS COMPLETE")
    print("=" * 60)


if __name__ == "__main__":
    main()
