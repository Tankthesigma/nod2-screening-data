#!/usr/bin/env python3
"""
PHASE 7: RMSD Analysis
Computes ligand and protein backbone RMSD for all trajectories.
PhD-level: Individual replicates + ensemble average with 95% CI
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
OUTPUT_DIR = Path(r"C:\Users\vasud\nod2-screening-data\PHASE_7_analysis\rmsd")

# Compound configuration: (name, type, num_reps, topology_to_use)
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
    'budesonide': '#1f77b4',   # blue
    'febuxostat': '#2ca02c',   # green
    'ursodiol': '#ff7f0e',     # orange
    'natural_top': '#9467bd',  # purple
    'decoy': '#d62728',        # red
    'apo': '#7f7f7f',          # gray
}

# Plotting style
sns.set_context("paper", font_scale=1.5)
sns.set_style("whitegrid")


def compute_rmsd(universe, selection, ref_frame=0, stride=10):
    """
    Compute RMSD for a selection over trajectory.
    Uses strided loading for memory efficiency.
    """
    # Set reference to first frame
    universe.trajectory[ref_frame]
    ref_positions = universe.select_atoms(selection).positions.copy()
    ref_com = ref_positions.mean(axis=0)
    ref_positions -= ref_com

    rmsd_values = []
    times = []

    for ts in universe.trajectory[::stride]:
        current_positions = universe.select_atoms(selection).positions.copy()
        current_com = current_positions.mean(axis=0)
        current_positions -= current_com

        # Simple RMSD calculation
        diff = current_positions - ref_positions
        rmsd = np.sqrt((diff ** 2).sum() / len(diff))
        rmsd_values.append(rmsd)
        times.append(ts.time / 1000)  # Convert ps to ns

    return np.array(times), np.array(rmsd_values)


def compute_protein_rmsd(universe, stride=10):
    """Compute protein backbone RMSD using MDAnalysis built-in."""
    protein = universe.select_atoms("protein and backbone")

    R = rms.RMSD(protein, protein, select="backbone", ref_frame=0)
    R.run(step=stride, verbose=False)

    times = R.results.rmsd[:, 1] / 1000  # Convert to ns
    rmsd = R.results.rmsd[:, 2]  # RMSD in Angstroms

    return times, rmsd


def compute_ligand_rmsd(universe, stride=10):
    """Compute ligand RMSD after aligning protein backbone."""
    protein = universe.select_atoms("protein and backbone")

    # Try different ligand selections
    ligand = universe.select_atoms("not protein and not resname HOH WAT SOL NA CL")

    if len(ligand) == 0:
        print("  Warning: No ligand atoms found!")
        return None, None

    print(f"  Ligand atoms: {len(ligand)}")

    # Reference frame
    universe.trajectory[0]
    ref_protein = protein.positions.copy()
    ref_ligand = ligand.positions.copy()

    rmsd_values = []
    times = []

    for ts in universe.trajectory[::stride]:
        # Align protein to reference
        mobile_protein = protein.positions

        # Simple alignment (translate to COM)
        ref_com = ref_protein.mean(axis=0)
        mob_com = mobile_protein.mean(axis=0)

        # Get ligand positions and compute RMSD
        lig_positions = ligand.positions - mob_com + ref_com
        diff = lig_positions - ref_ligand
        rmsd = np.sqrt((diff ** 2).sum() / len(diff))

        rmsd_values.append(rmsd)
        times.append(ts.time / 1000)

    return np.array(times), np.array(rmsd_values)


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
    ax.axhline(y=3, color='green', linestyle='--', alpha=0.5, label='Stable threshold (3Å)')
    ax.axhline(y=5, color='red', linestyle='--', alpha=0.5, label='Unstable threshold (5Å)')

    ax.set_xlabel('Time (ns)', fontsize=14)
    ax.set_ylabel('Ligand RMSD (Å)', fontsize=14)
    ax.set_title('Ligand RMSD: Binding Stability Comparison', fontsize=16)
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
    ax.set_ylabel('Backbone RMSD (Å)', fontsize=14)
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
                'Mean RMSD (Å)': f"{mean_rmsd:.2f}",
                'Std Dev (Å)': f"{std_rmsd:.2f}",
                '% Frames <3Å': f"{pct_stable:.1f}",
                'Verdict': verdict
            })

    df = pd.DataFrame(stats)
    df.to_csv(output_dir / 'rmsd_statistics.csv', index=False)
    print(f"\nSaved: {output_dir / 'rmsd_statistics.csv'}")
    print("\n" + df.to_string(index=False))

    return df


def main():
    print("=" * 60)
    print("PHASE 7: RMSD ANALYSIS")
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
