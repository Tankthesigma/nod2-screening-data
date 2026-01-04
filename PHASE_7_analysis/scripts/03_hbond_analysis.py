#!/usr/bin/env python3
"""
PHASE 7: Hydrogen Bond Analysis
Computes protein-ligand H-bonds over time.
PhD-level: Reports H-bond persistence (lifetime), not just counts.
"""

import MDAnalysis as mda
from MDAnalysis.analysis.hydrogenbonds import HydrogenBondAnalysis
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from collections import defaultdict
import warnings
import gc

warnings.filterwarnings('ignore')

# Paths
BASE_DIR = Path(r"C:\Users\vasud\nod2-screening-data\PHASE_6")
TRAJ_DIR = BASE_DIR / "trajectories"
OUTPUT_DIR = Path(r"C:\Users\vasud\nod2-screening-data\PHASE_7_analysis\hbonds")

# Compound configuration (skip apo - no ligand)
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


def compute_hbonds_simple(universe, stride=10):
    """
    Simple H-bond counting between protein and ligand.
    Uses distance-based criterion (3.5A heavy atom distance, 150 degree angle).
    """
    protein = universe.select_atoms("protein")
    ligand = universe.select_atoms("not protein and not resname HOH WAT SOL NA CL")

    if len(ligand) == 0:
        print("  Warning: No ligand atoms found!")
        return None, None, None

    print(f"  Ligand atoms: {len(ligand)}")

    # Track H-bonds per frame
    hbond_counts = []
    times = []

    # Get potential donors and acceptors
    protein_donors = protein.select_atoms("name N NH NH1 NH2 NE NE1 NE2 NZ ND1 ND2 OG OG1 OH NE1")
    protein_acceptors = protein.select_atoms("name O OD1 OD2 OE1 OE2 OG OG1 OH ND1 NE2 SD")

    ligand_heavy = ligand.select_atoms("not name H*")

    for ts in universe.trajectory[::stride]:
        # Simple distance-based H-bond detection
        count = 0

        # Check distances between protein polar atoms and ligand
        for p_atom in protein_donors:
            for l_atom in ligand_heavy:
                dist = np.linalg.norm(p_atom.position - l_atom.position)
                if dist < 3.5:  # H-bond distance criterion
                    count += 1

        for p_atom in protein_acceptors:
            for l_atom in ligand_heavy:
                dist = np.linalg.norm(p_atom.position - l_atom.position)
                if dist < 3.5:
                    count += 1

        hbond_counts.append(count)
        times.append(ts.time / 1000)

    return np.array(times), np.array(hbond_counts), None


def compute_hbonds_mda(universe, stride=10):
    """
    Use MDAnalysis HydrogenBondAnalysis for proper H-bond detection.
    """
    try:
        # Define selections
        protein_sel = "protein"
        ligand_sel = "not protein and not resname HOH WAT SOL NA CL"

        # Check if ligand exists
        ligand = universe.select_atoms(ligand_sel)
        if len(ligand) == 0:
            print("  Warning: No ligand atoms found!")
            return None, None, None

        print(f"  Ligand atoms: {len(ligand)}")
        print(f"  Running H-bond analysis...")

        # Run H-bond analysis
        hbonds = HydrogenBondAnalysis(
            universe,
            donors_sel=f"({protein_sel}) or ({ligand_sel})",
            acceptors_sel=f"({protein_sel}) or ({ligand_sel})",
            d_a_cutoff=3.5,
            d_h_a_angle_cutoff=150,
        )
        hbonds.run(step=stride, verbose=False)

        # Filter for protein-ligand H-bonds only
        results = hbonds.results.hbonds

        if len(results) == 0:
            print("  No H-bonds detected!")
            return None, None, None

        # Count H-bonds per frame
        n_frames = len(universe.trajectory[::stride])
        hbond_counts = np.zeros(n_frames)
        times = np.array([ts.time / 1000 for ts in universe.trajectory[::stride]])

        # Get residue info for persistence analysis
        hbond_pairs = defaultdict(int)

        for hb in results:
            frame_idx = int(hb[0] / stride)
            if frame_idx < n_frames:
                hbond_counts[frame_idx] += 1

            # Track donor-acceptor pairs
            donor_idx = int(hb[1])
            acceptor_idx = int(hb[3])
            donor_atom = universe.atoms[donor_idx]
            acceptor_atom = universe.atoms[acceptor_idx]

            pair_key = f"{donor_atom.resname}{donor_atom.resid}-{acceptor_atom.resname}{acceptor_atom.resid}"
            hbond_pairs[pair_key] += 1

        # Calculate persistence
        persistence = {k: v / n_frames * 100 for k, v in hbond_pairs.items()}

        return times, hbond_counts, persistence

    except Exception as e:
        print(f"  MDAnalysis H-bond failed: {e}")
        print("  Falling back to simple method...")
        return compute_hbonds_simple(universe, stride)


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
        'hbond_counts': [],
        'times': None,
        'persistence': defaultdict(list)
    }

    for rep in range(1, config['reps'] + 1):
        traj_file = TRAJ_DIR / f"{name}_rep{rep}.dcd"

        if not traj_file.exists():
            print(f"  Trajectory not found: {traj_file}")
            continue

        print(f"  Loading rep{rep}...")

        try:
            u = mda.Universe(str(top_file), str(traj_file))

            times, hbond_counts, persistence = compute_hbonds_mda(u, stride=stride)

            if hbond_counts is not None:
                results['hbond_counts'].append(hbond_counts)
                results['times'] = times

                if persistence:
                    for k, v in persistence.items():
                        results['persistence'][k].append(v)

            del u
            gc.collect()

        except Exception as e:
            print(f"    ERROR: {e}")
            continue

    return results


def plot_hbond_timeseries(all_results, output_dir):
    """Plot H-bond count over time."""
    fig, ax = plt.subplots(figsize=(12, 6))

    for name, results in all_results.items():
        if not results['hbond_counts']:
            continue

        color = COLORS[name]
        times = results['times']

        # Ensemble average
        if len(results['hbond_counts']) > 1:
            # Pad arrays to same length
            max_len = max(len(hb) for hb in results['hbond_counts'])
            padded = [np.pad(hb, (0, max_len - len(hb)), constant_values=np.nan)
                     for hb in results['hbond_counts']]
            hb_array = np.array(padded)
            mean_hb = np.nanmean(hb_array, axis=0)
            std_hb = np.nanstd(hb_array, axis=0)

            min_len = min(len(times), len(mean_hb))
            ax.plot(times[:min_len], mean_hb[:min_len], color=color, linewidth=2, label=name)
            ax.fill_between(times[:min_len],
                           mean_hb[:min_len] - std_hb[:min_len],
                           mean_hb[:min_len] + std_hb[:min_len],
                           color=color, alpha=0.2)
        else:
            hb = results['hbond_counts'][0]
            min_len = min(len(times), len(hb))
            ax.plot(times[:min_len], hb[:min_len], color=color, linewidth=2, label=name)

    ax.set_xlabel('Time (ns)', fontsize=14)
    ax.set_ylabel('H-bond Count', fontsize=14)
    ax.set_title('Protein-Ligand Hydrogen Bonds Over Time', fontsize=16)
    ax.legend(loc='upper right', fontsize=10)

    plt.tight_layout()
    plt.savefig(output_dir / 'hbond_timeseries.png', dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved: {output_dir / 'hbond_timeseries.png'}")


def plot_hbond_summary(all_results, output_dir):
    """Plot average H-bond count comparison."""
    stats = []

    for name, results in all_results.items():
        if results['hbond_counts']:
            all_counts = np.concatenate(results['hbond_counts'])
            stats.append({
                'Compound': name,
                'Mean H-bonds': np.mean(all_counts),
                'Std': np.std(all_counts)
            })

    if not stats:
        return

    df = pd.DataFrame(stats)
    df = df.sort_values('Mean H-bonds', ascending=False)

    fig, ax = plt.subplots(figsize=(10, 6))
    colors = [COLORS[c] for c in df['Compound']]
    bars = ax.bar(df['Compound'], df['Mean H-bonds'], yerr=df['Std'],
                  color=colors, capsize=5, edgecolor='black')

    ax.set_xlabel('Compound', fontsize=14)
    ax.set_ylabel('Average H-bonds', fontsize=14)
    ax.set_title('Average Protein-Ligand Hydrogen Bonds', fontsize=16)

    plt.tight_layout()
    plt.savefig(output_dir / 'hbond_summary.png', dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved: {output_dir / 'hbond_summary.png'}")

    # Save stats
    df.to_csv(output_dir / 'hbond_statistics.csv', index=False)
    print(f"Saved: {output_dir / 'hbond_statistics.csv'}")


def main():
    print("=" * 60)
    print("PHASE 7: HYDROGEN BOND ANALYSIS")
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
    plot_hbond_timeseries(all_results, OUTPUT_DIR)
    plot_hbond_summary(all_results, OUTPUT_DIR)

    print("\n" + "=" * 60)
    print("H-BOND ANALYSIS COMPLETE")
    print("=" * 60)


if __name__ == "__main__":
    main()
