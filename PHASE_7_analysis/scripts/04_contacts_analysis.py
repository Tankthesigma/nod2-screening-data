#!/usr/bin/env python3
"""
PHASE 7: Contact Analysis
Identifies protein residues in contact with ligand.
PhD-level: Reports contact frequency across trajectory.
"""

import MDAnalysis as mda
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
OUTPUT_DIR = Path(r"C:\Users\vasud\nod2-screening-data\PHASE_7_analysis\contacts")

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

# Contact distance threshold
CONTACT_CUTOFF = 4.0  # Angstroms

sns.set_context("paper", font_scale=1.5)
sns.set_style("whitegrid")


def compute_contacts(universe, contact_cutoff=4.0, stride=10):
    """
    Compute residue-ligand contacts over trajectory.
    A contact is defined as any heavy atom within cutoff distance.
    """
    protein = universe.select_atoms("protein")
    ligand = universe.select_atoms("not protein and not resname HOH WAT SOL NA CL")

    if len(ligand) == 0:
        print("  Warning: No ligand atoms found!")
        return None, 0

    ligand_heavy = ligand.select_atoms("not name H*")
    print(f"  Ligand heavy atoms: {len(ligand_heavy)}")

    # Get unique residues
    protein_residues = protein.residues

    # Track contacts per residue
    contact_counts = defaultdict(int)
    n_frames = 0

    for ts in universe.trajectory[::stride]:
        n_frames += 1
        ligand_positions = ligand_heavy.positions

        for res in protein_residues:
            res_heavy = res.atoms.select_atoms("not name H*")
            if len(res_heavy) == 0:
                continue

            # Check for any contact
            for res_atom in res_heavy:
                min_dist = np.min(np.linalg.norm(ligand_positions - res_atom.position, axis=1))
                if min_dist < contact_cutoff:
                    contact_counts[f"{res.resname}{res.resid}"] += 1
                    break  # Count residue once per frame

    # Convert to frequency
    contact_freq = {k: v / n_frames * 100 for k, v in contact_counts.items()}

    return contact_freq, n_frames


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
        'contact_freq': defaultdict(list)
    }

    for rep in range(1, config['reps'] + 1):
        traj_file = TRAJ_DIR / f"{name}_rep{rep}.dcd"

        if not traj_file.exists():
            print(f"  Trajectory not found: {traj_file}")
            continue

        print(f"  Loading rep{rep}...")

        try:
            u = mda.Universe(str(top_file), str(traj_file))

            contact_freq, n_frames = compute_contacts(u, CONTACT_CUTOFF, stride=stride)

            if contact_freq:
                print(f"    Found {len(contact_freq)} contacting residues")
                for res, freq in contact_freq.items():
                    results['contact_freq'][res].append(freq)

            del u
            gc.collect()

        except Exception as e:
            print(f"    ERROR: {e}")
            continue

    return results


def plot_contact_heatmap(all_results, output_dir, min_freq=30):
    """Create contact frequency heatmap across compounds."""

    # Collect all contacting residues with >min_freq in any compound
    all_residues = set()
    for name, results in all_results.items():
        for res, freqs in results['contact_freq'].items():
            mean_freq = np.mean(freqs)
            if mean_freq >= min_freq:
                all_residues.add(res)

    if not all_residues:
        print("No residues with significant contact frequency found!")
        return

    # Sort residues by residue number
    def get_resid(res_str):
        import re
        match = re.search(r'\d+', res_str)
        return int(match.group()) if match else 0

    sorted_residues = sorted(all_residues, key=get_resid)

    # Build data matrix
    data = {}
    for name, results in all_results.items():
        data[name] = []
        for res in sorted_residues:
            if res in results['contact_freq']:
                mean_freq = np.mean(results['contact_freq'][res])
            else:
                mean_freq = 0
            data[name].append(mean_freq)

    df = pd.DataFrame(data, index=sorted_residues)

    # Plot heatmap
    fig, ax = plt.subplots(figsize=(12, max(8, len(sorted_residues) * 0.3)))

    sns.heatmap(df, cmap='YlOrRd', annot=True, fmt='.0f', ax=ax,
                cbar_kws={'label': 'Contact Frequency (%)'})

    ax.set_xlabel('Compound', fontsize=14)
    ax.set_ylabel('Residue', fontsize=14)
    ax.set_title(f'Protein-Ligand Contacts (>{min_freq}% frequency)', fontsize=16)

    plt.tight_layout()
    plt.savefig(output_dir / 'contact_heatmap.png', dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved: {output_dir / 'contact_heatmap.png'}")

    # Save data
    df.to_csv(output_dir / 'contact_frequencies.csv')
    print(f"Saved: {output_dir / 'contact_frequencies.csv'}")


def plot_contact_comparison(all_results, output_dir):
    """Compare number of stable contacts per compound."""

    stats = []
    for name, results in all_results.items():
        if results['contact_freq']:
            # Count residues with >50% contact
            stable_contacts = sum(1 for freqs in results['contact_freq'].values()
                                 if np.mean(freqs) > 50)
            total_contacts = len(results['contact_freq'])

            stats.append({
                'Compound': name,
                'Stable Contacts (>50%)': stable_contacts,
                'Total Contacts': total_contacts
            })

    if not stats:
        return

    df = pd.DataFrame(stats)
    df = df.sort_values('Stable Contacts (>50%)', ascending=False)

    fig, ax = plt.subplots(figsize=(10, 6))
    x = np.arange(len(df))
    width = 0.35

    colors = [COLORS[c] for c in df['Compound']]
    ax.bar(x - width/2, df['Stable Contacts (>50%)'], width, label='Stable (>50%)',
           color=colors, edgecolor='black')
    ax.bar(x + width/2, df['Total Contacts'], width, label='Total',
           color=colors, alpha=0.4, edgecolor='black')

    ax.set_xlabel('Compound', fontsize=14)
    ax.set_ylabel('Number of Residues', fontsize=14)
    ax.set_title('Protein-Ligand Contact Comparison', fontsize=16)
    ax.set_xticks(x)
    ax.set_xticklabels(df['Compound'])
    ax.legend()

    plt.tight_layout()
    plt.savefig(output_dir / 'contact_comparison.png', dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved: {output_dir / 'contact_comparison.png'}")

    # Save stats
    df.to_csv(output_dir / 'contact_statistics.csv', index=False)
    print(f"Saved: {output_dir / 'contact_statistics.csv'}")


def identify_key_residues(all_results, output_dir):
    """Identify residues that are common binding anchors."""

    # Find residues that contact multiple compounds with >50% freq
    residue_compounds = defaultdict(list)

    for name, results in all_results.items():
        for res, freqs in results['contact_freq'].items():
            if np.mean(freqs) > 50:
                residue_compounds[res].append((name, np.mean(freqs)))

    # Sort by number of compounds contacting
    key_residues = [(res, contacts) for res, contacts in residue_compounds.items()
                    if len(contacts) >= 2]
    key_residues.sort(key=lambda x: len(x[1]), reverse=True)

    print("\n" + "=" * 50)
    print("KEY BINDING RESIDUES (contact >50% in multiple compounds)")
    print("=" * 50)

    data = []
    for res, contacts in key_residues:
        compounds = ", ".join([f"{c[0]}({c[1]:.0f}%)" for c in contacts])
        print(f"  {res}: {compounds}")
        data.append({'Residue': res, 'Compounds': len(contacts), 'Details': compounds})

    if data:
        df = pd.DataFrame(data)
        df.to_csv(output_dir / 'key_binding_residues.csv', index=False)
        print(f"\nSaved: {output_dir / 'key_binding_residues.csv'}")


def main():
    print("=" * 60)
    print("PHASE 7: CONTACT ANALYSIS")
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
    plot_contact_heatmap(all_results, OUTPUT_DIR, min_freq=30)
    plot_contact_comparison(all_results, OUTPUT_DIR)
    identify_key_residues(all_results, OUTPUT_DIR)

    print("\n" + "=" * 60)
    print("CONTACT ANALYSIS COMPLETE")
    print("=" * 60)


if __name__ == "__main__":
    main()
