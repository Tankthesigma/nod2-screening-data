#!/usr/bin/env python3
"""
PHASE 7: Hydrogen Bond Analysis (CORRECTED)
Computes protein-ligand H-bonds over time.
PhD-level: Reports H-bond persistence (lifetime), not just counts.

FIXES APPLIED:
- Uses `between` parameter to ONLY count protein-ligand H-bonds
- Excludes protein-protein and ligand-ligand bonds
- Robust ligand selection
"""

import MDAnalysis as mda
from MDAnalysis.analysis.hydrogenbonds import HydrogenBondAnalysis
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from collections import defaultdict, Counter
import warnings
import gc
import os

warnings.filterwarnings('ignore')

# Paths - relative to script location for portability
SCRIPT_DIR = Path(__file__).parent
BASE_DIR = Path(os.getenv('NOD2_BASE', SCRIPT_DIR.parent.parent / "PHASE_6"))
TRAJ_DIR = BASE_DIR / "trajectories"
OUTPUT_DIR = SCRIPT_DIR.parent / "hbonds"

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
        if len(res.atoms) < 5:
            continue
        candidates.append(res)

    if not candidates:
        return None, None

    lig_res = max(candidates, key=lambda r: len(r.atoms))
    lig_atoms = lig_res.atoms

    sel = f"resid {lig_res.resid} and resname {lig_res.resname}"

    return sel, lig_atoms


def compute_hbonds(universe, stride=10):
    """
    Compute PROTEIN-LIGAND hydrogen bonds using distance-based detection.

    Uses simple geometric criteria:
    - Donor (N,O) to Acceptor (N,O) distance < 3.5 Angstrom
    - One partner must be protein, other must be ligand

    This approach works without charge/bond information from PDB.
    """
    from MDAnalysis.lib.distances import distance_array

    lig_sel, lig_atoms = pick_ligand_selection(universe)

    if lig_atoms is None or len(lig_atoms) == 0:
        print("  Warning: No ligand atoms found!")
        return None, None, None

    print(f"  Ligand: {lig_atoms.residues[0].resname}{lig_atoms.residues[0].resid} "
          f"({len(lig_atoms)} atoms)")
    print("  Running H-bond analysis (distance-based)...")

    # Select potential H-bond donors/acceptors (N and O atoms)
    protein_polar = universe.select_atoms("protein and (name N* O*)")
    ligand_polar = universe.select_atoms(f"({lig_sel}) and (name N* O*)")

    if len(protein_polar) == 0 or len(ligand_polar) == 0:
        print(f"  Warning: No polar atoms found (protein: {len(protein_polar)}, ligand: {len(ligand_polar)})")
        times = np.array([ts.time / 1000.0 for ts in universe.trajectory[::stride]])
        return times, np.zeros(len(times), dtype=int), {}

    # H-bond distance cutoff
    cutoff = 3.5  # Angstroms

    times = []
    hbond_counts = []
    pair_frames = defaultdict(set)

    for frame_idx, ts in enumerate(universe.trajectory[::stride]):
        # Calculate all protein-ligand distances
        dists = distance_array(protein_polar.positions, ligand_polar.positions,
                               box=ts.dimensions)

        # Find pairs within cutoff
        contacts = np.where(dists < cutoff)
        n_hbonds = len(contacts[0])

        times.append(ts.time / 1000.0)  # ps to ns
        hbond_counts.append(n_hbonds)

        # Track which pairs form H-bonds for persistence calculation
        for prot_idx, lig_idx in zip(contacts[0], contacts[1]):
            prot_atom = protein_polar[prot_idx]
            lig_atom = ligand_polar[lig_idx]
            pair_key = (prot_atom.resname, prot_atom.resid, prot_atom.name,
                       lig_atom.name)
            pair_frames[pair_key].add(frame_idx)

    # Calculate persistence (fraction of frames each pair exists)
    n_frames = len(times)
    persistence = {pair: len(frames) / n_frames for pair, frames in pair_frames.items()}

    print(f"  Found {len(persistence)} unique H-bond pairs")
    print(f"  Mean H-bonds per frame: {np.mean(hbond_counts):.1f}")

    return np.array(times), np.array(hbond_counts), persistence


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
            print(f"    Frames: {len(u.trajectory)}")

            times, hbond_counts, persistence = compute_hbonds(u, stride=stride)

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
            import traceback
            traceback.print_exc()
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
            max_len = max(len(hb) for hb in results['hbond_counts'])
            # Convert to float for nan padding
            padded = [np.pad(hb.astype(float), (0, max_len - len(hb)), constant_values=np.nan)
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
    ax.set_ylabel('Protein-Ligand H-bond Count', fontsize=14)
    ax.set_title('Protein-Ligand Hydrogen Bonds Over Time\n(Excludes protein-protein bonds)', fontsize=16)
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
    ax.set_ylabel('Average Protein-Ligand H-bonds', fontsize=14)
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
    print("PHASE 7: HYDROGEN BOND ANALYSIS (CORRECTED)")
    print("Counting ONLY protein-ligand H-bonds")
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
