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
    Compute PROTEIN-LIGAND hydrogen bonds ONLY.

    Uses the `between` parameter to ensure we only count:
    - Protein donor -> Ligand acceptor
    - Ligand donor -> Protein acceptor

    Does NOT count:
    - Protein-protein H-bonds (structural)
    - Ligand-ligand H-bonds (intramolecular)
    - Water H-bonds
    """
    protein_sel = "protein"
    lig_sel, lig_atoms = pick_ligand_selection(universe)

    if lig_atoms is None or len(lig_atoms) == 0:
        print("  Warning: No ligand atoms found!")
        return None, None, None

    print(f"  Ligand: {lig_atoms.residues[0].resname}{lig_atoms.residues[0].resid} "
          f"({len(lig_atoms)} atoms)")
    print("  Running H-bond analysis (protein-ligand ONLY)...")

    # Use the `between` parameter - this is the KEY FIX
    # between=[sel1, sel2] finds H-bonds where one partner is in sel1 and other in sel2
    try:
        hbonds = HydrogenBondAnalysis(
            universe,
            between=[protein_sel, lig_sel],
            d_a_cutoff=3.5,
            d_h_a_angle_cutoff=150,
            update_selections=True
        )
        hbonds.run(step=stride, verbose=False)
    except TypeError:
        # Older MDAnalysis versions may not support `between`
        # Fall back to manual filtering
        print("  Note: Using fallback H-bond detection (older MDAnalysis)")
        return compute_hbonds_fallback(universe, lig_sel, stride)

    results = hbonds.results.hbonds

    # Build time axis
    times = np.array([ts.time / 1000.0 for ts in universe.trajectory[::stride]])
    n_frames = len(times)

    if results is None or len(results) == 0:
        print("  No H-bonds detected!")
        return times, np.zeros(n_frames, dtype=int), {}

    # Map frame indices
    sampled_frames = np.arange(0, len(universe.trajectory), stride, dtype=int)
    frame_to_idx = {f: i for i, f in enumerate(sampled_frames)}

    # Count H-bonds per frame
    counts = np.zeros(n_frames, dtype=int)
    pair_frames = defaultdict(set)

    for row in results:
        frame = int(row[0])
        if frame not in frame_to_idx:
            continue

        idx = frame_to_idx[frame]
        counts[idx] += 1

        # Track donor-acceptor pairs for persistence
        donor_idx = int(row[1])
        acceptor_idx = int(row[3])

        d = universe.atoms[donor_idx]
        a = universe.atoms[acceptor_idx]

        pair_key = f"{d.resname}{d.resid}:{d.name}->{a.resname}{a.resid}:{a.name}"
        pair_frames[pair_key].add(idx)

    # Calculate persistence (% of frames each H-bond exists)
    persistence = {k: len(v) / n_frames * 100.0 for k, v in pair_frames.items()}

    print(f"  Found {len(persistence)} unique H-bond pairs")
    print(f"  Mean H-bonds per frame: {np.mean(counts):.1f}")

    return times, counts, persistence


def compute_hbonds_fallback(universe, lig_sel, stride=10):
    """
    Fallback H-bond detection for older MDAnalysis versions.
    Manually filters for protein-ligand bonds.
    """
    protein = universe.select_atoms("protein")
    ligand = universe.select_atoms(lig_sel)

    # Run H-bond analysis on combined selection
    combined_sel = f"protein or ({lig_sel})"

    hbonds = HydrogenBondAnalysis(
        universe,
        donors_sel=combined_sel,
        acceptors_sel=combined_sel,
        d_a_cutoff=3.5,
        d_h_a_angle_cutoff=150,
    )
    hbonds.run(step=stride, verbose=False)

    results = hbonds.results.hbonds
    times = np.array([ts.time / 1000.0 for ts in universe.trajectory[::stride]])
    n_frames = len(times)

    if results is None or len(results) == 0:
        return times, np.zeros(n_frames, dtype=int), {}

    sampled_frames = np.arange(0, len(universe.trajectory), stride, dtype=int)
    frame_to_idx = {f: i for i, f in enumerate(sampled_frames)}

    counts = np.zeros(n_frames, dtype=int)
    pair_frames = defaultdict(set)

    # Get atom indices for filtering
    protein_indices = set(protein.indices)
    ligand_indices = set(ligand.indices)

    for row in results:
        frame = int(row[0])
        donor_idx = int(row[1])
        acceptor_idx = int(row[3])

        if frame not in frame_to_idx:
            continue

        # FILTER: One must be protein, other must be ligand (XOR logic)
        donor_is_protein = donor_idx in protein_indices
        donor_is_ligand = donor_idx in ligand_indices
        acceptor_is_protein = acceptor_idx in protein_indices
        acceptor_is_ligand = acceptor_idx in ligand_indices

        is_protein_ligand = (
            (donor_is_protein and acceptor_is_ligand) or
            (donor_is_ligand and acceptor_is_protein)
        )

        if not is_protein_ligand:
            continue  # Skip protein-protein or ligand-ligand bonds

        idx = frame_to_idx[frame]
        counts[idx] += 1

        d = universe.atoms[donor_idx]
        a = universe.atoms[acceptor_idx]
        pair_key = f"{d.resname}{d.resid}:{d.name}->{a.resname}{a.resid}:{a.name}"
        pair_frames[pair_key].add(idx)

    persistence = {k: len(v) / n_frames * 100.0 for k, v in pair_frames.items()}

    return times, counts, persistence


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
