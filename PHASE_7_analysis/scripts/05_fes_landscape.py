#!/usr/bin/env python3
"""
PHASE 7: Free Energy Surface (FES) Proxy (CORRECTED)
Creates 2D density plot of Ligand RMSD vs Radius of Gyration.
PhD-level: Shows conformational stability landscape.

FIXES APPLIED:
- Proper backbone-aligned ligand RMSD using Kabsch algorithm
- Ligand Rg (measures ligand compactness/conformational state)
- Robust ligand selection
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
import os

warnings.filterwarnings('ignore')

# Paths - relative to script location for portability
SCRIPT_DIR = Path(__file__).parent
BASE_DIR = Path(os.getenv('NOD2_BASE', SCRIPT_DIR.parent.parent / "PHASE_6"))
TRAJ_DIR = BASE_DIR / "trajectories"
OUTPUT_DIR = SCRIPT_DIR.parent / "fes"

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


def compute_ligand_rmsd_rg(universe, stride=10):
    """
    Compute ligand RMSD (properly aligned) and ligand radius of gyration.

    Method:
    1. Align protein backbone to reference (Kabsch algorithm)
    2. Measure ligand RMSD in aligned frame
    3. Compute ligand Rg (intrinsic property, doesn't need alignment)
    """
    lig_sel, lig_atoms = pick_ligand_selection(universe)

    if lig_atoms is None or len(lig_atoms) == 0:
        print("  Warning: No ligand atoms found!")
        return None, None

    # Get heavy atoms for RMSD
    lig_heavy_sel = f"({lig_sel}) and not name H*"
    lig_heavy = universe.select_atoms(lig_heavy_sel)

    if len(lig_heavy) == 0:
        lig_heavy_sel = lig_sel
        lig_heavy = lig_atoms

    print(f"  Ligand: {lig_atoms.residues[0].resname}{lig_atoms.residues[0].resid} "
          f"({len(lig_heavy)} heavy atoms)")

    # Use MDAnalysis RMSD with proper alignment
    # Align on protein backbone, measure ligand RMSD
    R = rms.RMSD(
        universe,
        universe,
        select="protein and backbone",
        groupselections=[lig_heavy_sel],
        ref_frame=0
    )
    R.run(step=stride, verbose=False)

    # Extract ligand RMSD (column 3)
    rmsd_values = R.results.rmsd[:, 3]

    # Compute ligand Rg for each frame
    # Rg is an intrinsic property - doesn't depend on alignment
    rg_values = []
    ligand = universe.select_atoms(lig_sel)

    frame_idx = 0
    for ts in universe.trajectory[::stride]:
        # Radius of gyration
        rg = ligand.radius_of_gyration()
        rg_values.append(rg)
        frame_idx += 1

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
            print(f"    Frames: {len(u.trajectory)}")

            rmsd, rg = compute_ligand_rmsd_rg(u, stride=stride)

            if rmsd is not None:
                results['rmsd'].extend(rmsd)
                results['rg'].extend(rg)
                print(f"    Mean RMSD: {np.mean(rmsd):.2f} Ang, Mean Rg: {np.mean(rg):.2f} Ang")

            del u
            gc.collect()

        except Exception as e:
            print(f"    ERROR: {e}")
            import traceback
            traceback.print_exc()
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

        ax.axvline(x=3, color='orange', linestyle='--', alpha=0.7, label='Stable (3 Ang)')
        ax.axvline(x=5, color='red', linestyle='--', alpha=0.7, label='Unstable (5 Ang)')

        ax.set_xlabel('Ligand RMSD (Ang)')
        ax.set_ylabel('Ligand Rg (Ang)')
        ax.set_title(f'{name}', fontsize=14, fontweight='bold')
        ax.legend(loc='upper right', fontsize=8)

    # Hide unused axes
    for idx in range(len(all_results), 6):
        axes[idx].set_visible(False)

    plt.suptitle('Free Energy Surface Proxy: RMSD vs Rg\n(Backbone-aligned RMSD)', fontsize=16, fontweight='bold')
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

    ax.set_xlabel('Ligand RMSD (Ang)', fontsize=14)
    ax.set_ylabel('Ligand Rg (Ang)', fontsize=14)
    ax.set_title('Conformational Landscape: All Compounds\n(X = centroid)', fontsize=16)
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
                'Mean RMSD (Ang)': f"{np.mean(rmsd):.2f}",
                'Std RMSD (Ang)': f"{np.std(rmsd):.2f}",
                'Mean Rg (Ang)': f"{np.mean(rg):.2f}",
                'Std Rg (Ang)': f"{np.std(rg):.2f}",
                '% RMSD < 3Ang': f"{np.mean(rmsd < 3) * 100:.1f}",
            })

    df = pd.DataFrame(stats)
    df.to_csv(output_dir / 'fes_statistics.csv', index=False)
    print(f"\nSaved: {output_dir / 'fes_statistics.csv'}")
    print("\n" + df.to_string(index=False))


def main():
    print("=" * 60)
    print("PHASE 7: FREE ENERGY SURFACE ANALYSIS (CORRECTED)")
    print("Using backbone-aligned ligand RMSD")
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
