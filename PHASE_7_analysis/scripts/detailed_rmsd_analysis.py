#!/usr/bin/env python3
"""
Detailed RMSD Analysis per User Request:
1. Individual replicate RMSD plots for budesonide
2. Bound fraction calculation (% frames with RMSD < 3 Å)
3. Ligand COM distance to binding pocket center over time
"""

import MDAnalysis as mda
from MDAnalysis.analysis import align
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import warnings

warnings.filterwarnings('ignore')

TRAJ_DIR = Path(r'C:\Users\vasud\nod2-screening-data\PHASE_6\trajectories')
OUTPUT_DIR = Path(r'C:\Users\vasud\nod2-screening-data\PHASE_7_analysis\detailed_rmsd')
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

COMPOUNDS = {
    'budesonide': {'reps': 3, 'top': 'budesonide_rep1_solvated.pdb'},
    'febuxostat': {'reps': 3, 'top': 'febuxostat_rep1_solvated.pdb'},
    'ursodiol': {'reps': 3, 'top': 'ursodiol_rep3_solvated.pdb'},
    'natural_top': {'reps': 3, 'top': 'natural_top_rep1_solvated.pdb'},
    'decoy': {'reps': 1, 'top': 'decoy_rep1_solvated.pdb'},
}

WATER_RESNAMES = {"HOH", "WAT", "SOL", "TIP3", "TIP3P"}
ION_RESNAMES = {"NA", "CL", "K", "CA", "MG", "NA+", "CL-"}


def pick_ligand(universe):
    """Find the ligand (largest non-protein, non-water residue)."""
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
    sel = f"resid {lig_res.resid} and resname {lig_res.resname}"
    return sel, lig_res.atoms


def wrap_ligand_to_protein(lig_pos, prot_com, box):
    """Apply minimum image convention."""
    lig_com = lig_pos.mean(axis=0)
    delta = lig_com - prot_com
    shift = np.zeros(3)
    for i in range(3):
        if delta[i] > box[i] / 2:
            shift[i] = -box[i]
        elif delta[i] < -box[i] / 2:
            shift[i] = box[i]
    return lig_pos + shift


def compute_rmsd_and_com_distance(universe, stride=10):
    """
    Compute:
    1. Ligand RMSD (with proper backbone alignment)
    2. Ligand COM distance from initial position (binding pocket proxy)
    """
    lig_sel, lig_atoms = pick_ligand(universe)
    if lig_atoms is None:
        return None, None, None, None

    lig_heavy_sel = f"({lig_sel}) and not name H*"
    lig_heavy = universe.select_atoms(lig_heavy_sel)
    if len(lig_heavy) == 0:
        lig_heavy = lig_atoms

    protein_bb = universe.select_atoms("protein and backbone")

    # Get reference (frame 0)
    universe.trajectory[0]
    box = universe.dimensions[:3] if universe.dimensions is not None else np.array([100, 100, 100])
    ref_prot_pos = protein_bb.positions.copy()
    ref_prot_com = ref_prot_pos.mean(axis=0)
    ref_lig_pos = wrap_ligand_to_protein(lig_heavy.positions.copy(), ref_prot_com, box)

    # Binding pocket center = initial ligand COM position (after alignment)
    binding_pocket_center = ref_lig_pos.mean(axis=0)

    times = []
    rmsd_values = []
    com_distances = []

    for ts in universe.trajectory[::stride]:
        box = ts.dimensions[:3] if ts.dimensions is not None else np.array([100, 100, 100])

        cur_prot_pos = protein_bb.positions.copy()
        cur_prot_com = cur_prot_pos.mean(axis=0)
        cur_lig_pos = wrap_ligand_to_protein(lig_heavy.positions.copy(), cur_prot_com, box)

        # Kabsch alignment
        cur_prot_centered = cur_prot_pos - cur_prot_com
        ref_prot_centered = ref_prot_pos - ref_prot_com
        R, _ = align.rotation_matrix(cur_prot_centered, ref_prot_centered)

        # Transform ligand
        lig_centered = cur_lig_pos - cur_prot_com
        lig_rotated = np.dot(lig_centered, R.T)
        lig_aligned = lig_rotated + ref_prot_com

        # RMSD
        diff = lig_aligned - ref_lig_pos
        rmsd = np.sqrt((diff ** 2).sum() / len(lig_heavy))
        rmsd_values.append(rmsd)

        # COM distance from binding pocket
        cur_lig_com = lig_aligned.mean(axis=0)
        com_dist = np.linalg.norm(cur_lig_com - binding_pocket_center)
        com_distances.append(com_dist)

        times.append(ts.time / 1000)  # ps to ns

    return np.array(times), np.array(rmsd_values), np.array(com_distances), len(lig_heavy)


def main():
    print("=" * 70)
    print("DETAILED RMSD ANALYSIS")
    print("=" * 70)

    # =========================================================================
    # PART 1: Individual RMSD plots for budesonide
    # =========================================================================
    print("\n### PART 1: Individual Budesonide Replicates ###\n")

    fig, axes = plt.subplots(3, 1, figsize=(12, 10), sharex=True)

    budesonide_data = {}

    for rep in range(1, 4):
        print(f"Processing budesonide rep{rep}...")
        u = mda.Universe(
            str(TRAJ_DIR / 'budesonide_rep1_solvated.pdb'),
            str(TRAJ_DIR / f'budesonide_rep{rep}.dcd')
        )

        times, rmsd, com_dist, n_atoms = compute_rmsd_and_com_distance(u, stride=10)
        budesonide_data[rep] = {'times': times, 'rmsd': rmsd, 'com_dist': com_dist}

        ax = axes[rep - 1]
        ax.plot(times, rmsd, 'b-', linewidth=1.5, label=f'Rep{rep}')
        ax.axhline(y=3, color='green', linestyle='--', alpha=0.7, label='3Å threshold')
        ax.axhline(y=5, color='red', linestyle='--', alpha=0.7, label='5Å threshold')
        ax.set_ylabel('RMSD (Å)', fontsize=12)
        ax.set_title(f'Budesonide Replicate {rep}', fontsize=14)
        ax.legend(loc='upper right')
        ax.set_ylim(0, max(25, rmsd.max() * 1.1))
        ax.grid(True, alpha=0.3)

        # Statistics
        bound_frac = (rmsd < 3).sum() / len(rmsd) * 100
        print(f"  Rep{rep}: Mean={rmsd.mean():.2f}Å, Max={rmsd.max():.2f}Å, Bound(<3Å)={bound_frac:.1f}%")

    axes[2].set_xlabel('Time (ns)', fontsize=12)
    plt.tight_layout()
    plt.savefig(OUTPUT_DIR / 'budesonide_individual_rmsd.png', dpi=300)
    plt.close()
    print(f"\nSaved: {OUTPUT_DIR / 'budesonide_individual_rmsd.png'}")

    # =========================================================================
    # PART 2: Bound Fraction for ALL compounds
    # =========================================================================
    print("\n### PART 2: Bound Fraction (% frames with RMSD < 3Å) ###\n")

    print(f"{'Compound':<15} {'Rep':>4} {'Mean RMSD':>12} {'Max RMSD':>10} {'Bound Frac':>12}")
    print("-" * 60)

    all_bound_fractions = {}

    for name, config in COMPOUNDS.items():
        all_bound_fractions[name] = []

        for rep in range(1, config['reps'] + 1):
            traj_file = TRAJ_DIR / f"{name}_rep{rep}.dcd"
            if not traj_file.exists():
                continue

            u = mda.Universe(str(TRAJ_DIR / config['top']), str(traj_file))
            times, rmsd, com_dist, n_atoms = compute_rmsd_and_com_distance(u, stride=10)

            if rmsd is None:
                continue

            bound_frac = (rmsd < 3).sum() / len(rmsd) * 100
            all_bound_fractions[name].append(bound_frac)

            print(f"{name:<15} {rep:>4} {rmsd.mean():>12.2f} {rmsd.max():>10.2f} {bound_frac:>11.1f}%")

    print("\n" + "=" * 60)
    print("SUMMARY: Average Bound Fraction per Compound")
    print("=" * 60)
    for name, fracs in all_bound_fractions.items():
        if fracs:
            avg_frac = np.mean(fracs)
            print(f"{name:<15}: {avg_frac:>6.1f}% of frames with RMSD < 3Å")

    # =========================================================================
    # PART 3: Ligand COM Distance to Binding Pocket
    # =========================================================================
    print("\n### PART 3: Ligand COM Distance to Binding Pocket ###\n")

    fig, axes = plt.subplots(2, 3, figsize=(15, 8))
    axes = axes.flatten()

    plot_idx = 0
    for name, config in COMPOUNDS.items():
        ax = axes[plot_idx]

        for rep in range(1, config['reps'] + 1):
            traj_file = TRAJ_DIR / f"{name}_rep{rep}.dcd"
            if not traj_file.exists():
                continue

            u = mda.Universe(str(TRAJ_DIR / config['top']), str(traj_file))
            times, rmsd, com_dist, n_atoms = compute_rmsd_and_com_distance(u, stride=10)

            if com_dist is None:
                continue

            ax.plot(times, com_dist, linewidth=1, alpha=0.8, label=f'Rep{rep}')

        ax.set_title(f'{name}', fontsize=12)
        ax.set_xlabel('Time (ns)')
        ax.set_ylabel('COM Distance (Å)')
        ax.axhline(y=5, color='green', linestyle='--', alpha=0.5, label='5Å')
        ax.axhline(y=10, color='orange', linestyle='--', alpha=0.5, label='10Å')
        ax.legend(fontsize=8)
        ax.grid(True, alpha=0.3)

        plot_idx += 1

    # Hide unused subplot
    if plot_idx < len(axes):
        axes[plot_idx].set_visible(False)

    plt.suptitle('Ligand COM Distance from Initial Binding Position', fontsize=14)
    plt.tight_layout()
    plt.savefig(OUTPUT_DIR / 'ligand_com_distance.png', dpi=300)
    plt.close()
    print(f"Saved: {OUTPUT_DIR / 'ligand_com_distance.png'}")

    # =========================================================================
    # PART 4: Combined budesonide RMSD + COM distance
    # =========================================================================
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 8), sharex=True)

    colors = ['#1f77b4', '#ff7f0e', '#2ca02c']
    for rep in range(1, 4):
        data = budesonide_data[rep]
        ax1.plot(data['times'], data['rmsd'], color=colors[rep-1],
                 linewidth=1.5, label=f'Rep{rep}', alpha=0.8)
        ax2.plot(data['times'], data['com_dist'], color=colors[rep-1],
                 linewidth=1.5, label=f'Rep{rep}', alpha=0.8)

    ax1.axhline(y=3, color='green', linestyle='--', alpha=0.7)
    ax1.set_ylabel('Ligand RMSD (Å)', fontsize=12)
    ax1.set_title('Budesonide: RMSD and COM Distance', fontsize=14)
    ax1.legend()
    ax1.grid(True, alpha=0.3)

    ax2.axhline(y=5, color='green', linestyle='--', alpha=0.7)
    ax2.set_ylabel('COM Distance from\nBinding Pocket (Å)', fontsize=12)
    ax2.set_xlabel('Time (ns)', fontsize=12)
    ax2.legend()
    ax2.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(OUTPUT_DIR / 'budesonide_rmsd_and_com.png', dpi=300)
    plt.close()
    print(f"Saved: {OUTPUT_DIR / 'budesonide_rmsd_and_com.png'}")

    print("\n" + "=" * 70)
    print("DETAILED ANALYSIS COMPLETE")
    print("=" * 70)


if __name__ == "__main__":
    main()
