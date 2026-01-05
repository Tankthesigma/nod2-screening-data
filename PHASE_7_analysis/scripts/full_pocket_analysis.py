#!/usr/bin/env python3
"""
FULL POCKET OCCUPANCY ANALYSIS - ALL COMPOUNDS, ALL REPLICATES
"""

import MDAnalysis as mda
from MDAnalysis.lib.distances import distance_array
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from pathlib import Path
import warnings

warnings.filterwarnings('ignore')

TRAJ_DIR = Path(r'C:\Users\vasud\nod2-screening-data\PHASE_6\trajectories')
OUTPUT_DIR = Path(r'C:\Users\vasud\nod2-screening-data\PHASE_7_analysis\pocket_analysis')
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

COMPOUNDS = {
    'budesonide': {'type': 'Positive Control', 'reps': 3, 'top': 'budesonide_rep1_solvated.pdb', 'color': '#1f77b4'},
    'febuxostat': {'type': 'Lead Candidate', 'reps': 3, 'top': 'febuxostat_rep1_solvated.pdb', 'color': '#2ca02c'},
    'ursodiol': {'type': 'Secondary', 'reps': 3, 'top': 'ursodiol_rep3_solvated.pdb', 'color': '#ff7f0e'},
    'natural_top': {'type': 'Natural Product', 'reps': 3, 'top': 'natural_top_rep1_solvated.pdb', 'color': '#9467bd'},
    'decoy': {'type': 'Negative Control', 'reps': 1, 'top': 'decoy_rep1_solvated.pdb', 'color': '#d62728'},
}

WATER_RESNAMES = {"HOH", "WAT", "SOL", "TIP3", "TIP3P"}
ION_RESNAMES = {"NA", "CL", "K", "CA", "MG", "NA+", "CL-"}


def pick_ligand(universe):
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
    lig_com = lig_pos.mean(axis=0)
    delta = lig_com - prot_com
    shift = np.zeros(3)
    for i in range(3):
        if delta[i] > box[i] / 2:
            shift[i] = -box[i]
        elif delta[i] < -box[i] / 2:
            shift[i] = box[i]
    return lig_pos + shift


def get_pocket_residues(universe, lig_sel, cutoff=5.0):
    universe.trajectory[0]
    lig_heavy = universe.select_atoms(f"({lig_sel}) and not name H*")
    protein_heavy = universe.select_atoms("protein and not name H*")
    box = universe.dimensions
    dists = distance_array(lig_heavy.positions, protein_heavy.positions, box=box)
    contact_mask = np.any(dists < cutoff, axis=0)
    contact_atoms = protein_heavy[contact_mask]
    pocket_resids = set(contact_atoms.resids)
    return pocket_resids


def compute_full_metrics(universe, lig_sel, pocket_resids, stride=10):
    lig_heavy = universe.select_atoms(f"({lig_sel}) and not name H*")
    if len(lig_heavy) == 0:
        lig_heavy = universe.select_atoms(lig_sel)

    pocket_sel = " or ".join([f"resid {r}" for r in pocket_resids])
    pocket_atoms = universe.select_atoms(f"protein and ({pocket_sel}) and not name H*")
    protein_bb = universe.select_atoms("protein and backbone")

    times = []
    n_contacts = []
    min_distances = []

    for ts in universe.trajectory[::stride]:
        box = ts.dimensions
        cur_prot_com = protein_bb.positions.mean(axis=0)
        lig_pos = wrap_ligand_to_protein(lig_heavy.positions.copy(), cur_prot_com, box)
        pocket_pos = pocket_atoms.positions
        dists = distance_array(lig_pos, pocket_pos, box=box)

        min_dist = dists.min()
        min_distances.append(min_dist)

        contact_mask = np.any(dists < 5.0, axis=0)
        contact_atoms_frame = pocket_atoms[contact_mask]
        contacting_resids = set(contact_atoms_frame.resids)
        n_contacts.append(len(contacting_resids))
        times.append(ts.time / 1000)

    pocket_occ_5 = sum(1 for n in n_contacts if n >= 5) / len(n_contacts) * 100
    pocket_occ_3 = sum(1 for n in n_contacts if n >= 3) / len(n_contacts) * 100
    pocket_occ_1 = sum(1 for n in n_contacts if n >= 1) / len(n_contacts) * 100

    return {
        'times': np.array(times),
        'n_contacts': np.array(n_contacts),
        'min_distances': np.array(min_distances),
        'pocket_occ_5': pocket_occ_5,
        'pocket_occ_3': pocket_occ_3,
        'pocket_occ_1': pocket_occ_1,
        'mean_contacts': np.mean(n_contacts),
        'mean_min_dist': np.mean(min_distances),
        'max_contacts': np.max(n_contacts),
    }


def main():
    print("=" * 80)
    print("FULL POCKET OCCUPANCY ANALYSIS - ALL COMPOUNDS, ALL REPLICATES")
    print("=" * 80)

    # Define pocket from budesonide frame 0
    u_bud = mda.Universe(
        str(TRAJ_DIR / 'budesonide_rep1_solvated.pdb'),
        str(TRAJ_DIR / 'budesonide_rep1.dcd')
    )
    lig_sel_bud, _ = pick_ligand(u_bud)
    reference_pocket_resids = get_pocket_residues(u_bud, lig_sel_bud, cutoff=5.0)

    print(f"\nReference pocket: {len(reference_pocket_resids)} residues from budesonide frame 0")
    print(f"Residues: {sorted(reference_pocket_resids)}\n")

    # Collect all data
    all_data = []

    print("=" * 100)
    print(f"{'Compound':<12} {'Type':<18} {'Rep':>4} {'Occ>=5':>8} {'Occ>=3':>8} {'Occ>=1':>8} {'MeanCont':>9} {'MinDist':>9} {'MaxCont':>8}")
    print("=" * 100)

    for name, config in COMPOUNDS.items():
        for rep in range(1, config['reps'] + 1):
            traj_file = TRAJ_DIR / f"{name}_rep{rep}.dcd"
            if not traj_file.exists():
                continue

            u = mda.Universe(str(TRAJ_DIR / config['top']), str(traj_file))
            lig_sel, _ = pick_ligand(u)

            if lig_sel is None:
                continue

            metrics = compute_full_metrics(u, lig_sel, reference_pocket_resids, stride=10)

            row = {
                'Compound': name,
                'Type': config['type'],
                'Rep': rep,
                'Pocket_Occ_5': metrics['pocket_occ_5'],
                'Pocket_Occ_3': metrics['pocket_occ_3'],
                'Pocket_Occ_1': metrics['pocket_occ_1'],
                'Mean_Contacts': metrics['mean_contacts'],
                'Mean_Min_Dist': metrics['mean_min_dist'],
                'Max_Contacts': metrics['max_contacts'],
                'times': metrics['times'],
                'n_contacts': metrics['n_contacts'],
                'min_distances': metrics['min_distances'],
                'color': config['color'],
            }
            all_data.append(row)

            print(f"{name:<12} {config['type']:<18} {rep:>4} {metrics['pocket_occ_5']:>7.1f}% {metrics['pocket_occ_3']:>7.1f}% {metrics['pocket_occ_1']:>7.1f}% {metrics['mean_contacts']:>9.1f} {metrics['mean_min_dist']:>8.2f}Å {metrics['max_contacts']:>8}")

    # Create summary by compound
    print("\n" + "=" * 80)
    print("SUMMARY BY COMPOUND (Averaged across replicates)")
    print("=" * 80)

    df = pd.DataFrame([{k: v for k, v in d.items() if k not in ['times', 'n_contacts', 'min_distances', 'color']} for d in all_data])

    summary = df.groupby(['Compound', 'Type']).agg({
        'Pocket_Occ_5': ['mean', 'std', 'min', 'max'],
        'Mean_Contacts': 'mean',
        'Mean_Min_Dist': 'mean',
    }).round(1)

    print("\n" + summary.to_string())

    # Save to CSV
    df_save = df.drop(columns=['Type'])
    df_save.to_csv(OUTPUT_DIR / 'full_pocket_metrics.csv', index=False)
    print(f"\nSaved: {OUTPUT_DIR / 'full_pocket_metrics.csv'}")

    # =========================================================================
    # MEGA VISUALIZATION
    # =========================================================================
    fig, axes = plt.subplots(5, 3, figsize=(18, 20))

    row_idx = 0
    for name, config in COMPOUNDS.items():
        compound_data = [d for d in all_data if d['Compound'] == name]

        for rep_idx, data in enumerate(compound_data):
            if row_idx >= 5:
                break
            ax = axes[row_idx, rep_idx]

            # Plot contacts over time
            ax.fill_between(data['times'], 0, data['n_contacts'],
                           where=np.array(data['n_contacts']) >= 5,
                           color='green', alpha=0.4, label='>=5 contacts')
            ax.fill_between(data['times'], 0, data['n_contacts'],
                           where=(np.array(data['n_contacts']) >= 1) & (np.array(data['n_contacts']) < 5),
                           color='yellow', alpha=0.4, label='1-4 contacts')
            ax.plot(data['times'], data['n_contacts'], color=config['color'], linewidth=1)
            ax.axhline(y=5, color='green', linestyle='--', alpha=0.7)
            ax.set_ylim(0, 10)
            ax.set_title(f"{name} Rep{data['Rep']}\nPocket Occ: {data['Pocket_Occ_5']:.1f}%", fontsize=10)
            ax.set_xlabel('Time (ns)')
            ax.set_ylabel('# Contacts')
            ax.grid(True, alpha=0.3)

        # Fill empty slots
        for rep_idx in range(len(compound_data), 3):
            axes[row_idx, rep_idx].set_visible(False)

        row_idx += 1

    plt.suptitle('POCKET OCCUPANCY: All Compounds, All Replicates\n(Green = 5+ contacts = IN POCKET)', fontsize=14, fontweight='bold')
    plt.tight_layout()
    plt.savefig(OUTPUT_DIR / 'ALL_pocket_occupancy.png', dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved: {OUTPUT_DIR / 'ALL_pocket_occupancy.png'}")

    # =========================================================================
    # BAR CHART COMPARISON
    # =========================================================================
    fig, ax = plt.subplots(figsize=(14, 8))

    compounds = []
    occ_values = []
    colors = []
    labels = []

    for name, config in COMPOUNDS.items():
        compound_data = [d for d in all_data if d['Compound'] == name]
        for data in compound_data:
            compounds.append(f"{name}\nRep{data['Rep']}")
            occ_values.append(data['Pocket_Occ_5'])
            colors.append(config['color'])
            labels.append(f"{data['Pocket_Occ_5']:.0f}%")

    bars = ax.bar(compounds, occ_values, color=colors, edgecolor='black', linewidth=1)

    # Add value labels on bars
    for bar, label in zip(bars, labels):
        height = bar.get_height()
        ax.annotate(label,
                    xy=(bar.get_x() + bar.get_width() / 2, height),
                    xytext=(0, 3),
                    textcoords="offset points",
                    ha='center', va='bottom', fontsize=10, fontweight='bold')

    ax.axhline(y=50, color='orange', linestyle='--', linewidth=2, label='50% threshold')
    ax.axhline(y=80, color='green', linestyle='--', linewidth=2, label='80% threshold (stable)')
    ax.set_ylabel('Pocket Occupancy (% frames with ≥5 contacts)', fontsize=12)
    ax.set_title('POCKET OCCUPANCY COMPARISON\nAll Compounds, All Replicates', fontsize=14, fontweight='bold')
    ax.legend(loc='upper right')
    ax.set_ylim(0, 110)
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    plt.savefig(OUTPUT_DIR / 'pocket_occupancy_barplot.png', dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved: {OUTPUT_DIR / 'pocket_occupancy_barplot.png'}")

    # =========================================================================
    # FINAL VERDICT TABLE
    # =========================================================================
    print("\n" + "=" * 80)
    print("*** FINAL VERDICT: POCKET OCCUPANCY RANKING ***")
    print("=" * 80)

    # Sort by pocket occupancy
    sorted_data = sorted(all_data, key=lambda x: x['Pocket_Occ_5'], reverse=True)

    print(f"\n{'Rank':<5} {'Compound':<15} {'Rep':<5} {'Pocket Occ':<12} {'Status':<15}")
    print("-" * 55)

    for i, d in enumerate(sorted_data, 1):
        occ = d['Pocket_Occ_5']
        if occ >= 90:
            status = "[EXCELLENT]"
        elif occ >= 70:
            status = "[GOOD]"
        elif occ >= 50:
            status = "[MODERATE]"
        elif occ >= 20:
            status = "[WEAK]"
        else:
            status = "[UNBOUND]"

        print(f"{i:<5} {d['Compound']:<15} {d['Rep']:<5} {occ:>10.1f}%   {status}")

    print("\n" + "=" * 80)
    print("INTERPRETATION:")
    print("  >=90% = Excellent binding (stays in pocket almost entire simulation)")
    print("  70-90% = Good binding")
    print("  50-70% = Moderate (explores but returns)")
    print("  20-50% = Weak binding")
    print("  <20% = Essentially unbound")
    print("=" * 80)


if __name__ == "__main__":
    main()
