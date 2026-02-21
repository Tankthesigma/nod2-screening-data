#!/usr/bin/env python
"""
RMSF Analysis - Compare protein flexibility with/without drug
Tests chaperone hypothesis: Does drug stabilize mutant NOD2?
"""

import os
import sys
import warnings
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from datetime import datetime
warnings.filterwarnings('ignore')

import MDAnalysis as mda
from MDAnalysis.analysis import rms, align


def calculate_rmsf(pdb_path, dcd_path, selection='protein and name CA',
                   start_frac=0.05, end_frac=0.25):
    """
    Calculate RMSF for C-alpha atoms with proper alignment
    Returns: residue numbers, RMSF values
    """
    print(f"  Loading: {os.path.basename(dcd_path)}")
    u = mda.Universe(pdb_path, dcd_path)

    total_frames = len(u.trajectory)
    start_frame = int(total_frames * start_frac)
    end_frame = int(total_frames * end_frac)

    print(f"  Frames: {start_frame}-{end_frame} of {total_frames}")

    # Select C-alpha atoms
    ca_atoms = u.select_atoms(selection)
    n_atoms = len(ca_atoms)
    print(f"  C-alpha atoms: {n_atoms}")

    # Get residue info
    resids = ca_atoms.resids
    resnames = ca_atoms.resnames

    # Use MDAnalysis RMSF with alignment
    from MDAnalysis.analysis.rms import RMSF

    # First align trajectory to reference (first frame in analysis window)
    u.trajectory[start_frame]
    ref_positions = ca_atoms.positions.copy()

    # Collect aligned positions
    positions = []

    for i, ts in enumerate(u.trajectory):
        if i < start_frame:
            continue
        if i > end_frame:
            break

        # Align current frame to reference using Kabsch algorithm
        current_pos = ca_atoms.positions.copy()

        # Center both
        ref_centered = ref_positions - ref_positions.mean(axis=0)
        cur_centered = current_pos - current_pos.mean(axis=0)

        # Compute optimal rotation (Kabsch)
        H = cur_centered.T @ ref_centered
        U, S, Vt = np.linalg.svd(H)
        R = Vt.T @ U.T

        # Handle reflection
        if np.linalg.det(R) < 0:
            Vt[-1, :] *= -1
            R = Vt.T @ U.T

        # Apply rotation and store
        aligned_pos = (cur_centered @ R) + ref_positions.mean(axis=0)
        positions.append(aligned_pos)

    positions = np.array(positions)
    n_frames = len(positions)
    print(f"  Collected and aligned {n_frames} frames")

    # Calculate mean structure
    mean_pos = positions.mean(axis=0)

    # Calculate RMSF
    diff = positions - mean_pos
    rmsf = np.sqrt((diff ** 2).sum(axis=2).mean(axis=0))

    print(f"  Mean RMSF: {np.mean(rmsf):.2f} A, Max: {np.max(rmsf):.2f} A")

    return resids, rmsf, resnames


def get_simulation_paths():
    """Define simulation paths for RMSF analysis - multiple replicates"""
    base = r'C:\Users\vasud\nod2-screening-data'

    simulations = {
        # WT controls
        'WT_Decoy': {
            'label': 'WT (Decoy/Unbound)',
            'pdb': os.path.join(base, 'PHASE_6', 'trajectories', 'decoy_rep1_solvated.pdb'),
            'dcd': os.path.join(base, 'PHASE_6', 'trajectories', 'decoy_rep1.dcd'),
            'color': 'gray',
            'linestyle': '--',
            'group': 'WT_apo'
        },
        # WT + Febuxostat (3 replicates)
        'WT_Febuxostat_rep1': {
            'label': 'WT + Feb rep1',
            'pdb': os.path.join(base, 'PHASE_6', 'trajectories', 'febuxostat_rep1_solvated.pdb'),
            'dcd': os.path.join(base, 'PHASE_6', 'trajectories', 'febuxostat_rep1.dcd'),
            'color': 'blue',
            'linestyle': ':',
            'group': 'WT_Feb'
        },
        'WT_Febuxostat_rep2': {
            'label': 'WT + Feb rep2',
            'pdb': os.path.join(base, 'PHASE_6', 'trajectories', 'febuxostat_rep2_solvated.pdb'),
            'dcd': os.path.join(base, 'PHASE_6', 'trajectories', 'febuxostat_rep2.dcd'),
            'color': 'blue',
            'linestyle': '-',
            'group': 'WT_Feb'
        },
        'WT_Febuxostat_rep3': {
            'label': 'WT + Feb rep3',
            'pdb': os.path.join(base, 'PHASE_6', 'trajectories', 'febuxostat_rep1_solvated.pdb'),
            'dcd': os.path.join(base, 'PHASE_6', 'trajectories', 'febuxostat_rep3.dcd'),
            'color': 'blue',
            'linestyle': '--',
            'group': 'WT_Feb'
        },
        # R702W + Febuxostat (3 replicates)
        'R702W_Febuxostat_rep1': {
            'label': 'R702W + Feb rep1',
            'pdb': os.path.join(base, 'PHASE_A1_mutant_MD', 'trajectories', 'R702W_febuxostat_rep1_solvated.pdb'),
            'dcd': os.path.join(base, 'PHASE_A1_mutant_MD', 'trajectories', 'R702W_febuxostat_rep1.dcd'),
            'color': 'red',
            'linestyle': ':',
            'group': 'R702W_Feb'
        },
        'R702W_Febuxostat_rep2': {
            'label': 'R702W + Feb rep2',
            'pdb': os.path.join(base, 'PHASE_A1_mutant_MD', 'trajectories', 'R702W_febuxostat_rep2_solvated.pdb'),
            'dcd': os.path.join(base, 'PHASE_A1_mutant_MD', 'trajectories', 'R702W_febuxostat_rep2.dcd'),
            'color': 'red',
            'linestyle': '-',
            'group': 'R702W_Feb'
        },
        'R702W_Febuxostat_rep3': {
            'label': 'R702W + Feb rep3',
            'pdb': os.path.join(base, 'PHASE_A1_mutant_MD', 'trajectories', 'R702W_febuxostat_rep3_solvated.pdb'),
            'dcd': os.path.join(base, 'PHASE_A1_mutant_MD', 'trajectories', 'R702W_febuxostat_rep3.dcd'),
            'color': 'red',
            'linestyle': '--',
            'group': 'R702W_Feb'
        },
        # G908R + Febuxostat (3 replicates)
        'G908R_Febuxostat_rep1': {
            'label': 'G908R + Feb rep1',
            'pdb': os.path.join(base, 'PHASE_A1_mutant_MD', 'trajectories', 'G908R_febuxostat_rep1_solvated.pdb'),
            'dcd': os.path.join(base, 'PHASE_A1_mutant_MD', 'trajectories', 'G908R_febuxostat_rep1.dcd'),
            'color': 'green',
            'linestyle': ':',
            'group': 'G908R_Feb'
        },
        'G908R_Febuxostat_rep2': {
            'label': 'G908R + Feb rep2',
            'pdb': os.path.join(base, 'PHASE_A1_mutant_MD', 'vast_downloads', '5080', 'G908R_febuxostat_rep2_solvated.pdb'),
            'dcd': os.path.join(base, 'PHASE_A1_mutant_MD', 'vast_downloads', '5080', 'G908R_febuxostat_rep2.dcd'),
            'color': 'green',
            'linestyle': '-',
            'group': 'G908R_Feb'
        },
        'G908R_Febuxostat_rep3': {
            'label': 'G908R + Feb rep3',
            'pdb': os.path.join(base, 'PHASE_A1_mutant_MD', 'vast_downloads', '5080', 'G908R_febuxostat_rep3_solvated.pdb'),
            'dcd': os.path.join(base, 'PHASE_A1_mutant_MD', 'vast_downloads', '5080', 'G908R_febuxostat_rep3.dcd'),
            'color': 'green',
            'linestyle': '--',
            'group': 'G908R_Feb'
        },
    }

    return simulations


def identify_lrr_domain(resids):
    """
    Identify LRR domain residues in NOD2
    NOD2 LRR domain is approximately residues 744-1040
    """
    # NOD2 domain boundaries (approximate):
    # CARD1: 1-93
    # CARD2: 104-191
    # NBD: 273-577
    # HD1: 578-628
    # HD2: 629-743
    # LRR: 744-1040

    lrr_start = 744
    lrr_end = 1040

    lrr_mask = (resids >= lrr_start) & (resids <= lrr_end)
    return lrr_mask, lrr_start, lrr_end


def identify_binding_site(resids):
    """
    Identify binding site residues (from docking/MD analysis)
    These are residues within 5A of the ligand
    """
    # Key binding site residues identified from previous analysis
    # (approximate - from typical NOD2-ligand interactions)
    binding_residues = [
        # LRR region binding residues (example set)
        750, 751, 752, 753, 775, 776, 777, 778, 800, 801, 802,
        825, 826, 827, 850, 851, 852, 875, 876, 877, 900, 901,
        925, 926, 950, 951, 975, 976, 1000, 1001
    ]

    binding_mask = np.isin(resids, binding_residues)
    return binding_mask, binding_residues


def main():
    base = r'C:\Users\vasud\nod2-screening-data'
    output_dir = os.path.join(base, 'PHASE_A1_mutant_MD', 'analysis', 'mmgbsa')

    print("=" * 80)
    print("RMSF ANALYSIS - PROTEIN FLEXIBILITY COMPARISON")
    print("=" * 80)
    print(f"Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("\nHypothesis: Drug binding reduces mutant NOD2 flexibility (chaperone effect)")

    simulations = get_simulation_paths()
    rmsf_data = {}

    # Calculate RMSF for each simulation
    for name, sim in simulations.items():
        print(f"\n{'='*60}")
        print(f"Analyzing: {sim['label']}")
        print(f"{'='*60}")

        if not os.path.exists(sim['dcd']):
            print(f"  ERROR: DCD not found")
            continue
        if not os.path.exists(sim['pdb']):
            print(f"  ERROR: PDB not found")
            continue

        try:
            resids, rmsf, resnames = calculate_rmsf(sim['pdb'], sim['dcd'])
            rmsf_data[name] = {
                'resids': resids,
                'rmsf': rmsf,
                'resnames': resnames,
                'label': sim['label'],
                'color': sim['color'],
                'linestyle': sim['linestyle']
            }
            print(f"  Mean RMSF: {np.mean(rmsf):.2f} A")
            print(f"  Max RMSF: {np.max(rmsf):.2f} A at residue {resids[np.argmax(rmsf)]}")
        except Exception as e:
            print(f"  ERROR: {e}")

    if len(rmsf_data) < 2:
        print("\nERROR: Not enough data for comparison")
        return

    # Get common residue range
    ref_resids = rmsf_data[list(rmsf_data.keys())[0]]['resids']

    # Identify LRR domain
    lrr_mask, lrr_start, lrr_end = identify_lrr_domain(ref_resids)

    # =========================================================================
    # SUMMARY STATISTICS (PER REPLICATE)
    # =========================================================================
    print("\n" + "=" * 80)
    print("RMSF SUMMARY STATISTICS (PER REPLICATE)")
    print("=" * 80)

    print(f"\n{'System':<25} | {'Mean RMSF':>10} | {'LRR RMSF':>10} | {'Max RMSF':>10}")
    print("-" * 70)

    summary_data = []
    group_data = {}

    for name, data in rmsf_data.items():
        mean_rmsf = np.mean(data['rmsf'])

        # LRR-specific RMSF
        lrr_mask_local, _, _ = identify_lrr_domain(data['resids'])
        if np.any(lrr_mask_local):
            lrr_rmsf = np.mean(data['rmsf'][lrr_mask_local])
        else:
            lrr_rmsf = np.nan

        max_rmsf = np.max(data['rmsf'])

        # Get group from simulations dict
        group = simulations.get(name, {}).get('group', 'unknown')

        summary_data.append({
            'system': name,
            'label': data['label'],
            'group': group,
            'mean_rmsf': mean_rmsf,
            'lrr_rmsf': lrr_rmsf,
            'max_rmsf': max_rmsf
        })

        # Collect for group averaging
        if group not in group_data:
            group_data[group] = {'mean': [], 'lrr': [], 'rmsf_arrays': []}
        group_data[group]['mean'].append(mean_rmsf)
        group_data[group]['lrr'].append(lrr_rmsf)
        group_data[group]['rmsf_arrays'].append(data['rmsf'])

        print(f"{data['label']:<25} | {mean_rmsf:>9.2f}A | {lrr_rmsf:>9.2f}A | {max_rmsf:>9.2f}A")

    # =========================================================================
    # GROUP AVERAGES
    # =========================================================================
    print("\n" + "=" * 80)
    print("RMSF GROUP AVERAGES (MEAN +/- SD ACROSS REPLICATES)")
    print("=" * 80)

    print(f"\n{'Group':<20} | {'Mean RMSF':>15} | {'LRR RMSF':>15} | {'N reps':>6}")
    print("-" * 70)

    group_summary = {}
    for group, gdata in sorted(group_data.items()):
        mean_avg = np.mean(gdata['mean'])
        mean_std = np.std(gdata['mean'])
        lrr_avg = np.mean([x for x in gdata['lrr'] if not np.isnan(x)])
        lrr_std = np.std([x for x in gdata['lrr'] if not np.isnan(x)])
        n_reps = len(gdata['mean'])

        group_summary[group] = {
            'mean_avg': mean_avg, 'mean_std': mean_std,
            'lrr_avg': lrr_avg, 'lrr_std': lrr_std,
            'n_reps': n_reps
        }

        print(f"{group:<20} | {mean_avg:>6.2f} +/- {mean_std:<5.2f} | {lrr_avg:>6.2f} +/- {lrr_std:<5.2f} | {n_reps:>6}")

    # =========================================================================
    # CHAPERONE HYPOTHESIS TEST (USING GROUP AVERAGES)
    # =========================================================================
    print("\n" + "=" * 80)
    print("CHAPERONE HYPOTHESIS TEST (GROUP AVERAGES)")
    print("=" * 80)

    # Compare WT+drug vs WT alone (Decoy)
    if 'WT_apo' in group_summary and 'WT_Feb' in group_summary:
        wt_alone = group_summary['WT_apo']['mean_avg']
        wt_drug = group_summary['WT_Feb']['mean_avg']
        wt_drug_std = group_summary['WT_Feb']['mean_std']
        delta_wt = wt_drug - wt_alone
        print(f"\nWT + Febuxostat vs WT alone (Decoy):")
        print(f"  WT alone (Decoy):  {wt_alone:.2f} A")
        print(f"  WT + Febuxostat:   {wt_drug:.2f} +/- {wt_drug_std:.2f} A (n={group_summary['WT_Feb']['n_reps']})")
        print(f"  Delta RMSF:        {delta_wt:+.2f} A")
        if delta_wt < -0.3:
            print(f"  --> Drug REDUCES flexibility (stabilizes)")
        elif delta_wt > 0.3:
            print(f"  --> Drug INCREASES flexibility")
        else:
            print(f"  --> No significant effect")

    # Compare mutant+drug vs WT+drug
    if 'WT_Feb' in group_summary:
        wt_drug = group_summary['WT_Feb']['mean_avg']
        wt_drug_std = group_summary['WT_Feb']['mean_std']

        for mut_group in ['R702W_Feb', 'G908R_Feb']:
            if mut_group in group_summary:
                mut_drug = group_summary[mut_group]['mean_avg']
                mut_drug_std = group_summary[mut_group]['mean_std']
                delta = mut_drug - wt_drug
                mut_name = mut_group.replace('_Feb', '')

                print(f"\n{mut_name} + Febuxostat vs WT + Febuxostat:")
                print(f"  WT + Febuxostat:     {wt_drug:.2f} +/- {wt_drug_std:.2f} A")
                print(f"  {mut_name} + Febuxostat:  {mut_drug:.2f} +/- {mut_drug_std:.2f} A")
                print(f"  Delta RMSF:          {delta:+.2f} A")

                if delta < -0.3:
                    print(f"  --> Mutant is MORE stable than WT with drug (supports chaperone)")
                elif delta > 0.3:
                    print(f"  --> Mutant is LESS stable than WT with drug")
                else:
                    print(f"  --> Similar stability")

    # =========================================================================
    # CREATE PLOT
    # =========================================================================
    print("\n" + "=" * 80)
    print("GENERATING RMSF PLOT")
    print("=" * 80)

    fig, axes = plt.subplots(2, 1, figsize=(14, 10))

    # Plot 1: Full protein RMSF
    ax1 = axes[0]
    for name, data in rmsf_data.items():
        ax1.plot(data['resids'], data['rmsf'],
                label=data['label'],
                color=data['color'],
                linestyle=data['linestyle'],
                linewidth=1.5, alpha=0.8)

    # Highlight LRR domain
    ax1.axvspan(lrr_start, lrr_end, alpha=0.15, color='yellow', label='LRR Domain')

    ax1.set_xlabel('Residue Number', fontsize=12)
    ax1.set_ylabel('RMSF (Angstrom)', fontsize=12)
    ax1.set_title('NOD2 C-alpha RMSF - Full Protein', fontsize=14, fontweight='bold')
    ax1.legend(loc='upper right', fontsize=10)
    ax1.set_xlim(min(ref_resids), max(ref_resids))
    ax1.set_ylim(0, None)
    ax1.grid(True, alpha=0.3)

    # Plot 2: LRR domain zoom
    ax2 = axes[1]
    for name, data in rmsf_data.items():
        lrr_mask_local, _, _ = identify_lrr_domain(data['resids'])
        if np.any(lrr_mask_local):
            ax2.plot(data['resids'][lrr_mask_local], data['rmsf'][lrr_mask_local],
                    label=data['label'],
                    color=data['color'],
                    linestyle=data['linestyle'],
                    linewidth=2, alpha=0.9)

    # Mark mutation sites
    ax2.axvline(x=702, color='red', linestyle=':', linewidth=2, alpha=0.7, label='R702W site')
    ax2.axvline(x=908, color='green', linestyle=':', linewidth=2, alpha=0.7, label='G908R site')

    ax2.set_xlabel('Residue Number', fontsize=12)
    ax2.set_ylabel('RMSF (Angstrom)', fontsize=12)
    ax2.set_title('NOD2 LRR Domain RMSF (residues 744-1040)', fontsize=14, fontweight='bold')
    ax2.legend(loc='upper right', fontsize=10)
    ax2.set_xlim(lrr_start, lrr_end)
    ax2.set_ylim(0, None)
    ax2.grid(True, alpha=0.3)

    plt.tight_layout()

    # Save plot
    plot_path = os.path.join(output_dir, 'RMSF_comparison.png')
    plt.savefig(plot_path, dpi=150, bbox_inches='tight')
    print(f"Saved: {plot_path}")

    # Also save as PDF for publication
    pdf_path = os.path.join(output_dir, 'RMSF_comparison.pdf')
    plt.savefig(pdf_path, bbox_inches='tight')
    print(f"Saved: {pdf_path}")

    plt.close()

    # =========================================================================
    # SAVE DATA TO CSV
    # =========================================================================
    # Create combined dataframe
    combined_data = []
    for name, data in rmsf_data.items():
        for i, (resid, rmsf_val) in enumerate(zip(data['resids'], data['rmsf'])):
            combined_data.append({
                'system': name,
                'resid': resid,
                'resname': data['resnames'][i],
                'rmsf': rmsf_val
            })

    df = pd.DataFrame(combined_data)
    csv_path = os.path.join(output_dir, 'rmsf_per_residue.csv')
    df.to_csv(csv_path, index=False)
    print(f"Saved: {csv_path}")

    # Save summary
    df_summary = pd.DataFrame(summary_data)
    summary_path = os.path.join(output_dir, 'rmsf_summary.csv')
    df_summary.to_csv(summary_path, index=False)
    print(f"Saved: {summary_path}")

    print("\n" + "=" * 80)
    print(f"Completed: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("=" * 80)

    return rmsf_data, summary_data


if __name__ == '__main__':
    rmsf_data, summary_data = main()
