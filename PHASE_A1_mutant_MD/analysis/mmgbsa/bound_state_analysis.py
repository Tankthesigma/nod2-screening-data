#!/usr/bin/env python
"""
Bound-State Analysis - Contact Statistics per Frame
Calculates min, median, max contacts and threshold percentages
"""

import os
import sys
import warnings
import numpy as np
import pandas as pd
from datetime import datetime
warnings.filterwarnings('ignore')

import MDAnalysis as mda
from MDAnalysis.analysis import distances


def count_contacts(universe, cutoff=5.0):
    """Count protein-ligand contacts for current frame"""
    protein = universe.select_atoms('protein')
    ligand = universe.select_atoms('resname UNK')

    if len(protein) == 0 or len(ligand) == 0:
        return 0

    box = universe.trajectory.ts.dimensions
    dist_matrix = distances.distance_array(
        ligand.positions, protein.positions, box=box
    )

    # Count contacts within cutoff
    contacts = np.sum(dist_matrix < cutoff)
    return contacts


def analyze_simulation_contacts(sim_info, base_dir, n_frames=50):
    """Analyze contact distribution for a single simulation"""
    pdb_path = os.path.join(base_dir, sim_info['pdb'])
    dcd_path = os.path.join(base_dir, sim_info['dcd'])

    u = mda.Universe(pdb_path, dcd_path)
    total_frames = len(u.trajectory)

    # Analyze frames 5-25% (same as energy analysis)
    start_frame = int(total_frames * 0.05)
    end_frame = int(total_frames * 0.25)
    stride = max(1, (end_frame - start_frame) // n_frames)

    frames_to_analyze = set(range(start_frame, end_frame, stride))

    # Collect contacts per frame
    contact_counts = []
    frame_count = 0

    for i, ts in enumerate(u.trajectory):
        if i > end_frame:
            break
        if i < start_frame or i not in frames_to_analyze:
            continue
        if frame_count >= n_frames:
            break

        contacts = count_contacts(u)
        contact_counts.append(contacts)
        frame_count += 1

    contact_counts = np.array(contact_counts)

    # Calculate statistics
    result = {
        'num': sim_info['num'],
        'system': sim_info['system'],
        'ligand': sim_info['ligand'],
        'rep': sim_info['rep'],
        'n_frames': len(contact_counts),
        'min_contacts': int(np.min(contact_counts)),
        'median_contacts': int(np.median(contact_counts)),
        'max_contacts': int(np.max(contact_counts)),
        'pct_gte_50': 100 * np.sum(contact_counts >= 50) / len(contact_counts),
        'pct_gte_100': 100 * np.sum(contact_counts >= 100) / len(contact_counts),
        'pct_gte_75': 100 * np.sum(contact_counts >= 75) / len(contact_counts),
        'pct_gte_125': 100 * np.sum(contact_counts >= 125) / len(contact_counts),
        'contact_counts': contact_counts  # Keep raw data
    }

    return result


def get_all_simulations():
    """Return all 25 simulation configurations - CORRECT PATHS"""
    sims = []

    # Phase 6 - Wild Type (13 sims)
    wt_sims = [
        (1, 'Febuxostat', 1, 'febuxostat_rep1', 'febuxostat_rep1_solvated.pdb', 'Lead'),
        (2, 'Febuxostat', 2, 'febuxostat_rep2', 'febuxostat_rep2_solvated.pdb', ''),
        (3, 'Febuxostat', 3, 'febuxostat_rep3', 'febuxostat_rep1_solvated.pdb', 'Using rep1 PDB'),
        (4, 'Natural', 1, 'natural_top_rep1', 'natural_top_rep1_solvated.pdb', 'Backup'),
        (5, 'Natural', 2, 'natural_top_rep2', 'natural_top_rep2_solvated_new.pdb', ''),
        (6, 'Natural', 3, 'natural_top_rep3', 'natural_top_rep3_solvated.pdb', ''),
        (7, 'Budesonide', 1, 'budesonide_rep1', 'budesonide_rep1_solvated.pdb', 'Pos Control'),
        (8, 'Budesonide', 2, 'budesonide_rep2', 'budesonide_rep2_solvated.pdb', ''),
        (9, 'Budesonide', 3, 'budesonide_rep3', 'budesonide_rep3_solvated.pdb', ''),
        (10, 'Ursodiol', 1, 'ursodiol_rep1', 'ursodiol_rep3_solvated.pdb', 'Using rep3 PDB'),
        (11, 'Ursodiol', 2, 'ursodiol_rep2', 'ursodiol_rep3_solvated.pdb', 'Using rep3 PDB'),
        (12, 'Ursodiol', 3, 'ursodiol_rep3', 'ursodiol_rep3_solvated.pdb', ''),
        (13, 'Decoy', 1, 'decoy_rep1', 'decoy_rep1_solvated.pdb', 'NEG CONTROL'),
    ]

    for num, ligand, rep, dcd_base, pdb_file, note in wt_sims:
        sims.append({
            'num': num,
            'phase': 'Phase6',
            'system': 'WT',
            'ligand': ligand,
            'rep': rep,
            'dcd': os.path.join('PHASE_6', 'trajectories', f'{dcd_base}.dcd'),
            'pdb': os.path.join('PHASE_6', 'trajectories', pdb_file),
            'note': note
        })

    # Phase A1 - Mutants (12 sims) - CORRECT PATHS
    mutant_sims = [
        # R702W Febuxostat
        (14, 'R702W', 'Febuxostat', 1,
         os.path.join('PHASE_A1_mutant_MD', 'trajectories', 'R702W_febuxostat_rep1.dcd'),
         os.path.join('PHASE_A1_mutant_MD', 'trajectories', 'R702W_febuxostat_rep1_solvated.pdb'), ''),
        (15, 'R702W', 'Febuxostat', 2,
         os.path.join('PHASE_A1_mutant_MD', 'trajectories', 'R702W_febuxostat_rep2.dcd'),
         os.path.join('PHASE_A1_mutant_MD', 'trajectories', 'R702W_febuxostat_rep2_solvated.pdb'), ''),
        (16, 'R702W', 'Febuxostat', 3,
         os.path.join('PHASE_A1_mutant_MD', 'trajectories', 'R702W_febuxostat_rep3.dcd'),
         os.path.join('PHASE_A1_mutant_MD', 'trajectories', 'R702W_febuxostat_rep3_solvated.pdb'), ''),
        # R702W Natural
        (17, 'R702W', 'Natural', 1,
         os.path.join('PHASE_A1_mutant_MD', 'vast_downloads', '5080', 'R702W_natural_rep1.dcd'),
         os.path.join('PHASE_A1_mutant_MD', 'vast_downloads', '5080', 'R702W_natural_rep1_solvated.pdb'), ''),
        (18, 'R702W', 'Natural', 2,
         os.path.join('PHASE_A1_mutant_MD', 'vast_downloads', '5080', 'R702W_natural_rep2.dcd'),
         os.path.join('PHASE_A1_mutant_MD', 'vast_downloads', '5080', 'R702W_natural_rep2_solvated.pdb'), ''),
        (19, 'R702W', 'Natural', 3,
         os.path.join('PHASE_A1_mutant_MD', 'vast_downloads', 'new_5090', 'R702W_natural_rep3.dcd'),
         os.path.join('PHASE_A1_mutant_MD', 'vast_downloads', 'new_5090', 'R702W_natural_rep3_solvated.pdb'), ''),
        # G908R Febuxostat
        (20, 'G908R', 'Febuxostat', 1,
         os.path.join('PHASE_A1_mutant_MD', 'trajectories', 'G908R_febuxostat_rep1.dcd'),
         os.path.join('PHASE_A1_mutant_MD', 'trajectories', 'G908R_febuxostat_rep1_solvated.pdb'), ''),
        (21, 'G908R', 'Febuxostat', 2,
         os.path.join('PHASE_A1_mutant_MD', 'vast_downloads', '5080', 'G908R_febuxostat_rep2.dcd'),
         os.path.join('PHASE_A1_mutant_MD', 'vast_downloads', '5080', 'G908R_febuxostat_rep2_solvated.pdb'), ''),
        (22, 'G908R', 'Febuxostat', 3,
         os.path.join('PHASE_A1_mutant_MD', 'vast_downloads', '5080', 'G908R_febuxostat_rep3.dcd'),
         os.path.join('PHASE_A1_mutant_MD', 'vast_downloads', '5080', 'G908R_febuxostat_rep3_solvated.pdb'), ''),
        # G908R Natural
        (23, 'G908R', 'Natural', 1,
         os.path.join('PHASE_A1_mutant_MD', 'vast_downloads', '5090', 'G908R_natural_rep1.dcd'),
         os.path.join('PHASE_A1_mutant_MD', 'vast_downloads', '5090', 'G908R_natural_rep1_solvated.pdb'), ''),
        (24, 'G908R', 'Natural', 2,
         os.path.join('PHASE_A1_mutant_MD', 'vast_downloads', 'new_5090', 'G908R_natural_rep2.dcd'),
         os.path.join('PHASE_A1_mutant_MD', 'vast_downloads', 'new_5090', 'G908R_natural_rep2_solvated.pdb'), ''),
        (25, 'G908R', 'Natural', 3,
         os.path.join('PHASE_A1_mutant_MD', 'vast_downloads', 'new_5090', 'G908R_natural_rep3.dcd'),
         os.path.join('PHASE_A1_mutant_MD', 'vast_downloads', 'new_5090', 'G908R_natural_rep3_solvated.pdb'), ''),
    ]

    for num, system, ligand, rep, dcd, pdb, note in mutant_sims:
        sims.append({
            'num': num,
            'phase': 'PhaseA1',
            'system': system,
            'ligand': ligand,
            'rep': rep,
            'dcd': dcd,
            'pdb': pdb,
            'note': note
        })

    return sims


def main():
    base = r'C:\Users\vasud\nod2-screening-data'
    output_dir = os.path.join(base, 'PHASE_A1_mutant_MD', 'analysis', 'mmgbsa')

    print("=" * 80)
    print("BOUND-STATE ANALYSIS - CONTACT STATISTICS")
    print("=" * 80)
    print(f"Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")

    all_sims = get_all_simulations()
    results = []
    all_contacts = {}  # Store raw contact data

    for sim in all_sims:
        name = f"{sim['system']}_{sim['ligand']}_rep{sim['rep']}"
        print(f"\n[{sim['num']:2d}/25] {name}")

        try:
            result = analyze_simulation_contacts(sim, base, n_frames=50)
            results.append(result)
            all_contacts[name] = result['contact_counts']

            print(f"  Min={result['min_contacts']:3d}  Med={result['median_contacts']:3d}  Max={result['max_contacts']:3d}  "
                  f"%>=50: {result['pct_gte_50']:5.1f}%  %>=100: {result['pct_gte_100']:5.1f}%")
        except Exception as e:
            print(f"  ERROR: {e}")

    # Create output table
    print("\n" + "=" * 80)
    print("CONTACT STATISTICS TABLE")
    print("=" * 80)
    print(f"\n{'#':>2} | {'System':<6} | {'Ligand':<10} | {'Rep':<3} | {'Min':>4} | {'Med':>4} | {'Max':>4} | {'%>=50':>6} | {'%>=100':>6}")
    print("-" * 80)

    for r in results:
        print(f"{r['num']:2d} | {r['system']:<6} | {r['ligand']:<10} | {r['rep']:<3} | "
              f"{r['min_contacts']:4d} | {r['median_contacts']:4d} | {r['max_contacts']:4d} | "
              f"{r['pct_gte_50']:5.1f}% | {r['pct_gte_100']:5.1f}%")

    # Save to CSV
    df = pd.DataFrame([{k: v for k, v in r.items() if k != 'contact_counts'} for r in results])
    csv_path = os.path.join(output_dir, 'bound_state_contact_stats.csv')
    df.to_csv(csv_path, index=False)
    print(f"\nSaved: {csv_path}")

    # Threshold analysis
    print("\n" + "=" * 80)
    print("THRESHOLD ANALYSIS")
    print("=" * 80)

    # Identify stable vs outlier replicates
    stable_reps = []
    outlier_reps = []

    for r in results:
        name = f"{r['system']}_{r['ligand']}_rep{r['rep']}"
        if r['ligand'] == 'Decoy':
            continue  # Skip decoy

        # Outliers: median < 100 contacts
        if r['median_contacts'] < 100:
            outlier_reps.append(name)
        else:
            stable_reps.append(name)

    print(f"\nStable replicates ({len(stable_reps)}): median >= 100 contacts")
    print(f"Outlier replicates ({len(outlier_reps)}): median < 100 contacts")
    print(f"  Outliers: {', '.join(outlier_reps)}")

    # Calculate threshold statistics
    thresholds = [25, 50, 75, 100, 125, 150]

    print(f"\n{'Threshold':<12} | {'Stable % frames':<18} | {'Outlier % frames':<18}")
    print("-" * 55)

    for thresh in thresholds:
        stable_pct = []
        outlier_pct = []

        for r in results:
            name = f"{r['system']}_{r['ligand']}_rep{r['rep']}"
            if r['ligand'] == 'Decoy':
                continue

            pct = 100 * np.sum(all_contacts[name] >= thresh) / len(all_contacts[name])

            if name in stable_reps:
                stable_pct.append(pct)
            else:
                outlier_pct.append(pct)

        stable_mean = np.mean(stable_pct) if stable_pct else 0
        outlier_mean = np.mean(outlier_pct) if outlier_pct else 0

        print(f">= {thresh:<9} | {stable_mean:>6.1f}% (mean)     | {outlier_mean:>6.1f}% (mean)")

    # Recommendation
    print("\n" + "=" * 80)
    print("RECOMMENDATION")
    print("=" * 80)

    # Find threshold where stable > 80% and outlier < 50%
    print("\nLooking for threshold where:")
    print("  - Stable replicates: >= 80% of frames pass")
    print("  - Outlier replicates: < 50% of frames pass")

    for thresh in [50, 75, 100, 125]:
        stable_pct = []
        outlier_pct = []

        for r in results:
            name = f"{r['system']}_{r['ligand']}_rep{r['rep']}"
            if r['ligand'] == 'Decoy':
                continue

            pct = 100 * np.sum(all_contacts[name] >= thresh) / len(all_contacts[name])

            if name in stable_reps:
                stable_pct.append(pct)
            else:
                outlier_pct.append(pct)

        stable_mean = np.mean(stable_pct)
        outlier_mean = np.mean(outlier_pct)

        if stable_mean >= 80 and outlier_mean < 50:
            print(f"\n  ** RECOMMENDED THRESHOLD: {thresh} contacts **")
            print(f"     Stable replicates: {stable_mean:.1f}% of frames pass")
            print(f"     Outlier replicates: {outlier_mean:.1f}% of frames pass")
            break

    print("\n" + "=" * 80)
    print(f"Completed: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("=" * 80)

    return results


if __name__ == '__main__':
    results = main()
