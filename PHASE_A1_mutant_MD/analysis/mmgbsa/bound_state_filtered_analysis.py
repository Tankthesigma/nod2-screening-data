#!/usr/bin/env python
"""
Bound-State Filtered Analysis
Apply contacts >= 75 filter and compare All-Frames vs Bound-Only
Includes Bootstrap 95% CI for DeltaDeltaE
"""

import os
import sys
import warnings
import numpy as np
import pandas as pd
from datetime import datetime
from collections import Counter
warnings.filterwarnings('ignore')

import MDAnalysis as mda
from MDAnalysis.analysis import distances

# Contact threshold
CONTACT_THRESHOLD = 75


class BoundStateAnalyzer:
    """Analyze binding energies with bound-state filtering"""

    def __init__(self, base_dir):
        self.base_dir = base_dir
        self.results = []

        # Contact energy parameters
        self.contact_energies = {
            'hydrophobic': -0.3,
            'hbond': -1.5,
            'ionic': -2.0,
            'clash': 0.5,
            'neutral': -0.1,
        }

    def classify_contact(self, lig_elem, prot_elem, distance):
        if distance > 5.0:
            return None
        if distance < 1.5:
            return 'clash'
        if lig_elem == 'C' and prot_elem == 'C':
            if distance < 4.5:
                return 'hydrophobic'
        if lig_elem in ['N', 'O'] and prot_elem in ['N', 'O']:
            if distance < 3.5:
                return 'hbond'
        if distance < 4.0:
            return 'neutral'
        return None

    def calculate_frame_energy(self, universe, cutoff=6.0):
        """Calculate energy and contacts for current frame"""
        protein = universe.select_atoms('protein')
        ligand = universe.select_atoms('resname UNK')

        if len(protein) == 0 or len(ligand) == 0:
            return {'energy': np.nan, 'contacts': 0}

        box = universe.trajectory.ts.dimensions
        dist_matrix = distances.distance_array(
            ligand.positions, protein.positions, box=box
        )

        lig_elems = [a.element if hasattr(a, 'element') and a.element else a.name[0]
                     for a in ligand]
        prot_elems = [a.element if hasattr(a, 'element') and a.element else a.name[0]
                      for a in protein]

        contact_counts = Counter()
        total_energy = 0.0

        for i, lig_elem in enumerate(lig_elems):
            for j, prot_elem in enumerate(prot_elems):
                r = dist_matrix[i, j]
                if r > cutoff:
                    continue
                contact_type = self.classify_contact(lig_elem, prot_elem, r)
                if contact_type:
                    contact_counts[contact_type] += 1
                    total_energy += self.contact_energies.get(contact_type, 0)

        return {
            'energy': total_energy,
            'contacts': sum(contact_counts.values()),
        }

    def analyze_simulation(self, sim_info, n_frames=50):
        """Analyze simulation with per-frame data"""
        name = f"{sim_info['system']}_{sim_info['ligand']}_rep{sim_info['rep']}"

        pdb_path = os.path.join(self.base_dir, sim_info['pdb'])
        dcd_path = os.path.join(self.base_dir, sim_info['dcd'])

        if not os.path.exists(dcd_path) or not os.path.exists(pdb_path):
            print(f"[{sim_info['num']:2d}] {name}: FILE NOT FOUND")
            return None

        try:
            u = mda.Universe(pdb_path, dcd_path)
            total_frames = len(u.trajectory)
        except Exception as e:
            print(f"[{sim_info['num']:2d}] {name}: ERROR - {e}")
            return None

        # Analyze frames 5-25%
        start_frame = int(total_frames * 0.05)
        end_frame = int(total_frames * 0.25)
        stride = max(1, (end_frame - start_frame) // n_frames)
        frames_to_analyze = set(range(start_frame, end_frame, stride))

        # Collect per-frame data
        frame_data = []
        frame_count = 0

        for i, ts in enumerate(u.trajectory):
            if i > end_frame:
                break
            if i < start_frame or i not in frames_to_analyze:
                continue
            if frame_count >= n_frames:
                break

            result = self.calculate_frame_energy(u)
            frame_data.append(result)
            frame_count += 1

        if not frame_data:
            return None

        # Convert to arrays
        energies = np.array([f['energy'] for f in frame_data])
        contacts = np.array([f['contacts'] for f in frame_data])

        # All-frames statistics
        all_mean = np.mean(energies)
        all_std = np.std(energies)
        all_n = len(energies)

        # Bound-only statistics (contacts >= threshold)
        bound_mask = contacts >= CONTACT_THRESHOLD
        bound_energies = energies[bound_mask]
        bound_contacts = contacts[bound_mask]

        if len(bound_energies) > 0:
            bound_mean = np.mean(bound_energies)
            bound_std = np.std(bound_energies)
            bound_n = len(bound_energies)
            bound_contacts_mean = np.mean(bound_contacts)
        else:
            bound_mean = np.nan
            bound_std = np.nan
            bound_n = 0
            bound_contacts_mean = 0

        pct_bound = 100 * bound_n / all_n

        result = {
            'num': sim_info['num'],
            'system': sim_info['system'],
            'ligand': sim_info['ligand'],
            'rep': sim_info['rep'],
            'all_mean': all_mean,
            'all_std': all_std,
            'all_n': all_n,
            'bound_mean': bound_mean,
            'bound_std': bound_std,
            'bound_n': bound_n,
            'pct_bound': pct_bound,
            'bound_contacts_mean': bound_contacts_mean,
            'bound_energies': bound_energies,  # Keep for bootstrap
        }

        print(f"[{sim_info['num']:2d}] {name}: All={all_mean:.1f}, Bound={bound_mean:.1f} ({pct_bound:.0f}% pass)")

        self.results.append(result)
        return result


def get_all_simulations():
    """Return all 25 simulation configurations"""
    sims = []

    # Phase 6 - Wild Type
    wt_sims = [
        (1, 'Febuxostat', 1, 'febuxostat_rep1', 'febuxostat_rep1_solvated.pdb'),
        (2, 'Febuxostat', 2, 'febuxostat_rep2', 'febuxostat_rep2_solvated.pdb'),
        (3, 'Febuxostat', 3, 'febuxostat_rep3', 'febuxostat_rep1_solvated.pdb'),
        (4, 'Natural', 1, 'natural_top_rep1', 'natural_top_rep1_solvated.pdb'),
        (5, 'Natural', 2, 'natural_top_rep2', 'natural_top_rep2_solvated_new.pdb'),
        (6, 'Natural', 3, 'natural_top_rep3', 'natural_top_rep3_solvated.pdb'),
        (7, 'Budesonide', 1, 'budesonide_rep1', 'budesonide_rep1_solvated.pdb'),
        (8, 'Budesonide', 2, 'budesonide_rep2', 'budesonide_rep2_solvated.pdb'),
        (9, 'Budesonide', 3, 'budesonide_rep3', 'budesonide_rep3_solvated.pdb'),
        (10, 'Ursodiol', 1, 'ursodiol_rep1', 'ursodiol_rep3_solvated.pdb'),
        (11, 'Ursodiol', 2, 'ursodiol_rep2', 'ursodiol_rep3_solvated.pdb'),
        (12, 'Ursodiol', 3, 'ursodiol_rep3', 'ursodiol_rep3_solvated.pdb'),
        (13, 'Decoy', 1, 'decoy_rep1', 'decoy_rep1_solvated.pdb'),
    ]

    for num, ligand, rep, dcd_base, pdb_file in wt_sims:
        sims.append({
            'num': num, 'system': 'WT', 'ligand': ligand, 'rep': rep,
            'dcd': os.path.join('PHASE_6', 'trajectories', f'{dcd_base}.dcd'),
            'pdb': os.path.join('PHASE_6', 'trajectories', pdb_file),
        })

    # Phase A1 - Mutants
    mutant_sims = [
        (14, 'R702W', 'Febuxostat', 1, 'PHASE_A1_mutant_MD/trajectories/R702W_febuxostat_rep1.dcd',
         'PHASE_A1_mutant_MD/trajectories/R702W_febuxostat_rep1_solvated.pdb'),
        (15, 'R702W', 'Febuxostat', 2, 'PHASE_A1_mutant_MD/trajectories/R702W_febuxostat_rep2.dcd',
         'PHASE_A1_mutant_MD/trajectories/R702W_febuxostat_rep2_solvated.pdb'),
        (16, 'R702W', 'Febuxostat', 3, 'PHASE_A1_mutant_MD/trajectories/R702W_febuxostat_rep3.dcd',
         'PHASE_A1_mutant_MD/trajectories/R702W_febuxostat_rep3_solvated.pdb'),
        (17, 'R702W', 'Natural', 1, 'PHASE_A1_mutant_MD/vast_downloads/5080/R702W_natural_rep1.dcd',
         'PHASE_A1_mutant_MD/vast_downloads/5080/R702W_natural_rep1_solvated.pdb'),
        (18, 'R702W', 'Natural', 2, 'PHASE_A1_mutant_MD/vast_downloads/5080/R702W_natural_rep2.dcd',
         'PHASE_A1_mutant_MD/vast_downloads/5080/R702W_natural_rep2_solvated.pdb'),
        (19, 'R702W', 'Natural', 3, 'PHASE_A1_mutant_MD/vast_downloads/new_5090/R702W_natural_rep3.dcd',
         'PHASE_A1_mutant_MD/vast_downloads/new_5090/R702W_natural_rep3_solvated.pdb'),
        (20, 'G908R', 'Febuxostat', 1, 'PHASE_A1_mutant_MD/trajectories/G908R_febuxostat_rep1.dcd',
         'PHASE_A1_mutant_MD/trajectories/G908R_febuxostat_rep1_solvated.pdb'),
        (21, 'G908R', 'Febuxostat', 2, 'PHASE_A1_mutant_MD/vast_downloads/5080/G908R_febuxostat_rep2.dcd',
         'PHASE_A1_mutant_MD/vast_downloads/5080/G908R_febuxostat_rep2_solvated.pdb'),
        (22, 'G908R', 'Febuxostat', 3, 'PHASE_A1_mutant_MD/vast_downloads/5080/G908R_febuxostat_rep3.dcd',
         'PHASE_A1_mutant_MD/vast_downloads/5080/G908R_febuxostat_rep3_solvated.pdb'),
        (23, 'G908R', 'Natural', 1, 'PHASE_A1_mutant_MD/vast_downloads/5090/G908R_natural_rep1.dcd',
         'PHASE_A1_mutant_MD/vast_downloads/5090/G908R_natural_rep1_solvated.pdb'),
        (24, 'G908R', 'Natural', 2, 'PHASE_A1_mutant_MD/vast_downloads/new_5090/G908R_natural_rep2.dcd',
         'PHASE_A1_mutant_MD/vast_downloads/new_5090/G908R_natural_rep2_solvated.pdb'),
        (25, 'G908R', 'Natural', 3, 'PHASE_A1_mutant_MD/vast_downloads/new_5090/G908R_natural_rep3.dcd',
         'PHASE_A1_mutant_MD/vast_downloads/new_5090/G908R_natural_rep3_solvated.pdb'),
    ]

    for num, system, ligand, rep, dcd, pdb in mutant_sims:
        sims.append({
            'num': num, 'system': system, 'ligand': ligand, 'rep': rep,
            'dcd': dcd, 'pdb': pdb,
        })

    return sims


def bootstrap_ci(data1, data2, n_bootstrap=1000, ci=95):
    """Calculate bootstrap CI for difference of means"""
    np.random.seed(42)
    diffs = []

    for _ in range(n_bootstrap):
        sample1 = np.random.choice(data1, size=len(data1), replace=True)
        sample2 = np.random.choice(data2, size=len(data2), replace=True)
        diffs.append(np.mean(sample1) - np.mean(sample2))

    lower = np.percentile(diffs, (100 - ci) / 2)
    upper = np.percentile(diffs, 100 - (100 - ci) / 2)
    return np.mean(diffs), lower, upper


def main():
    base = r'C:\Users\vasud\nod2-screening-data'
    output_dir = os.path.join(base, 'PHASE_A1_mutant_MD', 'analysis', 'mmgbsa')

    print("=" * 80)
    print(f"BOUND-STATE FILTERED ANALYSIS (threshold >= {CONTACT_THRESHOLD} contacts)")
    print("=" * 80)
    print(f"Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")

    analyzer = BoundStateAnalyzer(base)
    all_sims = get_all_simulations()

    # Run all analyses
    for sim in all_sims:
        analyzer.analyze_simulation(sim, n_frames=50)

    results = analyzer.results

    # =========================================================================
    # TABLE 1: Per-Replicate Comparison
    # =========================================================================
    print("\n" + "=" * 100)
    print("TABLE 1: PER-REPLICATE COMPARISON")
    print("=" * 100)
    print(f"{'#':>2} | {'System':<6} | {'Ligand':<10} | {'Rep':<3} | {'All-Frames':<12} | {'Bound-Only':<12} | {'%Bound':>7} | {'Contacts':>8}")
    print("-" * 100)

    for r in results:
        all_str = f"{r['all_mean']:.1f} +/- {r['all_std']:.1f}"
        if r['bound_n'] > 0:
            bound_str = f"{r['bound_mean']:.1f} +/- {r['bound_std']:.1f}"
        else:
            bound_str = "N/A"
        print(f"{r['num']:2d} | {r['system']:<6} | {r['ligand']:<10} | {r['rep']:<3} | "
              f"{all_str:<12} | {bound_str:<12} | {r['pct_bound']:6.1f}% | {r['bound_contacts_mean']:7.0f}")

    # =========================================================================
    # TABLE 2: Condition Means
    # =========================================================================
    print("\n" + "=" * 100)
    print("TABLE 2: CONDITION MEANS (for DeltaDeltaE calculation)")
    print("=" * 100)

    # Group by system + ligand
    conditions = {}
    for r in results:
        if r['ligand'] == 'Decoy':
            continue
        key = (r['system'], r['ligand'])
        if key not in conditions:
            conditions[key] = {'all': [], 'bound': [], 'bound_energies': []}
        conditions[key]['all'].append(r['all_mean'])
        if r['bound_n'] > 0:
            conditions[key]['bound'].append(r['bound_mean'])
            conditions[key]['bound_energies'].extend(r['bound_energies'])

    print(f"{'System':<8} | {'Ligand':<10} | {'All-Frames Mean':<16} | {'Bound-Only Mean':<16} | {'SD':<8} | {'SEM':<8} | {'N reps':>6}")
    print("-" * 100)

    condition_stats = {}
    for (system, ligand), data in sorted(conditions.items()):
        all_mean = np.mean(data['all'])
        all_sd = np.std(data['all'])

        if len(data['bound']) > 0:
            bound_mean = np.mean(data['bound'])
            bound_sd = np.std(data['bound'])
            bound_sem = bound_sd / np.sqrt(len(data['bound']))
        else:
            bound_mean = np.nan
            bound_sd = np.nan
            bound_sem = np.nan

        condition_stats[(system, ligand)] = {
            'all_mean': all_mean,
            'bound_mean': bound_mean,
            'bound_sd': bound_sd,
            'bound_sem': bound_sem,
            'bound_energies': np.array(data['bound_energies']),
            'n_reps': len(data['bound'])
        }

        print(f"{system:<8} | {ligand:<10} | {all_mean:>7.1f} +/- {all_sd:<5.1f}  | "
              f"{bound_mean:>7.1f} +/- {bound_sd:<5.1f}  | {bound_sd:>6.1f}  | {bound_sem:>6.2f}  | {len(data['bound']):>6}")

    # =========================================================================
    # TABLE 3: DeltaDeltaE with Both Methods
    # =========================================================================
    print("\n" + "=" * 100)
    print("TABLE 3: DeltaDeltaE (MUTATION EFFECTS)")
    print("=" * 100)
    print(f"{'Mutation':<10} | {'Ligand':<10} | {'DDE All-Frames':>14} | {'DDE Bound-Only':>14} | {'Interpretation':<20}")
    print("-" * 100)

    comparisons = [
        ('R702W', 'Febuxostat'),
        ('R702W', 'Natural'),
        ('G908R', 'Febuxostat'),
        ('G908R', 'Natural'),
    ]

    dde_results = []
    for mutation, ligand in comparisons:
        wt_stats = condition_stats.get(('WT', ligand))
        mut_stats = condition_stats.get((mutation, ligand))

        if wt_stats and mut_stats:
            dde_all = mut_stats['all_mean'] - wt_stats['all_mean']
            dde_bound = mut_stats['bound_mean'] - wt_stats['bound_mean']

            # Interpretation
            if abs(dde_bound) < 2:
                interp = "Similar"
            elif dde_bound < -2:
                interp = "Stronger binding"
            else:
                interp = "Weaker binding"

            dde_results.append({
                'mutation': mutation,
                'ligand': ligand,
                'dde_all': dde_all,
                'dde_bound': dde_bound,
                'wt_energies': wt_stats['bound_energies'],
                'mut_energies': mut_stats['bound_energies'],
            })

            print(f"{mutation:<10} | {ligand:<10} | {dde_all:>+13.1f}  | {dde_bound:>+13.1f}  | {interp:<20}")

    # =========================================================================
    # BOOTSTRAP 95% CI for Bound-Only DeltaDeltaE
    # =========================================================================
    print("\n" + "=" * 100)
    print("BOOTSTRAP 95% CI FOR BOUND-ONLY DeltaDeltaE (1000 resamples)")
    print("=" * 100)
    print(f"{'Comparison':<25} | {'DDE Mean':>10} | {'95% CI Lower':>12} | {'95% CI Upper':>12} | {'Significant?':<12}")
    print("-" * 100)

    for dde in dde_results:
        if len(dde['wt_energies']) > 0 and len(dde['mut_energies']) > 0:
            mean_diff, ci_lower, ci_upper = bootstrap_ci(
                dde['mut_energies'], dde['wt_energies'], n_bootstrap=1000
            )

            # Check if CI excludes zero
            if ci_lower > 0 or ci_upper < 0:
                significant = "YES"
            else:
                significant = "NO"

            comparison = f"{dde['mutation']} vs WT ({dde['ligand']})"
            print(f"{comparison:<25} | {mean_diff:>+9.1f}  | {ci_lower:>+11.1f}  | {ci_upper:>+11.1f}  | {significant:<12}")

    # Save results
    df_results = pd.DataFrame([{k: v for k, v in r.items() if k != 'bound_energies'} for r in results])
    df_results.to_csv(os.path.join(output_dir, 'bound_state_filtered_results.csv'), index=False)

    print("\n" + "=" * 100)
    print(f"Completed: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"Saved: bound_state_filtered_results.csv")
    print("=" * 100)

    return results, condition_stats, dde_results


if __name__ == '__main__':
    results, condition_stats, dde_results = main()
