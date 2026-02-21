#!/usr/bin/env python
"""
Convergence Analysis using Block Averaging
Calculate DeltaE_GBSA in 4 equal time blocks for all 25 simulations
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

# Convergence threshold
CONVERGENCE_THRESHOLD = 6.0  # kcal/mol


class BlockAverageAnalyzer:
    """Analyze convergence using block averaging"""

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
        """Calculate energy for current frame"""
        protein = universe.select_atoms('protein')
        ligand = universe.select_atoms('resname UNK')

        if len(protein) == 0 or len(ligand) == 0:
            return np.nan

        box = universe.trajectory.ts.dimensions
        dist_matrix = distances.distance_array(
            ligand.positions, protein.positions, box=box
        )

        lig_elems = [a.element if hasattr(a, 'element') and a.element else a.name[0]
                     for a in ligand]
        prot_elems = [a.element if hasattr(a, 'element') and a.element else a.name[0]
                      for a in protein]

        total_energy = 0.0
        for i, lig_elem in enumerate(lig_elems):
            for j, prot_elem in enumerate(prot_elems):
                r = dist_matrix[i, j]
                if r > cutoff:
                    continue
                contact_type = self.classify_contact(lig_elem, prot_elem, r)
                if contact_type:
                    total_energy += self.contact_energies.get(contact_type, 0)

        return total_energy

    def analyze_simulation_blocks(self, sim_info, n_frames_per_block=25):
        """Analyze simulation in 4 blocks"""
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

        # Define 4 blocks spanning entire trajectory
        block_frames = total_frames // 4
        blocks = [
            (0, block_frames),
            (block_frames, 2 * block_frames),
            (2 * block_frames, 3 * block_frames),
            (3 * block_frames, total_frames)
        ]

        block_energies = {1: [], 2: [], 3: [], 4: []}

        # Sample frames from each block
        for block_idx, (start, end) in enumerate(blocks, 1):
            stride = max(1, (end - start) // n_frames_per_block)
            frames_to_analyze = list(range(start, end, stride))[:n_frames_per_block]

            frame_count = 0
            for i, ts in enumerate(u.trajectory):
                if i in frames_to_analyze:
                    energy = self.calculate_frame_energy(u)
                    if not np.isnan(energy):
                        block_energies[block_idx].append(energy)
                    frame_count += 1
                    if frame_count >= n_frames_per_block:
                        break

        # Calculate block means
        block_means = {}
        for b in range(1, 5):
            if block_energies[b]:
                block_means[b] = np.mean(block_energies[b])
            else:
                block_means[b] = np.nan

        # Calculate range and convergence
        valid_means = [v for v in block_means.values() if not np.isnan(v)]
        if len(valid_means) >= 2:
            energy_range = max(valid_means) - min(valid_means)
            converged = energy_range < CONVERGENCE_THRESHOLD
        else:
            energy_range = np.nan
            converged = False

        # Detect patterns
        pattern = self.detect_pattern(block_means)

        result = {
            'num': sim_info['num'],
            'system': sim_info['system'],
            'ligand': sim_info['ligand'],
            'rep': sim_info['rep'],
            'block1': block_means[1],
            'block2': block_means[2],
            'block3': block_means[3],
            'block4': block_means[4],
            'range': energy_range,
            'converged': 'YES' if converged else 'NO',
            'pattern': pattern
        }

        status = "OK" if converged else "DRIFT"
        print(f"[{sim_info['num']:2d}] {name}: B1={block_means[1]:.1f} B2={block_means[2]:.1f} "
              f"B3={block_means[3]:.1f} B4={block_means[4]:.1f} | Range={energy_range:.1f} | {status}")

        self.results.append(result)
        return result

    def detect_pattern(self, block_means):
        """Detect drift or equilibration patterns"""
        b1, b2, b3, b4 = block_means[1], block_means[2], block_means[3], block_means[4]

        if any(np.isnan(v) for v in [b1, b2, b3, b4]):
            return "INCOMPLETE"

        # Check for drift (Block4 very different from average of Block1-3)
        avg_123 = np.mean([b1, b2, b3])
        if abs(b4 - avg_123) > 5:
            return "DRIFT"

        # Check for equilibration (Block1 different, Block2-4 stable)
        avg_234 = np.mean([b2, b3, b4])
        std_234 = np.std([b2, b3, b4])
        if abs(b1 - avg_234) > 5 and std_234 < 3:
            return "EQUILIBRATION"

        # Check for oscillation
        if (b1 > b2 and b2 < b3 and b3 > b4) or (b1 < b2 and b2 > b3 and b3 < b4):
            return "OSCILLATION"

        return "STABLE"


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


def main():
    base = r'C:\Users\vasud\nod2-screening-data'
    output_dir = os.path.join(base, 'PHASE_A1_mutant_MD', 'analysis', 'convergence')
    os.makedirs(output_dir, exist_ok=True)

    print("=" * 100)
    print("CONVERGENCE ANALYSIS - BLOCK AVERAGING")
    print("=" * 100)
    print(f"Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"Convergence threshold: Range < {CONVERGENCE_THRESHOLD} kcal/mol")
    print("\n4 blocks: 0-25%, 25-50%, 50-75%, 75-100% of trajectory")

    analyzer = BlockAverageAnalyzer(base)
    all_sims = get_all_simulations()

    print("\n" + "=" * 100)
    print("ANALYZING ALL 25 SIMULATIONS")
    print("=" * 100 + "\n")

    for sim in all_sims:
        analyzer.analyze_simulation_blocks(sim, n_frames_per_block=25)

    results = analyzer.results

    # =========================================================================
    # TABLE 1: All Simulations
    # =========================================================================
    print("\n" + "=" * 120)
    print("TABLE 1: ALL SIMULATIONS - BLOCK AVERAGES")
    print("=" * 120)
    print(f"{'#':>2} | {'System':<6} | {'Ligand':<10} | {'Rep':<3} | {'Block1':>8} | {'Block2':>8} | "
          f"{'Block3':>8} | {'Block4':>8} | {'Range':>7} | {'Conv?':<5} | {'Pattern':<12}")
    print("-" * 120)

    for r in results:
        print(f"{r['num']:2d} | {r['system']:<6} | {r['ligand']:<10} | {r['rep']:<3} | "
              f"{r['block1']:>8.1f} | {r['block2']:>8.1f} | {r['block3']:>8.1f} | {r['block4']:>8.1f} | "
              f"{r['range']:>7.1f} | {r['converged']:<5} | {r['pattern']:<12}")

    # =========================================================================
    # TABLE 2: Summary by Condition
    # =========================================================================
    print("\n" + "=" * 100)
    print("TABLE 2: SUMMARY BY CONDITION")
    print("=" * 100)

    # Group by system + ligand
    conditions = {}
    for r in results:
        key = (r['system'], r['ligand'])
        if key not in conditions:
            conditions[key] = {'ranges': [], 'converged': []}
        conditions[key]['ranges'].append(r['range'])
        conditions[key]['converged'].append(r['converged'] == 'YES')

    print(f"\n{'System':<8} | {'Ligand':<10} | {'Avg Range':>10} | {'% Converged':>12} | {'Notes':<30}")
    print("-" * 80)

    condition_summary = []
    for (system, ligand), data in sorted(conditions.items()):
        avg_range = np.mean(data['ranges'])
        pct_converged = 100 * sum(data['converged']) / len(data['converged'])

        if pct_converged == 100:
            notes = "All replicates converged"
        elif pct_converged == 0:
            notes = "NO replicates converged"
        else:
            notes = f"{sum(data['converged'])}/{len(data['converged'])} converged"

        condition_summary.append({
            'system': system, 'ligand': ligand,
            'avg_range': avg_range, 'pct_converged': pct_converged,
            'notes': notes
        })

        print(f"{system:<8} | {ligand:<10} | {avg_range:>9.1f}  | {pct_converged:>11.0f}% | {notes:<30}")

    # =========================================================================
    # TABLE 3: Problem Simulations
    # =========================================================================
    print("\n" + "=" * 100)
    print("TABLE 3: PROBLEM SIMULATIONS (Range >= 6 kcal/mol)")
    print("=" * 100)

    problems = [r for r in results if r['converged'] == 'NO']

    if problems:
        print(f"\n{'Simulation':<25} | {'Range':>8} | {'Issue':<15} | {'Details':<40}")
        print("-" * 100)

        for r in problems:
            name = f"{r['system']}_{r['ligand']}_rep{r['rep']}"

            # Determine which block is problematic
            blocks = [r['block1'], r['block2'], r['block3'], r['block4']]
            mean_all = np.mean(blocks)

            if r['pattern'] == 'DRIFT':
                issue = "DRIFT"
                # Find which block drifted
                deviations = [abs(b - mean_all) for b in blocks]
                worst_block = deviations.index(max(deviations)) + 1
                details = f"Block {worst_block} deviates most ({blocks[worst_block-1]:.1f} vs avg {mean_all:.1f})"
            elif r['pattern'] == 'EQUILIBRATION':
                issue = "EQUILIBRATION"
                details = f"Block1={r['block1']:.1f}, avg(B2-4)={np.mean(blocks[1:]):.1f}"
            elif r['pattern'] == 'OSCILLATION':
                issue = "OSCILLATION"
                details = f"Blocks oscillate: {r['block1']:.0f}->{r['block2']:.0f}->{r['block3']:.0f}->{r['block4']:.0f}"
            else:
                issue = "HIGH VARIANCE"
                details = f"Range {r['range']:.1f} exceeds threshold"

            print(f"{name:<25} | {r['range']:>7.1f}  | {issue:<15} | {details:<40}")
    else:
        print("\n  No problem simulations detected - all converged!")

    # =========================================================================
    # PATTERN SUMMARY
    # =========================================================================
    print("\n" + "=" * 100)
    print("PATTERN DETECTION SUMMARY")
    print("=" * 100)

    patterns = Counter(r['pattern'] for r in results)
    print(f"\n  STABLE:        {patterns.get('STABLE', 0):2d} simulations")
    print(f"  DRIFT:         {patterns.get('DRIFT', 0):2d} simulations")
    print(f"  EQUILIBRATION: {patterns.get('EQUILIBRATION', 0):2d} simulations")
    print(f"  OSCILLATION:   {patterns.get('OSCILLATION', 0):2d} simulations")
    print(f"  INCOMPLETE:    {patterns.get('INCOMPLETE', 0):2d} simulations")

    # =========================================================================
    # SAVE RESULTS
    # =========================================================================
    df = pd.DataFrame(results)
    csv_path = os.path.join(output_dir, 'convergence_results.csv')
    df.to_csv(csv_path, index=False)
    print(f"\nSaved: {csv_path}")

    df_summary = pd.DataFrame(condition_summary)
    summary_path = os.path.join(output_dir, 'convergence_summary.csv')
    df_summary.to_csv(summary_path, index=False)
    print(f"Saved: {summary_path}")

    # =========================================================================
    # FINAL STATS
    # =========================================================================
    total_converged = sum(1 for r in results if r['converged'] == 'YES')
    pct_total = 100 * total_converged / len(results)

    print("\n" + "=" * 100)
    print("FINAL CONVERGENCE STATISTICS")
    print("=" * 100)
    print(f"\n  Total simulations: {len(results)}")
    print(f"  Converged:         {total_converged} ({pct_total:.0f}%)")
    print(f"  Not converged:     {len(results) - total_converged}")
    print(f"\n  Convergence criterion: Range < {CONVERGENCE_THRESHOLD} kcal/mol")

    print("\n" + "=" * 100)
    print(f"Completed: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("=" * 100)

    return results


if __name__ == '__main__':
    results = main()
