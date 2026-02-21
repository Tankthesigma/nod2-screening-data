#!/usr/bin/env python
"""
Run Contact-Based Binding Energy Analysis on ALL 25 Simulations
Phase 6 (WT): 13 simulations
Phase A1 (Mutants): 12 simulations
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
from collections import Counter


class ContactEnergyCalculator:
    """Calculate contact-based binding energies from MD trajectories"""

    def __init__(self, base_dir):
        self.base_dir = base_dir
        self.results = []

        # Contact energy parameters (empirical, in kcal/mol per contact)
        self.contact_energies = {
            'hydrophobic': -0.3,
            'hbond': -1.5,
            'ionic': -2.0,
            'clash': 0.5,
            'neutral': -0.1,
        }

    def classify_contact(self, lig_elem, prot_elem, distance):
        """Classify a contact between ligand and protein atoms"""
        if distance > 5.0:
            return None
        if distance < 1.5:
            return 'clash'

        # Hydrophobic contact
        if lig_elem == 'C' and prot_elem == 'C':
            if distance < 4.5:
                return 'hydrophobic'

        # H-bond (simplified)
        if lig_elem in ['N', 'O'] and prot_elem in ['N', 'O']:
            if distance < 3.5:
                return 'hbond'

        # Generic contact
        if distance < 4.0:
            return 'neutral'

        return None

    def calculate_contact_energy(self, universe, cutoff=6.0):
        """Calculate contact-based energy for current frame"""
        protein = universe.select_atoms('protein')
        ligand = universe.select_atoms('resname UNK')

        if len(protein) == 0 or len(ligand) == 0:
            return {'energy': np.nan, 'contacts': 0, 'hbonds': 0}

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
            'hbonds': contact_counts.get('hbond', 0),
            'hydrophobic': contact_counts.get('hydrophobic', 0),
        }

    def analyze_simulation(self, sim_info, n_frames=50):
        """Analyze a single simulation"""
        name = f"{sim_info['system']}_{sim_info['ligand']}_rep{sim_info['rep']}"
        print(f"\n[{sim_info['num']:2d}/25] {name}")

        pdb_path = os.path.join(self.base_dir, sim_info['pdb'])
        dcd_path = os.path.join(self.base_dir, sim_info['dcd'])

        if not os.path.exists(dcd_path):
            print(f"  ERROR: DCD not found")
            return None
        if not os.path.exists(pdb_path):
            print(f"  ERROR: PDB not found")
            return None

        try:
            u = mda.Universe(pdb_path, dcd_path)
            total_frames = len(u.trajectory)
        except Exception as e:
            print(f"  ERROR loading: {e}")
            return None

        # Analyze first 20% (equilibrated but still bound)
        start_frame = int(total_frames * 0.05)
        end_frame = int(total_frames * 0.25)
        stride = max(1, (end_frame - start_frame) // n_frames)

        ligand = u.select_atoms('resname UNK')
        if len(ligand) == 0:
            print(f"  ERROR: No ligand found")
            return None

        # Calculate energies
        energies = []
        frame_count = 0
        frames_to_analyze = set(range(start_frame, end_frame, stride))

        for i, ts in enumerate(u.trajectory):
            if i > end_frame:
                break
            if i < start_frame or i not in frames_to_analyze:
                continue
            if frame_count >= n_frames:
                break

            E = self.calculate_contact_energy(u)
            energies.append(E)
            frame_count += 1

        if not energies:
            print(f"  ERROR: No frames analyzed")
            return None

        energy_vals = [e['energy'] for e in energies if not np.isnan(e['energy'])]
        contact_vals = [e['contacts'] for e in energies]
        hbond_vals = [e['hbonds'] for e in energies]

        result = {
            'num': sim_info['num'],
            'phase': sim_info['phase'],
            'system': sim_info['system'],
            'ligand': sim_info['ligand'],
            'rep': sim_info['rep'],
            'n_frames': len(energy_vals),
            'DeltaE_mean': np.mean(energy_vals),
            'DeltaE_std': np.std(energy_vals),
            'DeltaE_min': np.min(energy_vals),
            'DeltaE_max': np.max(energy_vals),
            'contacts_mean': np.mean(contact_vals),
            'hbonds_mean': np.mean(hbond_vals),
            'note': sim_info.get('note', '')
        }

        print(f"  E={result['DeltaE_mean']:.1f} +/- {result['DeltaE_std']:.1f} kcal/mol, "
              f"contacts={result['contacts_mean']:.0f}, H-bonds={result['hbonds_mean']:.1f}")

        self.results.append(result)
        return result


def get_all_simulations():
    """Define all 25 simulations"""
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

    # Phase A1 - Mutants (12 sims)
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
         os.path.join('PHASE_A1_mutant_MD', 'trajectories', 'G908R_febuxostat_rep1_solvated.pdb'), 'OUTLIER 45.8%'),
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
    """Run all analyses"""
    base = r'C:\Users\vasud\nod2-screening-data'
    calc = ContactEnergyCalculator(base)

    all_sims = get_all_simulations()

    print("="*70)
    print("CONTACT-BASED BINDING ENERGY ANALYSIS")
    print("="*70)
    print(f"Total simulations: {len(all_sims)}")
    print(f"Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")

    # Run all simulations
    for sim in all_sims:
        try:
            calc.analyze_simulation(sim, n_frames=50)
        except Exception as e:
            print(f"  FAILED: {e}")

    # Create results DataFrame
    df = pd.DataFrame(calc.results)

    # Save raw results
    output_dir = os.path.join(base, 'PHASE_A1_mutant_MD', 'analysis', 'mmgbsa')
    os.makedirs(output_dir, exist_ok=True)

    df.to_csv(os.path.join(output_dir, 'contact_energy_results.csv'), index=False)
    print(f"\nSaved: contact_energy_results.csv")

    # Print summary tables
    print("\n" + "="*70)
    print("PHASE 6 (WILD-TYPE) RESULTS")
    print("="*70)

    wt_df = df[df['phase'] == 'Phase6']
    print(f"\n{'System':<25} {'Energy (kcal/mol)':<18} {'Contacts':<10} {'H-bonds':<10}")
    print("-"*65)
    for _, r in wt_df.iterrows():
        name = f"WT_{r['ligand']}_rep{r['rep']}"
        energy = f"{r['DeltaE_mean']:.1f} +/- {r['DeltaE_std']:.1f}"
        print(f"{name:<25} {energy:<18} {r['contacts_mean']:<10.0f} {r['hbonds_mean']:<10.1f}")

    # WT Summary by ligand
    print("\n" + "-"*65)
    print("WT SUMMARY BY LIGAND:")
    for ligand in ['Febuxostat', 'Natural', 'Budesonide', 'Ursodiol', 'Decoy']:
        lig_df = wt_df[wt_df['ligand'] == ligand]
        if len(lig_df) > 0:
            mean_e = lig_df['DeltaE_mean'].mean()
            std_e = lig_df['DeltaE_mean'].std()
            mean_c = lig_df['contacts_mean'].mean()
            print(f"  {ligand:<12}: {mean_e:6.1f} +/- {std_e:4.1f} kcal/mol, {mean_c:5.0f} contacts")

    print("\n" + "="*70)
    print("PHASE A1 (MUTANT) RESULTS")
    print("="*70)

    mut_df = df[df['phase'] == 'PhaseA1']
    print(f"\n{'System':<30} {'Energy (kcal/mol)':<18} {'Contacts':<10} {'H-bonds':<10}")
    print("-"*70)
    for _, r in mut_df.iterrows():
        name = f"{r['system']}_{r['ligand']}_rep{r['rep']}"
        energy = f"{r['DeltaE_mean']:.1f} +/- {r['DeltaE_std']:.1f}"
        note = f" *{r['note']}" if r['note'] else ""
        print(f"{name:<30} {energy:<18} {r['contacts_mean']:<10.0f} {r['hbonds_mean']:<10.1f}{note}")

    # Mutation effect summary
    print("\n" + "-"*70)
    print("MUTATION EFFECT SUMMARY (vs WT):")

    for ligand in ['Febuxostat', 'Natural']:
        wt_mean = wt_df[wt_df['ligand'] == ligand]['DeltaE_mean'].mean()
        for mut in ['R702W', 'G908R']:
            mut_mean = mut_df[(mut_df['system'] == mut) & (mut_df['ligand'] == ligand)]['DeltaE_mean'].mean()
            delta = mut_mean - wt_mean
            print(f"  {mut} + {ligand}: DeltaDeltaE = {delta:+.1f} kcal/mol")

    print("\n" + "="*70)
    print(f"Completed: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"Total successful: {len(df)}/25")
    print("="*70)

    return df


if __name__ == '__main__':
    df = main()
