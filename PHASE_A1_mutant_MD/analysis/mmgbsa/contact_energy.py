#!/usr/bin/env python
"""
Contact-Based Binding Energy Estimator
Uses a knowledge-based scoring approach for robust relative comparisons

This calculates:
- Favorable contacts (hydrophobic-hydrophobic, H-bond donors/acceptors)
- Unfavorable contacts (charge clashes)
- Total interaction score

More robust than simplified MM-GBSA for relative comparisons.
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

        # Atom type classifications
        self.hydrophobic = {'C'}  # Non-polar
        self.hbond_donor = {'N', 'O'}  # Can donate H-bonds (when bonded to H)
        self.hbond_acceptor = {'N', 'O', 'S', 'F'}  # Can accept H-bonds
        self.charged_pos = {'N'}  # Potentially positive (Lys, Arg)
        self.charged_neg = {'O'}  # Potentially negative (Asp, Glu)

        # Contact energy parameters (empirical, in kcal/mol per contact)
        self.contact_energies = {
            'hydrophobic': -0.3,   # Favorable hydrophobic contact
            'hbond': -1.5,         # Hydrogen bond
            'ionic': -2.0,         # Salt bridge
            'clash': 0.5,          # Steric clash
            'neutral': -0.1,       # Non-specific contact
        }

    def classify_contact(self, lig_elem, prot_elem, distance):
        """Classify a contact between ligand and protein atoms"""
        if distance > 5.0:
            return None
        if distance < 1.5:
            return 'clash'

        # Hydrophobic contact
        if lig_elem in self.hydrophobic and prot_elem in self.hydrophobic:
            if distance < 4.5:
                return 'hydrophobic'

        # H-bond (simplified: any N/O to N/O within 3.5 Å)
        if lig_elem in self.hbond_donor and prot_elem in self.hbond_acceptor:
            if distance < 3.5:
                return 'hbond'
        if lig_elem in self.hbond_acceptor and prot_elem in self.hbond_donor:
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

        # Get coordinates with PBC
        box = universe.trajectory.ts.dimensions
        dist_matrix = distances.distance_array(
            ligand.positions, protein.positions, box=box
        )

        # Get elements
        lig_elems = [a.element if hasattr(a, 'element') and a.element else a.name[0]
                     for a in ligand]
        prot_elems = [a.element if hasattr(a, 'element') and a.element else a.name[0]
                      for a in protein]

        # Count and classify contacts
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
            'clash': contact_counts.get('clash', 0),
            'neutral': contact_counts.get('neutral', 0),
        }

    def analyze_simulation(self, sim_info, n_frames=30):
        """Analyze a single simulation"""
        print(f"\n{'='*60}")
        print(f"Analyzing: {sim_info['system']}_{sim_info['ligand']}_rep{sim_info['rep']}")
        print(f"{'='*60}")

        pdb_path = os.path.join(self.base_dir, sim_info['pdb'])
        dcd_path = os.path.join(self.base_dir, sim_info['dcd'])

        # Load trajectory
        print(f"  Loading trajectory...")
        u = mda.Universe(pdb_path, dcd_path)
        total_frames = len(u.trajectory)

        # Analyze first 20% (equilibrated but still bound)
        start_frame = int(total_frames * 0.05)
        end_frame = int(total_frames * 0.25)
        stride = max(1, (end_frame - start_frame) // n_frames)

        print(f"    Frames {start_frame}-{end_frame}, stride={stride}")

        # Check ligand
        ligand = u.select_atoms('resname UNK')
        print(f"  Ligand atoms: {len(ligand)}")

        if len(ligand) == 0:
            print("  ERROR: No ligand found!")
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

            if frame_count % 10 == 0:
                print(f"  Frame {frame_count}: E={E['energy']:.1f}, contacts={E['contacts']}, H-bonds={E['hbonds']}")

        # Compile results
        energy_vals = [e['energy'] for e in energies if not np.isnan(e['energy'])]
        contact_vals = [e['contacts'] for e in energies]
        hbond_vals = [e['hbonds'] for e in energies]

        result = {
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

        print(f"\n  Results:")
        print(f"    Contact Energy: {result['DeltaE_mean']:.1f} ± {result['DeltaE_std']:.1f} kcal/mol")
        print(f"    Contacts: {result['contacts_mean']:.0f}, H-bonds: {result['hbonds_mean']:.1f}")

        # Sanity check
        if result['DeltaE_mean'] > 0:
            result['flag'] = 'POSITIVE'
        elif result['contacts_mean'] < 10:
            result['flag'] = 'LOW_CONTACTS'
        else:
            result['flag'] = 'OK'

        self.results.append(result)
        return result


def run_pilot_test():
    """Run pilot test on 3 simulations"""
    base = r'C:\Users\vasud\nod2-screening-data'
    calc = ContactEnergyCalculator(base)

    # Pilot simulations
    pilots = [
        {'num': 1, 'phase': 'Phase6', 'system': 'WT', 'ligand': 'Febuxostat', 'rep': 1,
         'dcd': os.path.join('PHASE_6', 'trajectories', 'febuxostat_rep1.dcd'),
         'pdb': os.path.join('PHASE_6', 'trajectories', 'febuxostat_rep1_solvated.pdb'),
         'note': 'Lead compound'},
        {'num': 2, 'phase': 'Phase6', 'system': 'WT', 'ligand': 'Natural', 'rep': 1,
         'dcd': os.path.join('PHASE_6', 'trajectories', 'natural_top_rep1.dcd'),
         'pdb': os.path.join('PHASE_6', 'trajectories', 'natural_top_rep1_solvated.pdb'),
         'note': 'Backup compound'},
        {'num': 3, 'phase': 'Phase6', 'system': 'WT', 'ligand': 'Decoy', 'rep': 1,
         'dcd': os.path.join('PHASE_6', 'trajectories', 'decoy_rep1.dcd'),
         'pdb': os.path.join('PHASE_6', 'trajectories', 'decoy_rep1_solvated.pdb'),
         'note': 'NEGATIVE CONTROL'},
    ]

    print("="*70)
    print("CONTACT-BASED BINDING ENERGY - PILOT TEST")
    print("="*70)
    print(f"Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")

    for sim in pilots:
        calc.analyze_simulation(sim, n_frames=30)

    # Summary
    print("\n" + "="*70)
    print("PILOT TEST SUMMARY")
    print("="*70)
    print(f"\n{'System':<25} {'Energy (kcal/mol)':<18} {'Contacts':<10} {'H-bonds':<10} {'Flag'}")
    print("-"*70)

    for r in calc.results:
        name = f"{r['system']}_{r['ligand']}_rep{r['rep']}"
        energy = f"{r['DeltaE_mean']:.1f} ± {r['DeltaE_std']:.1f}"
        print(f"{name:<25} {energy:<18} {r['contacts_mean']:<10.0f} {r['hbonds_mean']:<10.1f} {r['flag']}")

    print("\n" + "="*70)
    print("Expected ranking: Febuxostat ~ Natural >> Decoy")
    print("Lower (more negative) energy = better binding")
    print("="*70)

    return calc.results


if __name__ == '__main__':
    results = run_pilot_test()
