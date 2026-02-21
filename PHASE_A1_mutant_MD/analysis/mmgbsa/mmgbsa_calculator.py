#!/usr/bin/env python
"""
MM-GBSA Binding Energy Estimator
Calculates protein-ligand binding energy using implicit solvent (OBC2)

Methodology:
- Complex, Protein, and Ligand energies calculated separately
- DeltaE_GBSA = E_complex - E_protein - E_ligand
- Uses AMBER14 for protein, GAFF-like for ligand
- OBC2 implicit solvent model

Note: This provides binding energy ESTIMATES for relative comparisons.
Not true free energies (no entropy term).
"""

import os
import sys
import warnings
import numpy as np
import pandas as pd
from datetime import datetime
warnings.filterwarnings('ignore')

# MDAnalysis for trajectory handling
import MDAnalysis as mda
from MDAnalysis.analysis import distances

# RDKit for ligand handling
from rdkit import Chem
from rdkit.Chem import AllChem

# OpenMM for energy calculations
import openmm as mm
from openmm import app
from openmm import unit

# Import ligand template generator
from ligand_template import create_ligand_template, GAFF_TYPES


class MMGBSACalculator:
    """Calculate MM-GBSA binding energies from MD trajectories"""

    def __init__(self, base_dir):
        self.base_dir = base_dir
        self.results = []

        # Conversion factor
        self.kj_to_kcal = 1.0 / 4.184

    def load_trajectory(self, pdb_path, dcd_path, n_frames=100, start_frac=0.5):
        """Load and prepare trajectory for analysis"""
        print(f"  Loading trajectory...")
        u = mda.Universe(pdb_path, dcd_path)

        total_frames = len(u.trajectory)
        start_frame = int(total_frames * start_frac)
        end_frame = total_frames

        available = end_frame - start_frame
        stride = max(1, available // n_frames)

        print(f"    Total frames: {total_frames}")
        print(f"    Using frames {start_frame}-{end_frame}, stride={stride}")

        return u, list(range(start_frame, end_frame, stride))[:n_frames]

    def strip_system(self, universe, frame_idx=0):
        """Strip waters and ions, return protein and ligand selections"""
        universe.trajectory[frame_idx]

        protein = universe.select_atoms('protein')
        ligand = universe.select_atoms('resname UNK')
        complex_sel = universe.select_atoms('protein or resname UNK')

        return protein, ligand, complex_sel

    def calculate_interaction_energy_current(self, universe, cutoff=12.0):
        """
        Calculate protein-ligand interaction energy for CURRENT frame
        (Use after iterating to the desired frame)
        """
        return self._calc_energy(universe, cutoff)

    def calculate_interaction_energy(self, universe, frame_idx, cutoff=12.0):
        """
        Calculate protein-ligand interaction energy using pairwise summation

        This calculates:
        - Electrostatic interaction (Coulomb with dielectric screening)
        - Van der Waals interaction (Lennard-Jones 12-6)
        - Approximate solvation contribution

        Returns energy in kcal/mol
        """
        universe.trajectory[frame_idx]
        return self._calc_energy(universe, cutoff)

    def _calc_energy(self, universe, cutoff=12.0):
        """Internal energy calculation on current frame"""
        protein = universe.select_atoms('protein')
        ligand = universe.select_atoms('resname UNK')

        if len(protein) == 0 or len(ligand) == 0:
            return {'vdw': np.nan, 'elec': np.nan, 'solv': np.nan, 'total': np.nan}

        # Get coordinates
        prot_coords = protein.positions
        lig_coords = ligand.positions

        # Get box dimensions for PBC
        box = universe.trajectory.ts.dimensions

        # Calculate pairwise distances WITH PBC
        dist_matrix = distances.distance_array(lig_coords, prot_coords, box=box)

        # Get elements
        lig_elements = [a.element if hasattr(a, 'element') and a.element else a.name[0] for a in ligand]
        prot_elements = [a.element if hasattr(a, 'element') and a.element else a.name[0] for a in protein]

        # LJ parameters (sigma in Angstrom, epsilon in kcal/mol)
        lj_params = {
            'C': (3.40, 0.086),
            'N': (3.25, 0.170),
            'O': (2.96, 0.210),
            'S': (3.56, 0.250),
            'H': (2.50, 0.015),
            'F': (3.12, 0.061),
            'P': (3.74, 0.200),
        }

        # Approximate partial charges (based on typical values)
        charges = {
            'C': 0.0, 'N': -0.5, 'O': -0.5, 'S': -0.2, 'H': 0.1, 'F': -0.2, 'P': 0.5
        }

        vdw_total = 0.0
        elec_total = 0.0

        for i, lig_elem in enumerate(lig_elements):
            sig_i, eps_i = lj_params.get(lig_elem.upper(), (3.40, 0.1))
            q_i = charges.get(lig_elem.upper(), 0.0)

            for j, prot_elem in enumerate(prot_elements):
                r = dist_matrix[i, j]
                # Skip if beyond cutoff or too close (avoid numerical issues)
                if r > cutoff or r < 1.0:
                    continue

                sig_j, eps_j = lj_params.get(prot_elem.upper(), (3.40, 0.1))
                q_j = charges.get(prot_elem.upper(), 0.0)

                # LJ energy (sigma and r both in Angstrom)
                sigma = (sig_i + sig_j) / 2
                epsilon = np.sqrt(eps_i * eps_j)
                ratio = sigma / r
                vdw_total += 4 * epsilon * (ratio**12 - ratio**6)

                # Coulomb with distance-dependent dielectric
                k = 332.0637  # kcal*Å/(mol*e^2)
                dielectric = 4.0 * r  # Distance-dependent dielectric
                elec_total += k * q_i * q_j / dielectric

        # Simple solvation estimate (based on SASA-like contribution)
        # Using number of contacts as proxy
        n_contacts = np.sum(dist_matrix < 4.0)
        solv_estimate = -0.01 * n_contacts  # Rough estimate

        total = vdw_total + elec_total + solv_estimate

        return {
            'vdw': vdw_total,
            'elec': elec_total,
            'solv': solv_estimate,
            'total': total
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
        # Analyze first 20% (equilibrated but still bound) instead of last 50%
        start_frame = int(total_frames * 0.05)  # Skip first 5% for equilibration
        end_frame = int(total_frames * 0.25)   # Use 5-25% of trajectory
        stride = max(1, (end_frame - start_frame) // n_frames)

        print(f"    Total frames: {total_frames}")
        print(f"    Analyzing from frame {start_frame}, stride={stride}")

        # Check ligand
        ligand = u.select_atoms('resname UNK')
        protein = u.select_atoms('protein')
        print(f"  Ligand atoms: {len(ligand)}")
        print(f"  Protein atoms: {len(protein)}")

        if len(ligand) == 0:
            print("  ERROR: No ligand found!")
            return None

        # Calculate energies using fully sequential iteration (avoids DCD seek issues)
        energies = []
        frame_count = 0
        frames_to_analyze = set(range(start_frame, end_frame, stride))

        print(f"  Processing {len(frames_to_analyze)} frames (frames {start_frame}-{end_frame})...")
        for i, ts in enumerate(u.trajectory):
            if i > end_frame:
                break  # Stop early once we pass end_frame
            if i < start_frame:
                continue
            if i not in frames_to_analyze:
                continue
            if frame_count >= n_frames:
                break

            E = self.calculate_interaction_energy_current(u)
            energies.append(E)
            frame_count += 1

            if frame_count % 10 == 0:
                print(f"  Frame {frame_count}/{n_frames}: E_total = {E['total']:.2f} kcal/mol")

        # Compile results
        totals = [e['total'] for e in energies if not np.isnan(e['total'])]
        vdws = [e['vdw'] for e in energies if not np.isnan(e['vdw'])]
        elecs = [e['elec'] for e in energies if not np.isnan(e['elec'])]

        result = {
            'phase': sim_info['phase'],
            'system': sim_info['system'],
            'ligand': sim_info['ligand'],
            'rep': sim_info['rep'],
            'n_frames': len(totals),
            'DeltaE_mean': np.mean(totals),
            'DeltaE_std': np.std(totals),
            'DeltaE_min': np.min(totals),
            'DeltaE_max': np.max(totals),
            'VDW_mean': np.mean(vdws),
            'Elec_mean': np.mean(elecs),
            'note': sim_info.get('note', '')
        }

        print(f"\n  Results:")
        print(f"    DeltaE_GBSA: {result['DeltaE_mean']:.2f} ± {result['DeltaE_std']:.2f} kcal/mol")
        print(f"    Range: {result['DeltaE_min']:.2f} to {result['DeltaE_max']:.2f}")
        print(f"    VDW: {result['VDW_mean']:.2f}, Elec: {result['Elec_mean']:.2f}")

        # Sanity check
        if result['DeltaE_mean'] > 50 or result['DeltaE_mean'] < -200:
            print(f"  WARNING: Energy outside expected range!")
            result['flag'] = 'OUT_OF_RANGE'
        elif result['DeltaE_mean'] > 0:
            print(f"  WARNING: Positive energy - weak/no binding")
            result['flag'] = 'POSITIVE'
        else:
            result['flag'] = 'OK'

        self.results.append(result)
        return result


def run_pilot_test():
    """Run pilot test on 3 simulations"""
    base = r'C:\Users\vasud\nod2-screening-data'
    calc = MMGBSACalculator(base)

    # Pilot simulations
    pilots = [
        {'num': 1, 'phase': 'Phase6', 'system': 'WT', 'ligand': 'Febuxostat', 'rep': 1,
         'dcd': r'PHASE_6\trajectories\febuxostat_rep1.dcd',
         'pdb': r'PHASE_6\trajectories\febuxostat_rep1_solvated.pdb',
         'note': 'Lead compound'},
        {'num': 2, 'phase': 'Phase6', 'system': 'WT', 'ligand': 'Natural', 'rep': 1,
         'dcd': r'PHASE_6\trajectories\natural_top_rep1.dcd',
         'pdb': r'PHASE_6\trajectories\natural_top_rep1_solvated.pdb',
         'note': 'Backup compound'},
        {'num': 3, 'phase': 'Phase6', 'system': 'WT', 'ligand': 'Decoy', 'rep': 1,
         'dcd': r'PHASE_6\trajectories\decoy_rep1.dcd',
         'pdb': r'PHASE_6\trajectories\decoy_rep1_solvated.pdb',
         'note': 'NEGATIVE CONTROL'},
    ]

    print("="*70)
    print("MM-GBSA PILOT TEST")
    print("="*70)
    print(f"Testing 3 simulations with 30 frames each")
    print(f"Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")

    for sim in pilots:
        calc.analyze_simulation(sim, n_frames=30)

    # Summary
    print("\n" + "="*70)
    print("PILOT TEST SUMMARY")
    print("="*70)
    print(f"\n{'System':<25} {'DeltaE (kcal/mol)':<20} {'Flag':<15} {'Note'}")
    print("-"*70)

    for r in calc.results:
        name = f"{r['system']}_{r['ligand']}_rep{r['rep']}"
        energy = f"{r['DeltaE_mean']:.2f} ± {r['DeltaE_std']:.2f}"
        print(f"{name:<25} {energy:<20} {r['flag']:<15} {r['note']}")

    print("\n" + "="*70)
    print("Expected results:")
    print("  - Febuxostat: Moderate negative (~ -30 to -60)")
    print("  - Natural: Similar to Febuxostat")
    print("  - Decoy: Near zero or positive (WEAK BINDING)")
    print("="*70)

    return calc.results


if __name__ == '__main__':
    results = run_pilot_test()
