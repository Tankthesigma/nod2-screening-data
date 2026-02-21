#!/usr/bin/env python
"""
Simplified Protein-Ligand Interaction Energy Calculator
Uses MDAnalysis for trajectory handling and basic energy calculations

This calculates:
1. Electrostatic interaction energy (Coulomb)
2. Van der Waals interaction energy (Lennard-Jones approximation)
3. Combined interaction score

Note: This is a simplified approach for RELATIVE comparisons.
Not a full MM-GBSA calculation but provides useful binding insights.
"""

import MDAnalysis as mda
from MDAnalysis.analysis import distances
import numpy as np
import warnings
warnings.filterwarnings('ignore')

# Standard AMBER-like parameters for common atom types
# LJ parameters (sigma in Angstrom, epsilon in kcal/mol)
LJ_PARAMS = {
    'C': (1.908, 0.086),
    'N': (1.824, 0.170),
    'O': (1.661, 0.210),
    'S': (2.000, 0.250),
    'H': (1.100, 0.016),
    'F': (1.750, 0.061),
    'Cl': (1.948, 0.265),
    'Br': (2.020, 0.420),
    'I': (2.150, 0.500),
    'P': (2.100, 0.200),
}

# Default partial charges for common elements
DEFAULT_CHARGES = {
    'C': 0.0,
    'N': -0.5,
    'O': -0.5,
    'S': -0.2,
    'H': 0.1,
    'F': -0.2,
    'Cl': -0.1,
    'Br': -0.1,
    'I': -0.1,
    'P': 0.5,
}

def get_lj_params(element):
    """Get LJ parameters for an element"""
    elem = element.upper().strip()
    if elem in LJ_PARAMS:
        return LJ_PARAMS[elem]
    return (1.9, 0.1)  # Default

def get_charge(element):
    """Get approximate charge for an element"""
    elem = element.upper().strip()
    if elem in DEFAULT_CHARGES:
        return DEFAULT_CHARGES[elem]
    return 0.0

def calc_lj_energy(r, sigma_i, epsilon_i, sigma_j, epsilon_j):
    """Calculate Lennard-Jones energy between two atoms"""
    sigma = (sigma_i + sigma_j) / 2
    epsilon = np.sqrt(epsilon_i * epsilon_j)

    # Avoid division by zero
    r = np.maximum(r, 0.1)

    # Standard 12-6 LJ
    ratio = sigma / r
    return 4 * epsilon * (ratio**12 - ratio**6)

def calc_coulomb_energy(r, q_i, q_j, dielectric=4.0):
    """Calculate Coulomb energy between two atoms"""
    # k = 332.0637 kcal*Å/(mol*e^2)
    k = 332.0637
    r = np.maximum(r, 0.1)
    return k * q_i * q_j / (dielectric * r)

def calculate_interaction_energy(universe, protein_sel='protein', ligand_sel='resname UNK',
                                  cutoff=12.0):
    """
    Calculate protein-ligand interaction energy for current frame

    Parameters:
    -----------
    universe : MDAnalysis.Universe
        Loaded trajectory/topology
    protein_sel : str
        Selection string for protein
    ligand_sel : str
        Selection string for ligand
    cutoff : float
        Distance cutoff in Angstrom

    Returns:
    --------
    dict : Energy components in kcal/mol
    """
    protein = universe.select_atoms(protein_sel)
    ligand = universe.select_atoms(ligand_sel)

    if len(protein) == 0 or len(ligand) == 0:
        return {'vdw': np.nan, 'elec': np.nan, 'total': np.nan}

    # Get coordinates
    prot_coords = protein.positions
    lig_coords = ligand.positions

    # Calculate all pairwise distances
    dist_matrix = distances.distance_array(lig_coords, prot_coords)

    # Get atom info
    lig_elements = [a.element if hasattr(a, 'element') and a.element else a.name[0] for a in ligand]
    prot_elements = [a.element if hasattr(a, 'element') and a.element else a.name[0] for a in protein]

    # Calculate energies
    vdw_total = 0.0
    elec_total = 0.0

    for i, lig_elem in enumerate(lig_elements):
        sig_i, eps_i = get_lj_params(lig_elem)
        q_i = get_charge(lig_elem)

        for j, prot_elem in enumerate(prot_elements):
            r = dist_matrix[i, j]
            if r > cutoff:
                continue

            sig_j, eps_j = get_lj_params(prot_elem)
            q_j = get_charge(prot_elem)

            # LJ energy
            vdw_total += calc_lj_energy(r, sig_i, eps_i, sig_j, eps_j)

            # Coulomb energy
            elec_total += calc_coulomb_energy(r, q_i, q_j)

    return {
        'vdw': vdw_total,
        'elec': elec_total,
        'total': vdw_total + elec_total
    }


def analyze_trajectory(pdb_path, dcd_path, n_frames=30, start_frac=0.5):
    """
    Analyze trajectory and calculate interaction energies

    Parameters:
    -----------
    pdb_path : str
        Path to topology PDB
    dcd_path : str
        Path to trajectory DCD
    n_frames : int
        Number of frames to analyze
    start_frac : float
        Start analysis from this fraction of trajectory

    Returns:
    --------
    dict : Energy statistics
    """
    print(f"Loading trajectory...")
    u = mda.Universe(pdb_path, dcd_path)

    total_frames = len(u.trajectory)
    start_frame = int(total_frames * start_frac)
    end_frame = total_frames

    # Calculate stride
    available = end_frame - start_frame
    stride = max(1, available // n_frames)

    print(f"  Total frames: {total_frames}")
    print(f"  Analyzing frames {start_frame} to {end_frame} (stride={stride})")

    # Check ligand
    ligand = u.select_atoms('resname UNK')
    protein = u.select_atoms('protein')
    print(f"  Ligand atoms: {len(ligand)}")
    print(f"  Protein atoms: {len(protein)}")

    if len(ligand) == 0:
        print("ERROR: No ligand found!")
        return None

    # Calculate energies
    energies = {'vdw': [], 'elec': [], 'total': []}

    frame_indices = list(range(start_frame, end_frame, stride))[:n_frames]

    for i, frame_idx in enumerate(frame_indices):
        u.trajectory[frame_idx]
        E = calculate_interaction_energy(u)

        for key in energies:
            energies[key].append(E[key])

        if (i + 1) % 10 == 0:
            print(f"  Processed {i+1}/{len(frame_indices)} frames...")

    # Calculate statistics
    results = {}
    for key in energies:
        vals = np.array(energies[key])
        results[f'{key}_mean'] = np.nanmean(vals)
        results[f'{key}_std'] = np.nanstd(vals)
        results[f'{key}_min'] = np.nanmin(vals)
        results[f'{key}_max'] = np.nanmax(vals)

    return results


if __name__ == '__main__':
    import os

    base = r'C:\Users\vasud\nod2-screening-data'

    # Test on febuxostat_rep1
    pdb = os.path.join(base, r'PHASE_6\trajectories\febuxostat_rep1_solvated.pdb')
    dcd = os.path.join(base, r'PHASE_6\trajectories\febuxostat_rep1.dcd')

    print("="*60)
    print("Simplified Interaction Energy Analysis")
    print("="*60)
    print(f"\nAnalyzing: WT_Febuxostat_rep1")

    results = analyze_trajectory(pdb, dcd, n_frames=30)

    if results:
        print(f"\n{'='*60}")
        print("RESULTS:")
        print(f"{'='*60}")
        print(f"  VDW Energy:   {results['vdw_mean']:8.2f} ± {results['vdw_std']:6.2f} kcal/mol")
        print(f"  Elec Energy:  {results['elec_mean']:8.2f} ± {results['elec_std']:6.2f} kcal/mol")
        print(f"  Total:        {results['total_mean']:8.2f} ± {results['total_std']:6.2f} kcal/mol")
        print(f"\n  Range: {results['total_min']:.2f} to {results['total_max']:.2f}")
