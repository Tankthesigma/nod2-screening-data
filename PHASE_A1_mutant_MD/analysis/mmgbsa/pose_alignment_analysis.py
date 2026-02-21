#!/usr/bin/env python
"""
Ligand Pose Alignment Analysis
Compare binding modes of Febuxostat vs Natural compound
"""

import os
import sys
import warnings
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from datetime import datetime
warnings.filterwarnings('ignore')

import MDAnalysis as mda
from MDAnalysis.analysis import distances, align


def kabsch_align(mobile, reference):
    """
    Align mobile coordinates to reference using Kabsch algorithm
    Returns: aligned coordinates, rotation matrix, RMSD
    """
    # Center both
    mobile_center = mobile.mean(axis=0)
    ref_center = reference.mean(axis=0)

    mobile_centered = mobile - mobile_center
    ref_centered = reference - ref_center

    # Compute optimal rotation
    H = mobile_centered.T @ ref_centered
    U, S, Vt = np.linalg.svd(H)
    R = Vt.T @ U.T

    # Handle reflection
    if np.linalg.det(R) < 0:
        Vt[-1, :] *= -1
        R = Vt.T @ U.T

    # Apply rotation and translation
    aligned = (mobile_centered @ R) + ref_center

    # Calculate RMSD
    rmsd = np.sqrt(np.mean(np.sum((aligned - reference) ** 2, axis=1)))

    return aligned, R, rmsd, ref_center, mobile_center


def get_contacts(universe, ligand_sel='resname UNK', protein_sel='protein', cutoff=4.0):
    """
    Get protein residues in contact with ligand
    Returns: list of (resid, resname, min_distance)
    """
    ligand = universe.select_atoms(ligand_sel)
    protein = universe.select_atoms(protein_sel)

    if len(ligand) == 0 or len(protein) == 0:
        return []

    box = universe.trajectory.ts.dimensions
    dist_matrix = distances.distance_array(ligand.positions, protein.positions, box=box)

    contacts = {}
    for j, atom in enumerate(protein):
        min_dist = np.min(dist_matrix[:, j])
        if min_dist < cutoff:
            resid = atom.resid
            resname = atom.resname
            key = (resid, resname)
            if key not in contacts or contacts[key] > min_dist:
                contacts[key] = min_dist

    return [(k[0], k[1], v) for k, v in sorted(contacts.items())]


def load_representative_frame(pdb_path, dcd_path, frame_frac=0.15):
    """
    Load a representative frame from trajectory
    """
    u = mda.Universe(pdb_path, dcd_path)
    total_frames = len(u.trajectory)
    target_frame = int(total_frames * frame_frac)
    u.trajectory[target_frame]
    return u


def main():
    base = r'C:\Users\vasud\nod2-screening-data'
    output_dir = os.path.join(base, 'PHASE_A1_mutant_MD', 'analysis', 'mmgbsa')

    print("=" * 80)
    print("LIGAND POSE ALIGNMENT ANALYSIS")
    print("Febuxostat vs Natural Compound Binding Mode Comparison")
    print("=" * 80)
    print(f"Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")

    # Define simulations
    simulations = {
        'Febuxostat': {
            'pdb': os.path.join(base, 'PHASE_6', 'trajectories', 'febuxostat_rep2_solvated.pdb'),
            'dcd': os.path.join(base, 'PHASE_6', 'trajectories', 'febuxostat_rep2.dcd'),
            'color': 'blue',
            'label': 'Febuxostat'
        },
        'Natural': {
            'pdb': os.path.join(base, 'PHASE_6', 'trajectories', 'natural_top_rep1_solvated.pdb'),
            'dcd': os.path.join(base, 'PHASE_6', 'trajectories', 'natural_top_rep1.dcd'),
            'color': 'green',
            'label': 'Natural Compound'
        }
    }

    # Load representative frames
    print("\n" + "=" * 60)
    print("LOADING REPRESENTATIVE FRAMES")
    print("=" * 60)

    universes = {}
    for name, sim in simulations.items():
        print(f"\nLoading {name}...")
        u = load_representative_frame(sim['pdb'], sim['dcd'], frame_frac=0.15)
        universes[name] = u

        ligand = u.select_atoms('resname UNK')
        protein = u.select_atoms('protein and name CA')
        print(f"  Ligand atoms: {len(ligand)}")
        print(f"  Protein CA atoms: {len(protein)}")

    # =========================================================================
    # ALIGN PROTEINS
    # =========================================================================
    print("\n" + "=" * 60)
    print("ALIGNING PROTEIN BACKBONES")
    print("=" * 60)

    # Use Febuxostat as reference
    ref_u = universes['Febuxostat']
    mobile_u = universes['Natural']

    ref_ca = ref_u.select_atoms('protein and name CA')
    mobile_ca = mobile_u.select_atoms('protein and name CA')

    # Align using Kabsch
    aligned_ca, R, backbone_rmsd, ref_center, mobile_center = kabsch_align(
        mobile_ca.positions, ref_ca.positions
    )

    print(f"Backbone CA RMSD after alignment: {backbone_rmsd:.2f} A")

    # Apply transformation to entire Natural universe
    mobile_protein = mobile_u.select_atoms('protein or resname UNK')
    mobile_positions = mobile_protein.positions - mobile_center
    aligned_positions = (mobile_positions @ R) + ref_center

    # Get aligned ligand positions
    mobile_ligand = mobile_u.select_atoms('resname UNK')
    mobile_lig_idx = [i for i, atom in enumerate(mobile_protein) if atom.resname == 'UNK']
    aligned_ligand_pos = aligned_positions[mobile_lig_idx]

    # =========================================================================
    # CALCULATE LIGAND RMSD
    # =========================================================================
    print("\n" + "=" * 60)
    print("LIGAND COMPARISON")
    print("=" * 60)

    ref_ligand = ref_u.select_atoms('resname UNK')

    # Get heavy atoms only (non-hydrogen)
    ref_heavy = ref_u.select_atoms('resname UNK and not name H*')
    mobile_heavy = mobile_u.select_atoms('resname UNK and not name H*')

    print(f"\nFebuxostat heavy atoms: {len(ref_heavy)}")
    print(f"Natural heavy atoms: {len(mobile_heavy)}")

    # Get aligned heavy atom positions for Natural
    mobile_heavy_idx = [i for i, atom in enumerate(mobile_protein)
                        if atom.resname == 'UNK' and not atom.name.startswith('H')]
    aligned_heavy_pos = aligned_positions[mobile_heavy_idx]

    # Calculate center of mass for each ligand
    ref_lig_com = ref_ligand.positions.mean(axis=0)
    aligned_lig_com = aligned_ligand_pos.mean(axis=0)

    com_distance = np.linalg.norm(ref_lig_com - aligned_lig_com)
    print(f"\nLigand center-of-mass distance: {com_distance:.2f} A")

    # For RMSD, we need to compare similar atoms - use COM-based comparison
    # Since ligands are different, we compare binding site overlap instead

    # =========================================================================
    # CONTACT ANALYSIS
    # =========================================================================
    print("\n" + "=" * 60)
    print("PROTEIN-LIGAND CONTACT ANALYSIS")
    print("=" * 60)

    # Get contacts for each ligand
    feb_contacts = get_contacts(ref_u, cutoff=4.5)
    nat_contacts = get_contacts(mobile_u, cutoff=4.5)

    feb_residues = set((c[0], c[1]) for c in feb_contacts)
    nat_residues = set((c[0], c[1]) for c in nat_contacts)

    # Find overlapping and unique contacts
    shared_residues = feb_residues & nat_residues
    feb_unique = feb_residues - nat_residues
    nat_unique = nat_residues - feb_residues

    print(f"\nFebuxostat contacts: {len(feb_residues)} residues")
    print(f"Natural contacts: {len(nat_residues)} residues")
    print(f"Shared contacts: {len(shared_residues)} residues")
    print(f"Febuxostat unique: {len(feb_unique)} residues")
    print(f"Natural unique: {len(nat_unique)} residues")

    # Calculate overlap percentage
    overlap_pct = 100 * len(shared_residues) / max(len(feb_residues), len(nat_residues))
    print(f"\nBinding site overlap: {overlap_pct:.1f}%")

    # List shared contacts
    print("\n--- SHARED CONTACTS ---")
    shared_list = sorted(shared_residues)
    for i, (resid, resname) in enumerate(shared_list):
        if i < 20:  # Show first 20
            print(f"  {resname} {resid}")
    if len(shared_list) > 20:
        print(f"  ... and {len(shared_list) - 20} more")

    # List unique contacts
    print("\n--- FEBUXOSTAT UNIQUE CONTACTS ---")
    for resid, resname in sorted(feb_unique)[:10]:
        print(f"  {resname} {resid}")

    print("\n--- NATURAL UNIQUE CONTACTS ---")
    for resid, resname in sorted(nat_unique)[:10]:
        print(f"  {resname} {resid}")

    # =========================================================================
    # BINDING MODE CLASSIFICATION
    # =========================================================================
    print("\n" + "=" * 60)
    print("BINDING MODE CLASSIFICATION")
    print("=" * 60)

    if com_distance < 2.0 and overlap_pct > 70:
        classification = "SIMILAR BINDING MODE"
        interpretation = "Both ligands occupy the same binding pocket with similar orientation"
    elif com_distance < 4.0 and overlap_pct > 50:
        classification = "OVERLAPPING BINDING MODE"
        interpretation = "Ligands share significant binding site but with some differences"
    elif com_distance > 4.0 or overlap_pct < 30:
        classification = "DIFFERENT BINDING MODE"
        interpretation = "Ligands bind to different regions or with different orientations"
    else:
        classification = "PARTIALLY OVERLAPPING"
        interpretation = "Moderate overlap in binding site"

    print(f"\nClassification: {classification}")
    print(f"  COM distance: {com_distance:.2f} A")
    print(f"  Contact overlap: {overlap_pct:.1f}%")
    print(f"  Interpretation: {interpretation}")

    # =========================================================================
    # CREATE VISUALIZATION
    # =========================================================================
    print("\n" + "=" * 60)
    print("GENERATING VISUALIZATION")
    print("=" * 60)

    fig = plt.figure(figsize=(16, 6))

    # Plot 1: 3D overlay of ligand positions
    ax1 = fig.add_subplot(131, projection='3d')

    # Plot Febuxostat ligand
    feb_pos = ref_ligand.positions
    ax1.scatter(feb_pos[:, 0], feb_pos[:, 1], feb_pos[:, 2],
                c='blue', s=50, alpha=0.8, label='Febuxostat')

    # Plot Natural ligand (aligned)
    ax1.scatter(aligned_ligand_pos[:, 0], aligned_ligand_pos[:, 1], aligned_ligand_pos[:, 2],
                c='green', s=50, alpha=0.8, label='Natural')

    # Plot binding site residues
    binding_site = ref_u.select_atoms(f'protein and name CA and resid {" ".join(str(r[0]) for r in shared_residues)}')
    if len(binding_site) > 0:
        bs_pos = binding_site.positions
        ax1.scatter(bs_pos[:, 0], bs_pos[:, 1], bs_pos[:, 2],
                    c='gray', s=20, alpha=0.3, label='Binding site CA')

    ax1.set_xlabel('X (A)')
    ax1.set_ylabel('Y (A)')
    ax1.set_zlabel('Z (A)')
    ax1.set_title('Ligand Pose Overlay\n(aligned proteins)')
    ax1.legend()

    # Plot 2: Contact residue comparison (bar chart)
    ax2 = fig.add_subplot(132)

    categories = ['Shared', 'Febuxostat\nUnique', 'Natural\nUnique']
    counts = [len(shared_residues), len(feb_unique), len(nat_unique)]
    colors = ['purple', 'blue', 'green']

    bars = ax2.bar(categories, counts, color=colors, alpha=0.7, edgecolor='black')
    ax2.set_ylabel('Number of Residues')
    ax2.set_title('Binding Site Contact Analysis')

    for bar, count in zip(bars, counts):
        ax2.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.5,
                str(count), ha='center', va='bottom', fontweight='bold')

    # Plot 3: Summary stats
    ax3 = fig.add_subplot(133)
    ax3.axis('off')

    summary_text = f"""
    BINDING MODE COMPARISON
    ========================

    Backbone RMSD: {backbone_rmsd:.2f} A

    Ligand COM Distance: {com_distance:.2f} A

    Contact Analysis:
      Febuxostat:  {len(feb_residues)} residues
      Natural:     {len(nat_residues)} residues
      Shared:      {len(shared_residues)} residues
      Overlap:     {overlap_pct:.1f}%

    Classification:
    {classification}

    Threshold Criteria:
      Similar: COM < 2A, overlap > 70%
      Different: COM > 4A, overlap < 30%
    """

    ax3.text(0.1, 0.9, summary_text, transform=ax3.transAxes, fontsize=11,
             verticalalignment='top', fontfamily='monospace',
             bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

    plt.tight_layout()

    # Save plot
    plot_path = os.path.join(output_dir, 'pose_alignment.png')
    plt.savefig(plot_path, dpi=150, bbox_inches='tight')
    print(f"Saved: {plot_path}")

    pdf_path = os.path.join(output_dir, 'pose_alignment.pdf')
    plt.savefig(pdf_path, bbox_inches='tight')
    print(f"Saved: {pdf_path}")

    plt.close()

    # =========================================================================
    # SAVE CONTACT DATA
    # =========================================================================
    contact_data = []

    for resid, resname in feb_residues:
        contact_data.append({
            'resid': resid,
            'resname': resname,
            'febuxostat': 'YES',
            'natural': 'YES' if (resid, resname) in nat_residues else 'NO',
            'category': 'shared' if (resid, resname) in shared_residues else 'febuxostat_unique'
        })

    for resid, resname in nat_unique:
        contact_data.append({
            'resid': resid,
            'resname': resname,
            'febuxostat': 'NO',
            'natural': 'YES',
            'category': 'natural_unique'
        })

    df_contacts = pd.DataFrame(contact_data)
    df_contacts = df_contacts.sort_values('resid')

    csv_path = os.path.join(output_dir, 'binding_site_comparison.csv')
    df_contacts.to_csv(csv_path, index=False)
    print(f"Saved: {csv_path}")

    # =========================================================================
    # FINAL SUMMARY
    # =========================================================================
    print("\n" + "=" * 80)
    print("FINAL SUMMARY")
    print("=" * 80)

    print(f"""
    Ligand COM Distance:  {com_distance:.2f} A
    Contact Overlap:      {overlap_pct:.1f}%

    Classification: {classification}

    Interpretation: {interpretation}

    Mechanistic Implication:
    """)

    if overlap_pct > 50:
        print("    Both ligands likely act through SIMILAR MECHANISM")
        print("    They compete for the same binding site")
        print("    Similar pharmacological effects expected")
    else:
        print("    Ligands may act through DIFFERENT MECHANISMS")
        print("    They may have distinct pharmacological profiles")
        print("    Combination therapy could be synergistic")

    print("\n" + "=" * 80)
    print(f"Completed: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("=" * 80)

    return {
        'com_distance': com_distance,
        'overlap_pct': overlap_pct,
        'classification': classification,
        'shared_contacts': shared_residues,
        'feb_unique': feb_unique,
        'nat_unique': nat_unique
    }


if __name__ == '__main__':
    results = main()
