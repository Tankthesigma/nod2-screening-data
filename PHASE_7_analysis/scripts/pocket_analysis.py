#!/usr/bin/env python3
"""
CRITICAL ANALYSIS: Is the ligand in the pocket or floating away?

1. Extract frame 0 and frame 100 for visualization
2. Calculate pocket-based metrics (NOT just RMSD < 3Å)
3. Define pocket from frame 0 contacts
4. Report: pocket occupancy, min distance to pocket, contact overlap
"""

import MDAnalysis as mda
from MDAnalysis.analysis import align
from MDAnalysis.lib.distances import distance_array
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import warnings

warnings.filterwarnings('ignore')

TRAJ_DIR = Path(r'C:\Users\vasud\nod2-screening-data\PHASE_6\trajectories')
OUTPUT_DIR = Path(r'C:\Users\vasud\nod2-screening-data\PHASE_7_analysis\pocket_analysis')
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

COMPOUNDS = {
    'budesonide': {'reps': 3, 'top': 'budesonide_rep1_solvated.pdb'},
    'febuxostat': {'reps': 3, 'top': 'febuxostat_rep1_solvated.pdb'},
    'ursodiol': {'reps': 3, 'top': 'ursodiol_rep3_solvated.pdb'},
    'natural_top': {'reps': 3, 'top': 'natural_top_rep1_solvated.pdb'},
    'decoy': {'reps': 1, 'top': 'decoy_rep1_solvated.pdb'},
}

WATER_RESNAMES = {"HOH", "WAT", "SOL", "TIP3", "TIP3P"}
ION_RESNAMES = {"NA", "CL", "K", "CA", "MG", "NA+", "CL-"}


def pick_ligand(universe):
    """Find the ligand."""
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
    """Apply minimum image convention."""
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
    """
    Define pocket residues from frame 0 as any residue with heavy atoms
    within cutoff of any ligand heavy atom.
    """
    universe.trajectory[0]

    lig_heavy = universe.select_atoms(f"({lig_sel}) and not name H*")
    protein_heavy = universe.select_atoms("protein and not name H*")

    box = universe.dimensions
    dists = distance_array(lig_heavy.positions, protein_heavy.positions, box=box)

    # Find protein atoms within cutoff of any ligand atom
    contact_mask = np.any(dists < cutoff, axis=0)
    contact_atoms = protein_heavy[contact_mask]

    # Get unique residues
    pocket_resids = set(contact_atoms.resids)
    pocket_residues = [(r.resname, r.resid) for r in contact_atoms.residues]
    pocket_residues = list(set(pocket_residues))

    return pocket_resids, pocket_residues


def compute_pocket_metrics(universe, lig_sel, pocket_resids, stride=10):
    """
    Compute pocket-based metrics:
    1. Number of contacts per frame (residues within 5Å)
    2. Pocket occupancy (% frames with ≥5 contacts)
    3. Minimum ligand-pocket distance per frame
    """
    lig_heavy = universe.select_atoms(f"({lig_sel}) and not name H*")
    if len(lig_heavy) == 0:
        lig_heavy = universe.select_atoms(lig_sel)

    # Select pocket residue atoms
    pocket_sel = " or ".join([f"resid {r}" for r in pocket_resids])
    pocket_atoms = universe.select_atoms(f"protein and ({pocket_sel}) and not name H*")

    protein_bb = universe.select_atoms("protein and backbone")

    times = []
    n_contacts = []
    min_distances = []
    contact_residues_per_frame = []

    # Reference for alignment
    universe.trajectory[0]
    ref_prot_pos = protein_bb.positions.copy()
    ref_prot_com = ref_prot_pos.mean(axis=0)

    for ts in universe.trajectory[::stride]:
        box = ts.dimensions

        # PBC wrap ligand
        cur_prot_com = protein_bb.positions.mean(axis=0)
        lig_pos = wrap_ligand_to_protein(lig_heavy.positions.copy(), cur_prot_com, box)

        # Calculate distances to pocket
        pocket_pos = pocket_atoms.positions
        dists = distance_array(lig_pos, pocket_pos, box=box)

        # Min distance from any ligand atom to any pocket atom
        min_dist = dists.min()
        min_distances.append(min_dist)

        # Count residues in contact (any atom within 5Å)
        contact_mask = np.any(dists < 5.0, axis=0)
        contact_atoms_frame = pocket_atoms[contact_mask]
        contacting_resids = set(contact_atoms_frame.resids)
        n_contacts.append(len(contacting_resids))
        contact_residues_per_frame.append(contacting_resids)

        times.append(ts.time / 1000)

    # Calculate pocket occupancy (% frames with ≥5 contacts)
    pocket_occupancy = sum(1 for n in n_contacts if n >= 5) / len(n_contacts) * 100

    return {
        'times': np.array(times),
        'n_contacts': np.array(n_contacts),
        'min_distances': np.array(min_distances),
        'pocket_occupancy': pocket_occupancy,
        'contact_residues_per_frame': contact_residues_per_frame
    }


def extract_frames_for_visualization(name, rep, frames=[0, 100, 1000]):
    """Extract specific frames as PDB files for VMD/PyMOL visualization."""
    config = COMPOUNDS[name]
    u = mda.Universe(
        str(TRAJ_DIR / config['top']),
        str(TRAJ_DIR / f'{name}_rep{rep}.dcd')
    )

    lig_sel, _ = pick_ligand(u)
    protein = u.select_atoms("protein")
    ligand = u.select_atoms(lig_sel)

    # Select protein + ligand only (no water)
    complex_sel = protein + ligand

    for frame in frames:
        if frame >= len(u.trajectory):
            continue
        u.trajectory[frame]

        outfile = OUTPUT_DIR / f'{name}_rep{rep}_frame{frame}.pdb'
        complex_sel.write(str(outfile))
        print(f"  Saved: {outfile.name}")


def main():
    print("=" * 70)
    print("POCKET-BASED ANALYSIS: Is the ligand IN the pocket?")
    print("=" * 70)

    # =========================================================================
    # PART 1: Extract frames for visualization
    # =========================================================================
    print("\n### PART 1: Extracting Frames for Visualization ###\n")

    # Budesonide Rep1
    print("Extracting budesonide_rep1 frames 0, 100, 1000...")
    extract_frames_for_visualization('budesonide', 1, [0, 100, 1000])

    # Febuxostat Rep2 (the "hero" replicate)
    print("\nExtracting febuxostat_rep2 frames 0, 100, 1000...")
    extract_frames_for_visualization('febuxostat', 2, [0, 100, 1000])

    # Decoy for comparison
    print("\nExtracting decoy_rep1 frames 0, 100, 1000...")
    extract_frames_for_visualization('decoy', 1, [0, 100, 1000])

    # =========================================================================
    # PART 2: Define pocket from budesonide frame 0
    # =========================================================================
    print("\n### PART 2: Defining Binding Pocket from Budesonide Frame 0 ###\n")

    u_bud = mda.Universe(
        str(TRAJ_DIR / 'budesonide_rep1_solvated.pdb'),
        str(TRAJ_DIR / 'budesonide_rep1.dcd')
    )
    lig_sel_bud, _ = pick_ligand(u_bud)

    # Get pocket residues from budesonide frame 0
    reference_pocket_resids, reference_pocket_residues = get_pocket_residues(u_bud, lig_sel_bud, cutoff=5.0)

    print(f"Reference Pocket ({len(reference_pocket_residues)} residues from budesonide frame 0):")
    for resname, resid in sorted(reference_pocket_residues, key=lambda x: x[1]):
        print(f"  {resname}{resid}")

    # =========================================================================
    # PART 3: Compute pocket metrics for all compounds
    # =========================================================================
    print("\n### PART 3: Pocket-Based Metrics ###\n")

    all_results = {}

    print(f"{'Compound':<15} {'Rep':>4} {'Pocket Occ':>12} {'Mean Contacts':>14} {'Min Dist':>10}")
    print("-" * 65)

    for name, config in COMPOUNDS.items():
        all_results[name] = []

        for rep in range(1, config['reps'] + 1):
            traj_file = TRAJ_DIR / f"{name}_rep{rep}.dcd"
            if not traj_file.exists():
                continue

            u = mda.Universe(str(TRAJ_DIR / config['top']), str(traj_file))
            lig_sel, _ = pick_ligand(u)

            if lig_sel is None:
                continue

            # Use the REFERENCE pocket from budesonide
            metrics = compute_pocket_metrics(u, lig_sel, reference_pocket_resids, stride=10)

            all_results[name].append(metrics)

            print(f"{name:<15} {rep:>4} {metrics['pocket_occupancy']:>11.1f}% "
                  f"{metrics['n_contacts'].mean():>13.1f} "
                  f"{metrics['min_distances'].mean():>9.2f}Å")

    # =========================================================================
    # PART 4: Summary Table
    # =========================================================================
    print("\n" + "=" * 70)
    print("SUMMARY: POCKET-BASED METRICS (Better than RMSD < 3Å!)")
    print("=" * 70)
    print(f"\n{'Compound':<15} {'Avg Pocket Occ':>15} {'Avg Contacts':>13} {'Avg Min Dist':>13}")
    print("-" * 60)

    summary_data = {}
    for name, results in all_results.items():
        if results:
            avg_occ = np.mean([r['pocket_occupancy'] for r in results])
            avg_contacts = np.mean([r['n_contacts'].mean() for r in results])
            avg_min_dist = np.mean([r['min_distances'].mean() for r in results])
            summary_data[name] = {
                'pocket_occupancy': avg_occ,
                'avg_contacts': avg_contacts,
                'avg_min_dist': avg_min_dist
            }
            print(f"{name:<15} {avg_occ:>14.1f}% {avg_contacts:>12.1f} {avg_min_dist:>12.2f}Å")

    # =========================================================================
    # PART 5: Contact Overlap with Budesonide
    # =========================================================================
    print("\n### Contact Residue Overlap with Budesonide ###\n")

    # Get budesonide's most frequent contacts
    bud_all_contacts = set()
    for r in all_results['budesonide']:
        for frame_contacts in r['contact_residues_per_frame']:
            bud_all_contacts.update(frame_contacts)

    print(f"Budesonide contacts {len(bud_all_contacts)} unique residues across all frames")

    for name in ['febuxostat', 'ursodiol', 'natural_top', 'decoy']:
        if name in all_results and all_results[name]:
            other_contacts = set()
            for r in all_results[name]:
                for frame_contacts in r['contact_residues_per_frame']:
                    other_contacts.update(frame_contacts)

            overlap = bud_all_contacts.intersection(other_contacts)
            overlap_pct = len(overlap) / len(bud_all_contacts) * 100 if bud_all_contacts else 0
            print(f"{name}: {len(overlap)}/{len(bud_all_contacts)} residue overlap ({overlap_pct:.1f}%)")

    # =========================================================================
    # PART 6: Visualization Plots
    # =========================================================================
    print("\n### Generating Plots ###\n")

    # Plot: Contacts over time for budesonide and febuxostat
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    # Budesonide Rep1
    if all_results['budesonide']:
        ax = axes[0, 0]
        r = all_results['budesonide'][0]  # Rep1
        ax.plot(r['times'], r['n_contacts'], 'b-', linewidth=1)
        ax.axhline(y=5, color='green', linestyle='--', alpha=0.7, label='5 contacts threshold')
        ax.fill_between(r['times'], 0, r['n_contacts'], where=np.array(r['n_contacts']) >= 5,
                       color='green', alpha=0.3, label=f'In pocket ({r["pocket_occupancy"]:.1f}%)')
        ax.set_ylabel('# Pocket Contacts')
        ax.set_title(f'Budesonide Rep1 - Pocket Occupancy: {r["pocket_occupancy"]:.1f}%')
        ax.legend()
        ax.set_ylim(0, max(15, max(r['n_contacts']) + 2))

    # Febuxostat Rep2 (hero)
    if len(all_results['febuxostat']) >= 2:
        ax = axes[0, 1]
        r = all_results['febuxostat'][1]  # Rep2
        ax.plot(r['times'], r['n_contacts'], 'g-', linewidth=1)
        ax.axhline(y=5, color='green', linestyle='--', alpha=0.7)
        ax.fill_between(r['times'], 0, r['n_contacts'], where=np.array(r['n_contacts']) >= 5,
                       color='green', alpha=0.3, label=f'In pocket ({r["pocket_occupancy"]:.1f}%)')
        ax.set_ylabel('# Pocket Contacts')
        ax.set_title(f'Febuxostat Rep2 - Pocket Occupancy: {r["pocket_occupancy"]:.1f}%')
        ax.legend()
        ax.set_ylim(0, max(15, max(r['n_contacts']) + 2))

    # Decoy
    if all_results['decoy']:
        ax = axes[1, 0]
        r = all_results['decoy'][0]
        ax.plot(r['times'], r['n_contacts'], 'r-', linewidth=1)
        ax.axhline(y=5, color='green', linestyle='--', alpha=0.7)
        ax.fill_between(r['times'], 0, r['n_contacts'], where=np.array(r['n_contacts']) >= 5,
                       color='green', alpha=0.3, label=f'In pocket ({r["pocket_occupancy"]:.1f}%)')
        ax.set_xlabel('Time (ns)')
        ax.set_ylabel('# Pocket Contacts')
        ax.set_title(f'Decoy Rep1 - Pocket Occupancy: {r["pocket_occupancy"]:.1f}%')
        ax.legend()
        ax.set_ylim(0, max(15, max(r['n_contacts']) + 2))

    # Min distance plot
    ax = axes[1, 1]
    for name, color in [('budesonide', 'blue'), ('febuxostat', 'green'), ('decoy', 'red')]:
        if all_results[name]:
            if name == 'febuxostat' and len(all_results[name]) >= 2:
                r = all_results[name][1]  # Rep2
                label = f'{name} (Rep2)'
            else:
                r = all_results[name][0]
                label = f'{name} (Rep1)'
            ax.plot(r['times'], r['min_distances'], color=color, linewidth=1, label=label, alpha=0.8)

    ax.axhline(y=5, color='gray', linestyle='--', alpha=0.5, label='5Å threshold')
    ax.set_xlabel('Time (ns)')
    ax.set_ylabel('Min Distance to Pocket (Å)')
    ax.set_title('Ligand-Pocket Distance Over Time')
    ax.legend()
    ax.set_ylim(0, 20)

    plt.tight_layout()
    plt.savefig(OUTPUT_DIR / 'pocket_occupancy_analysis.png', dpi=300)
    plt.close()
    print(f"Saved: {OUTPUT_DIR / 'pocket_occupancy_analysis.png'}")

    # =========================================================================
    # PART 7: Frame 0 vs Frame 100 comparison
    # =========================================================================
    print("\n### Frame 0 vs Frame 100 Analysis ###\n")

    for name, rep in [('budesonide', 1), ('febuxostat', 2)]:
        config = COMPOUNDS[name]
        u = mda.Universe(
            str(TRAJ_DIR / config['top']),
            str(TRAJ_DIR / f'{name}_rep{rep}.dcd')
        )
        lig_sel, _ = pick_ligand(u)
        lig_heavy = u.select_atoms(f"({lig_sel}) and not name H*")
        protein_bb = u.select_atoms("protein and backbone")

        # Frame 0
        u.trajectory[0]
        box = u.dimensions
        prot_com_0 = protein_bb.positions.mean(axis=0)
        lig_pos_0 = wrap_ligand_to_protein(lig_heavy.positions.copy(), prot_com_0, box)
        lig_com_0 = lig_pos_0.mean(axis=0)

        # Frame 100
        u.trajectory[100]
        box = u.dimensions
        prot_com_100 = protein_bb.positions.mean(axis=0)
        lig_pos_100 = wrap_ligand_to_protein(lig_heavy.positions.copy(), prot_com_100, box)
        lig_com_100 = lig_pos_100.mean(axis=0)

        # COM displacement
        com_displacement = np.linalg.norm(lig_com_100 - lig_com_0)

        print(f"{name} Rep{rep}:")
        print(f"  Frame 0  Ligand COM: [{lig_com_0[0]:.1f}, {lig_com_0[1]:.1f}, {lig_com_0[2]:.1f}]")
        print(f"  Frame 100 Ligand COM: [{lig_com_100[0]:.1f}, {lig_com_100[1]:.1f}, {lig_com_100[2]:.1f}]")
        print(f"  COM Displacement: {com_displacement:.2f} Å")
        print()

    print("\n" + "=" * 70)
    print("CONCLUSION: Check the extracted PDB files in VMD/PyMOL!")
    print(f"Files saved in: {OUTPUT_DIR}")
    print("=" * 70)


if __name__ == "__main__":
    main()
