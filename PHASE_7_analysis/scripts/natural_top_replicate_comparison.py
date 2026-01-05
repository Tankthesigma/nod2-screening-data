#!/usr/bin/env python3
"""
NATURAL_TOP REPLICATE COMPARISON
Rep1 (99.5% occ) vs Rep2 (46.5% occ) vs Rep3 (93% occ)

Analysis:
1. Compare starting poses (frame 0) - are they identical?
2. Compare H-bonds at frame 0
3. Check pocket contacts and position (center vs edge)
"""

import MDAnalysis as mda
from MDAnalysis.lib.distances import distance_array
from MDAnalysis.analysis import align
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import warnings

warnings.filterwarnings('ignore')

TRAJ_DIR = Path(r'C:\Users\vasud\nod2-screening-data\PHASE_6\trajectories')
OUTPUT_DIR = Path(r'C:\Users\vasud\nod2-screening-data\PHASE_7_analysis\natural_top_comparison')
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

WATER_RESNAMES = {"HOH", "WAT", "SOL", "TIP3", "TIP3P"}
ION_RESNAMES = {"NA", "CL", "K", "CA", "MG", "NA+", "CL-"}

# Pocket occupancy from full_pocket_metrics.csv
POCKET_OCC = {1: 99.5, 2: 46.5, 3: 93.0}


def pick_ligand(universe):
    """Find the ligand (largest non-protein, non-water residue)."""
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


def get_frame0_data(universe, lig_sel):
    """Get ligand positions and protein data at frame 0."""
    lig_heavy = universe.select_atoms(f"({lig_sel}) and not name H*")
    if len(lig_heavy) == 0:
        lig_heavy = universe.select_atoms(lig_sel)

    protein_bb = universe.select_atoms("protein and backbone")

    universe.trajectory[0]
    box = universe.dimensions[:3]
    prot_com = protein_bb.positions.mean(axis=0)
    lig_pos = wrap_ligand_to_protein(lig_heavy.positions.copy(), prot_com, box)
    lig_com = lig_pos.mean(axis=0)

    return {
        'lig_pos': lig_pos,
        'lig_com': lig_com,
        'prot_bb_pos': protein_bb.positions.copy(),
        'prot_com': prot_com,
        'box': box,
        'n_lig_atoms': len(lig_heavy)
    }


def count_hbonds_frame0(universe, lig_sel, cutoff=3.5):
    """Count H-bonds at frame 0 using distance-based detection."""
    universe.trajectory[0]

    protein_polar = universe.select_atoms("protein and (name N* O*)")
    ligand_polar = universe.select_atoms(f"({lig_sel}) and (name N* O*)")

    if len(ligand_polar) == 0:
        return 0, []

    box = universe.dimensions
    dists = distance_array(protein_polar.positions, ligand_polar.positions, box=box)

    contacts = np.where(dists < cutoff)

    # Get residue details with distances
    hbond_details = []
    for idx, (prot_idx, lig_idx) in enumerate(zip(contacts[0], contacts[1])):
        atom = protein_polar[prot_idx]
        dist = dists[prot_idx, lig_idx]
        hbond_details.append({
            'residue': f"{atom.resname}{atom.resid}",
            'atom': atom.name,
            'distance': dist
        })

    return len(contacts[0]), hbond_details


def get_pocket_contacts_frame0(universe, lig_sel, pocket_resids, cutoff=5.0):
    """Get detailed pocket contacts at frame 0."""
    universe.trajectory[0]

    lig_heavy = universe.select_atoms(f"({lig_sel}) and not name H*")
    if len(lig_heavy) == 0:
        lig_heavy = universe.select_atoms(lig_sel)

    pocket_sel = " or ".join([f"resid {r}" for r in pocket_resids])
    pocket_atoms = universe.select_atoms(f"protein and ({pocket_sel}) and not name H*")

    box = universe.dimensions
    prot_bb = universe.select_atoms("protein and backbone")
    prot_com = prot_bb.positions.mean(axis=0)
    lig_pos = wrap_ligand_to_protein(lig_heavy.positions.copy(), prot_com, box[:3])

    dists = distance_array(lig_pos, pocket_atoms.positions, box=box)

    # Get minimum distance to each pocket residue
    residue_distances = {}
    for i, atom in enumerate(pocket_atoms):
        resid = atom.resid
        min_dist = dists[:, i].min()
        if resid not in residue_distances or min_dist < residue_distances[resid]:
            residue_distances[resid] = min_dist

    contact_resids = [r for r, d in residue_distances.items() if d < cutoff]

    return len(contact_resids), residue_distances


def compute_pocket_center(universe, pocket_resids):
    """Compute the geometric center of the binding pocket."""
    universe.trajectory[0]
    pocket_sel = " or ".join([f"resid {r}" for r in pocket_resids])
    pocket_atoms = universe.select_atoms(f"protein and ({pocket_sel}) and name CA")
    return pocket_atoms.positions.mean(axis=0)


def main():
    print("=" * 80)
    print("NATURAL_TOP REPLICATE COMPARISON")
    print("Rep1 (99.5%) vs Rep2 (46.5%) vs Rep3 (93%)")
    print("=" * 80)

    # Load all three replicates
    reps = {}
    for rep in [1, 2, 3]:
        u = mda.Universe(
            str(TRAJ_DIR / 'natural_top_rep1_solvated.pdb'),
            str(TRAJ_DIR / f'natural_top_rep{rep}.dcd')
        )
        lig_sel, _ = pick_ligand(u)
        reps[rep] = {'universe': u, 'lig_sel': lig_sel}

    # Define reference pocket from budesonide (same as other analyses)
    u_bud = mda.Universe(
        str(TRAJ_DIR / 'budesonide_rep1_solvated.pdb'),
        str(TRAJ_DIR / 'budesonide_rep1.dcd')
    )
    lig_sel_bud, _ = pick_ligand(u_bud)
    u_bud.trajectory[0]

    lig_heavy_bud = u_bud.select_atoms(f"({lig_sel_bud}) and not name H*")
    protein_heavy = u_bud.select_atoms("protein and not name H*")
    box = u_bud.dimensions
    dists = distance_array(lig_heavy_bud.positions, protein_heavy.positions, box=box)
    contact_mask = np.any(dists < 5.0, axis=0)
    contact_atoms = protein_heavy[contact_mask]
    pocket_resids = set(contact_atoms.resids)

    print(f"\nReference pocket (from budesonide frame 0): {len(pocket_resids)} residues")
    print(f"Pocket residues: {sorted(pocket_resids)}")

    # Compute pocket center using natural_top universe
    pocket_center = compute_pocket_center(reps[1]['universe'], pocket_resids)
    print(f"Pocket center (CA atoms): ({pocket_center[0]:.2f}, {pocket_center[1]:.2f}, {pocket_center[2]:.2f})")

    # =========================================================================
    # PART 1: Compare starting poses
    # =========================================================================
    print("\n" + "=" * 80)
    print("PART 1: STARTING POSE COMPARISON (Frame 0)")
    print("=" * 80)

    frame0_data = {}
    for rep in [1, 2, 3]:
        u = reps[rep]['universe']
        lig_sel = reps[rep]['lig_sel']
        data = get_frame0_data(u, lig_sel)
        frame0_data[rep] = data

        print(f"\nRep{rep} Frame 0 (Pocket Occ: {POCKET_OCC[rep]}%):")
        print(f"  Ligand COM: ({data['lig_com'][0]:.2f}, {data['lig_com'][1]:.2f}, {data['lig_com'][2]:.2f})")

        # Distance from pocket center
        dist_to_pocket = np.linalg.norm(data['lig_com'] - pocket_center)
        print(f"  Distance from pocket center: {dist_to_pocket:.2f} A")
        frame0_data[rep]['dist_to_pocket_center'] = dist_to_pocket

    # =========================================================================
    # PART 2: RMSD between starting poses
    # =========================================================================
    print("\n" + "=" * 80)
    print("PART 2: RMSD BETWEEN STARTING POSES")
    print("=" * 80)

    print("\n(After backbone alignment)")

    pose_rmsd = {}
    for i in [1, 2, 3]:
        for j in [1, 2, 3]:
            if i >= j:
                continue

            # Use rep i as reference
            ref_prot = frame0_data[i]['prot_bb_pos']
            ref_lig = frame0_data[i]['lig_pos']

            mob_prot = frame0_data[j]['prot_bb_pos']
            mob_lig = frame0_data[j]['lig_pos']

            # Align proteins
            ref_com = ref_prot.mean(axis=0)
            mob_com = mob_prot.mean(axis=0)

            ref_centered = ref_prot - ref_com
            mob_centered = mob_prot - mob_com

            R, _ = align.rotation_matrix(mob_centered, ref_centered)

            # Transform mobile ligand
            lig_centered = mob_lig - mob_com
            lig_rotated = np.dot(lig_centered, R.T)
            lig_aligned = lig_rotated + ref_com

            # RMSD
            diff = lig_aligned - ref_lig
            rmsd = np.sqrt((diff ** 2).sum() / len(ref_lig))

            pose_rmsd[(i, j)] = rmsd
            print(f"\nRep{i} ({POCKET_OCC[i]}%) vs Rep{j} ({POCKET_OCC[j]}%) ligand RMSD: {rmsd:.3f} A")

            if rmsd < 0.5:
                print("  --> IDENTICAL starting poses")
            elif rmsd < 2.0:
                print("  --> SIMILAR starting poses (minor differences)")
            elif rmsd < 5.0:
                print("  --> MODERATELY DIFFERENT poses")
            else:
                print("  --> SIGNIFICANTLY DIFFERENT starting poses!")

    # =========================================================================
    # PART 3: H-bonds at frame 0
    # =========================================================================
    print("\n" + "=" * 80)
    print("PART 3: H-BONDS AT FRAME 0")
    print("=" * 80)

    hbond_data = {}
    for rep in [1, 2, 3]:
        u = reps[rep]['universe']
        lig_sel = reps[rep]['lig_sel']

        n_hbonds, hbond_details = count_hbonds_frame0(u, lig_sel, cutoff=3.5)
        hbond_data[rep] = {'count': n_hbonds, 'details': hbond_details}

        print(f"\nRep{rep} ({POCKET_OCC[rep]}% occ) Frame 0 H-bonds: {n_hbonds}")
        if hbond_details:
            # Group by residue
            residue_counts = {}
            for hb in hbond_details:
                res = hb['residue']
                if res not in residue_counts:
                    residue_counts[res] = []
                residue_counts[res].append(hb['distance'])

            for res, dists in sorted(residue_counts.items()):
                print(f"  {res}: {len(dists)} contact(s), min dist = {min(dists):.2f} A")
        else:
            print("  (No H-bonds detected)")

    # =========================================================================
    # PART 4: Pocket contacts at frame 0 (center vs edge)
    # =========================================================================
    print("\n" + "=" * 80)
    print("PART 4: POCKET CONTACTS AT FRAME 0 (Center vs Edge Analysis)")
    print("=" * 80)

    contact_data = {}
    for rep in [1, 2, 3]:
        u = reps[rep]['universe']
        lig_sel = reps[rep]['lig_sel']

        n_contacts, residue_distances = get_pocket_contacts_frame0(u, lig_sel, pocket_resids)
        contact_data[rep] = {'count': n_contacts, 'distances': residue_distances}

        print(f"\nRep{rep} ({POCKET_OCC[rep]}% occ) Frame 0 Pocket Contacts: {n_contacts}/{len(pocket_resids)}")
        print(f"  Distance to pocket center: {frame0_data[rep]['dist_to_pocket_center']:.2f} A")

        # Show distances to each pocket residue
        print(f"  Per-residue distances:")
        for resid in sorted(pocket_resids):
            dist = residue_distances.get(resid, float('inf'))
            status = "CONTACT" if dist < 5.0 else "no contact"
            print(f"    Res {resid}: {dist:.2f} A ({status})")

    # =========================================================================
    # PART 5: Visualization
    # =========================================================================
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    # Plot 1: Distance to pocket center
    ax1 = axes[0, 0]
    reps_list = [1, 2, 3]
    pocket_dists = [frame0_data[r]['dist_to_pocket_center'] for r in reps_list]
    colors = ['#2ca02c', '#d62728', '#1f77b4']  # Green for best, red for worst, blue for good

    bars = ax1.bar([f'Rep1\n({POCKET_OCC[1]}% occ)', f'Rep2\n({POCKET_OCC[2]}% occ)', f'Rep3\n({POCKET_OCC[3]}% occ)'],
                   pocket_dists, color=colors, edgecolor='black')
    ax1.set_ylabel('Distance from Pocket Center (A)', fontsize=11)
    ax1.set_title('Starting Position: Distance from Pocket Center', fontsize=12)
    for bar, val in zip(bars, pocket_dists):
        ax1.annotate(f'{val:.1f}A',
                     xy=(bar.get_x() + bar.get_width() / 2, bar.get_height()),
                     xytext=(0, 3), textcoords="offset points",
                     ha='center', va='bottom', fontweight='bold')
    ax1.grid(True, alpha=0.3, axis='y')

    # Plot 2: H-bonds at frame 0
    ax2 = axes[0, 1]
    hbond_counts = [hbond_data[r]['count'] for r in reps_list]
    bars = ax2.bar([f'Rep1\n({POCKET_OCC[1]}% occ)', f'Rep2\n({POCKET_OCC[2]}% occ)', f'Rep3\n({POCKET_OCC[3]}% occ)'],
                   hbond_counts, color=colors, edgecolor='black')
    ax2.set_ylabel('Number of H-bonds', fontsize=11)
    ax2.set_title('H-bonds at Frame 0', fontsize=12)
    for bar, val in zip(bars, hbond_counts):
        ax2.annotate(f'{val}',
                     xy=(bar.get_x() + bar.get_width() / 2, bar.get_height()),
                     xytext=(0, 3), textcoords="offset points",
                     ha='center', va='bottom', fontweight='bold', fontsize=14)
    ax2.grid(True, alpha=0.3, axis='y')

    # Plot 3: Pocket contacts at frame 0
    ax3 = axes[1, 0]
    contact_counts = [contact_data[r]['count'] for r in reps_list]
    bars = ax3.bar([f'Rep1\n({POCKET_OCC[1]}% occ)', f'Rep2\n({POCKET_OCC[2]}% occ)', f'Rep3\n({POCKET_OCC[3]}% occ)'],
                   contact_counts, color=colors, edgecolor='black')
    ax3.axhline(y=len(pocket_resids), color='green', linestyle='--', alpha=0.7,
                label=f'Max possible ({len(pocket_resids)})')
    ax3.set_ylabel('Number of Pocket Contacts', fontsize=11)
    ax3.set_title('Pocket Contacts at Frame 0', fontsize=12)
    ax3.legend()
    for bar, val in zip(bars, contact_counts):
        ax3.annotate(f'{val}',
                     xy=(bar.get_x() + bar.get_width() / 2, bar.get_height()),
                     xytext=(0, 3), textcoords="offset points",
                     ha='center', va='bottom', fontweight='bold', fontsize=14)
    ax3.grid(True, alpha=0.3, axis='y')

    # Plot 4: Pose RMSD comparison
    ax4 = axes[1, 1]
    rmsd_labels = ['Rep1 vs Rep2', 'Rep1 vs Rep3', 'Rep2 vs Rep3']
    rmsd_values = [pose_rmsd[(1, 2)], pose_rmsd[(1, 3)], pose_rmsd[(2, 3)]]
    bars = ax4.bar(rmsd_labels, rmsd_values, color=['#1f77b4', '#9467bd', '#8c564b'],
                   edgecolor='black')
    ax4.axhline(y=2.0, color='orange', linestyle='--', alpha=0.7, label='2A threshold')
    ax4.set_ylabel('Ligand RMSD (A)', fontsize=11)
    ax4.set_title('Starting Pose Differences (after backbone alignment)', fontsize=12)
    ax4.legend()
    for bar, val in zip(bars, rmsd_values):
        ax4.annotate(f'{val:.2f}A',
                     xy=(bar.get_x() + bar.get_width() / 2, bar.get_height()),
                     xytext=(0, 3), textcoords="offset points",
                     ha='center', va='bottom', fontweight='bold')
    ax4.grid(True, alpha=0.3, axis='y')

    plt.suptitle('NATURAL_TOP: Why did Rep1 (99.5%) and Rep3 (93%) outperform Rep2 (46.5%)?',
                 fontsize=14, fontweight='bold')
    plt.tight_layout()
    plt.savefig(OUTPUT_DIR / 'natural_top_replicate_comparison.png', dpi=300)
    plt.close()
    print(f"\nSaved: {OUTPUT_DIR / 'natural_top_replicate_comparison.png'}")

    # =========================================================================
    # PART 6: Extract frame 0 PDBs for visualization
    # =========================================================================
    print("\n" + "=" * 80)
    print("PART 6: EXTRACTING FRAME 0 PDBs FOR VISUALIZATION")
    print("=" * 80)

    for rep in [1, 2, 3]:
        u = reps[rep]['universe']
        u.trajectory[0]

        # Select protein + ligand
        lig_sel = reps[rep]['lig_sel']
        system = u.select_atoms(f"protein or ({lig_sel})")

        pdb_file = OUTPUT_DIR / f'natural_top_rep{rep}_frame0.pdb'
        system.write(str(pdb_file))
        print(f"Saved: {pdb_file}")

    # =========================================================================
    # FINAL VERDICT
    # =========================================================================
    print("\n" + "=" * 80)
    print("*** FINAL VERDICT: WHY DID REP2 UNDERPERFORM? ***")
    print("=" * 80)

    print("\nSUMMARY TABLE:")
    print("-" * 70)
    print(f"{'Replicate':<12} {'Pocket Occ':<12} {'H-bonds':<10} {'Contacts':<10} {'Dist to Center':<15}")
    print("-" * 70)
    for rep in [1, 2, 3]:
        print(f"Rep{rep:<10} {POCKET_OCC[rep]:<11}% {hbond_data[rep]['count']:<10} {contact_data[rep]['count']}/8{'':<6} {frame0_data[rep]['dist_to_pocket_center']:.2f} A")
    print("-" * 70)

    print("\nFINDINGS:")
    print("-" * 60)

    # Check if poses are identical
    max_pose_rmsd = max(pose_rmsd.values())
    if max_pose_rmsd < 0.5:
        print("1. STARTING POSES: All replicates started with IDENTICAL poses")
        print("   --> Differences arose from random velocity assignment")
    elif max_pose_rmsd < 2.0:
        print(f"1. STARTING POSES: Similar poses (max RMSD = {max_pose_rmsd:.2f}A)")
    else:
        print(f"1. STARTING POSES: Different poses detected (max RMSD = {max_pose_rmsd:.2f}A)")

    # H-bond comparison
    print(f"\n2. H-BONDS AT FRAME 0:")
    for rep in [1, 2, 3]:
        if POCKET_OCC[rep] == max(POCKET_OCC.values()):
            status = "BEST"
        elif POCKET_OCC[rep] == min(POCKET_OCC.values()):
            status = "WORST"
        else:
            status = "GOOD"
        print(f"   Rep{rep} ({status}): {hbond_data[rep]['count']} H-bonds")

    # Pocket contacts
    print(f"\n3. POCKET CONTACTS AT FRAME 0:")
    for rep in [1, 2, 3]:
        print(f"   Rep{rep}: {contact_data[rep]['count']}/{len(pocket_resids)} residues")

    # Position analysis
    print(f"\n4. DISTANCE FROM POCKET CENTER:")
    for rep in [1, 2, 3]:
        dist = frame0_data[rep]['dist_to_pocket_center']
        print(f"   Rep{rep}: {dist:.2f} A")

    print("\n" + "=" * 80)
    print("CONCLUSION:")
    print("=" * 80)

    # Determine main factors
    best_rep = max(POCKET_OCC, key=POCKET_OCC.get)
    worst_rep = min(POCKET_OCC, key=POCKET_OCC.get)

    conclusions = []
    if hbond_data[best_rep]['count'] > hbond_data[worst_rep]['count']:
        conclusions.append(f"  Rep{best_rep} had MORE H-BONDS ({hbond_data[best_rep]['count']} vs {hbond_data[worst_rep]['count']})")
    elif hbond_data[best_rep]['count'] < hbond_data[worst_rep]['count']:
        conclusions.append(f"  Rep{worst_rep} had more H-bonds but still failed - H-bonds not sufficient")

    if contact_data[best_rep]['count'] > contact_data[worst_rep]['count']:
        conclusions.append(f"  Rep{best_rep} had BETTER pocket coverage ({contact_data[best_rep]['count']}/8 vs {contact_data[worst_rep]['count']}/8)")

    if frame0_data[best_rep]['dist_to_pocket_center'] < frame0_data[worst_rep]['dist_to_pocket_center']:
        conclusions.append(f"  Rep{best_rep} started CLOSER to pocket center")

    if conclusions:
        for c in conclusions:
            print(c)
    else:
        print("  Differences may be due to subtle orientation or velocity effects")

    print("""
For your ISEF report:
  - Natural_top shows excellent binding in 2/3 replicates (99.5%, 93%)
  - Rep2 (46.5%) still shows moderate binding
  - This natural product demonstrates promising binding characteristics
""")
    print("=" * 80)


if __name__ == "__main__":
    main()
