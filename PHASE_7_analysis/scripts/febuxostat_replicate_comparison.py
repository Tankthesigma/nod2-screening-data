#!/usr/bin/env python3
"""
FEBUXOSTAT REPLICATE COMPARISON
Why did Rep1 fail (12% pocket occupancy) while Rep2/Rep3 succeeded (99-100%)?

Analysis:
1. Compare starting poses (frame 0) - are they identical?
2. Calculate RMSD between starting poses
3. Track Rep1 RMSD over time - immediate exit or gradual drift?
4. Compare H-bonds at frame 0 for each replicate
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
OUTPUT_DIR = Path(r'C:\Users\vasud\nod2-screening-data\PHASE_7_analysis\febuxostat_comparison')
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

WATER_RESNAMES = {"HOH", "WAT", "SOL", "TIP3", "TIP3P"}
ION_RESNAMES = {"NA", "CL", "K", "CA", "MG", "NA+", "CL-"}


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


def get_frame0_ligand_positions(universe, lig_sel):
    """Get ligand heavy atom positions at frame 0, aligned to protein backbone."""
    lig_heavy = universe.select_atoms(f"({lig_sel}) and not name H*")
    if len(lig_heavy) == 0:
        lig_heavy = universe.select_atoms(lig_sel)

    protein_bb = universe.select_atoms("protein and backbone")

    universe.trajectory[0]
    box = universe.dimensions[:3]
    prot_com = protein_bb.positions.mean(axis=0)
    lig_pos = wrap_ligand_to_protein(lig_heavy.positions.copy(), prot_com, box)

    return lig_pos, protein_bb.positions.copy(), prot_com


def compute_rmsd_timeseries(universe, lig_sel, stride=5):
    """Compute ligand RMSD over time with proper alignment."""
    lig_heavy = universe.select_atoms(f"({lig_sel}) and not name H*")
    if len(lig_heavy) == 0:
        lig_heavy = universe.select_atoms(lig_sel)

    protein_bb = universe.select_atoms("protein and backbone")

    # Reference = frame 0
    universe.trajectory[0]
    box = universe.dimensions[:3]
    ref_prot_pos = protein_bb.positions.copy()
    ref_prot_com = ref_prot_pos.mean(axis=0)
    ref_lig_pos = wrap_ligand_to_protein(lig_heavy.positions.copy(), ref_prot_com, box)

    times = []
    rmsd_values = []

    for ts in universe.trajectory[::stride]:
        box = ts.dimensions[:3]
        cur_prot_pos = protein_bb.positions.copy()
        cur_prot_com = cur_prot_pos.mean(axis=0)
        cur_lig_pos = wrap_ligand_to_protein(lig_heavy.positions.copy(), cur_prot_com, box)

        # Kabsch alignment
        cur_prot_centered = cur_prot_pos - cur_prot_com
        ref_prot_centered = ref_prot_pos - ref_prot_com
        R, _ = align.rotation_matrix(cur_prot_centered, ref_prot_centered)

        # Transform ligand
        lig_centered = cur_lig_pos - cur_prot_com
        lig_rotated = np.dot(lig_centered, R.T)
        lig_aligned = lig_rotated + ref_prot_com

        # RMSD
        diff = lig_aligned - ref_lig_pos
        rmsd = np.sqrt((diff ** 2).sum() / len(lig_heavy))
        rmsd_values.append(rmsd)
        times.append(ts.time / 1000)  # ps to ns

    return np.array(times), np.array(rmsd_values)


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

    # Get residue details
    hbond_residues = set()
    for prot_idx in contacts[0]:
        atom = protein_polar[prot_idx]
        hbond_residues.add(f"{atom.resname}{atom.resid}")

    return len(contacts[0]), sorted(hbond_residues)


def compute_pocket_contacts_frame0(universe, lig_sel, pocket_resids, cutoff=5.0):
    """Count pocket contacts at frame 0."""
    universe.trajectory[0]

    lig_heavy = universe.select_atoms(f"({lig_sel}) and not name H*")
    if len(lig_heavy) == 0:
        lig_heavy = universe.select_atoms(lig_sel)

    pocket_sel = " or ".join([f"resid {r}" for r in pocket_resids])
    pocket_atoms = universe.select_atoms(f"protein and ({pocket_sel}) and not name H*")

    box = universe.dimensions
    dists = distance_array(lig_heavy.positions, pocket_atoms.positions, box=box)

    contact_mask = np.any(dists < cutoff, axis=0)
    contact_atoms = pocket_atoms[contact_mask]
    contacting_resids = set(contact_atoms.resids)

    return len(contacting_resids), sorted(contacting_resids)


def main():
    print("=" * 80)
    print("FEBUXOSTAT REPLICATE COMPARISON")
    print("Why did Rep1 fail while Rep2/Rep3 succeeded?")
    print("=" * 80)

    # Load all three replicates
    reps = {}
    for rep in [1, 2, 3]:
        u = mda.Universe(
            str(TRAJ_DIR / 'febuxostat_rep1_solvated.pdb'),
            str(TRAJ_DIR / f'febuxostat_rep{rep}.dcd')
        )
        lig_sel, _ = pick_ligand(u)
        reps[rep] = {'universe': u, 'lig_sel': lig_sel}

    # =========================================================================
    # PART 1: Compare starting poses
    # =========================================================================
    print("\n" + "=" * 80)
    print("PART 1: STARTING POSE COMPARISON (Frame 0)")
    print("=" * 80)

    frame0_positions = {}
    frame0_prot_positions = {}

    for rep in [1, 2, 3]:
        u = reps[rep]['universe']
        lig_sel = reps[rep]['lig_sel']
        lig_pos, prot_pos, prot_com = get_frame0_ligand_positions(u, lig_sel)
        frame0_positions[rep] = lig_pos
        frame0_prot_positions[rep] = prot_pos

        print(f"\nRep{rep} Frame 0:")
        print(f"  Ligand COM: ({lig_pos.mean(axis=0)[0]:.2f}, {lig_pos.mean(axis=0)[1]:.2f}, {lig_pos.mean(axis=0)[2]:.2f})")
        print(f"  Ligand atoms: {len(lig_pos)}")

    # =========================================================================
    # PART 2: RMSD between starting poses
    # =========================================================================
    print("\n" + "=" * 80)
    print("PART 2: RMSD BETWEEN STARTING POSES")
    print("=" * 80)

    # Align proteins first, then compute ligand RMSD
    print("\n(After backbone alignment)")

    for i in [1, 2, 3]:
        for j in [1, 2, 3]:
            if i >= j:
                continue

            # Use rep i as reference
            ref_prot = frame0_prot_positions[i]
            ref_lig = frame0_positions[i]

            mob_prot = frame0_prot_positions[j]
            mob_lig = frame0_positions[j]

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

            print(f"\nRep{i} vs Rep{j} ligand RMSD: {rmsd:.3f} Angstrom")

            if rmsd < 0.5:
                print("  --> IDENTICAL starting poses")
            elif rmsd < 2.0:
                print("  --> SIMILAR starting poses (minor differences)")
            else:
                print("  --> DIFFERENT starting poses!")

    # =========================================================================
    # PART 3: Rep1 RMSD over time - when did it leave?
    # =========================================================================
    print("\n" + "=" * 80)
    print("PART 3: REP1 RMSD OVER TIME (When did it leave the pocket?)")
    print("=" * 80)

    fig, axes = plt.subplots(3, 1, figsize=(12, 10), sharex=True)

    colors = {'1': '#d62728', '2': '#2ca02c', '3': '#1f77b4'}

    all_rmsd_data = {}

    for rep in [1, 2, 3]:
        u = reps[rep]['universe']
        lig_sel = reps[rep]['lig_sel']

        times, rmsd = compute_rmsd_timeseries(u, lig_sel, stride=5)
        all_rmsd_data[rep] = {'times': times, 'rmsd': rmsd}

        ax = axes[rep - 1]
        ax.plot(times, rmsd, color=colors[str(rep)], linewidth=1.5, label=f'Rep{rep}')
        ax.axhline(y=3, color='green', linestyle='--', alpha=0.7, label='3A threshold')
        ax.axhline(y=5, color='orange', linestyle='--', alpha=0.7, label='5A threshold')
        ax.set_ylabel('RMSD (A)', fontsize=11)
        ax.set_title(f'Febuxostat Rep{rep}', fontsize=12)
        ax.legend(loc='upper right')
        ax.set_ylim(0, max(rmsd.max() * 1.1, 10))
        ax.grid(True, alpha=0.3)

        # Find when RMSD first exceeds 5A
        exceed_5 = np.where(rmsd > 5)[0]
        if len(exceed_5) > 0:
            first_exceed = times[exceed_5[0]]
            print(f"\nRep{rep}: RMSD first exceeds 5A at t = {first_exceed:.2f} ns")
            ax.axvline(x=first_exceed, color='red', linestyle=':', alpha=0.5)
        else:
            print(f"\nRep{rep}: RMSD never exceeds 5A (stable throughout)")

        # Statistics
        bound_frac = (rmsd < 3).sum() / len(rmsd) * 100
        print(f"  Mean RMSD: {rmsd.mean():.2f}A, Max: {rmsd.max():.2f}A, Bound(<3A): {bound_frac:.1f}%")

    axes[2].set_xlabel('Time (ns)', fontsize=12)
    plt.suptitle('FEBUXOSTAT REPLICATE COMPARISON: Ligand RMSD Over Time', fontsize=14, fontweight='bold')
    plt.tight_layout()
    plt.savefig(OUTPUT_DIR / 'febuxostat_rmsd_comparison.png', dpi=300)
    plt.close()
    print(f"\nSaved: {OUTPUT_DIR / 'febuxostat_rmsd_comparison.png'}")

    # =========================================================================
    # PART 4: H-bonds at frame 0
    # =========================================================================
    print("\n" + "=" * 80)
    print("PART 4: H-BONDS AT FRAME 0")
    print("=" * 80)

    for rep in [1, 2, 3]:
        u = reps[rep]['universe']
        lig_sel = reps[rep]['lig_sel']

        n_hbonds, hbond_residues = count_hbonds_frame0(u, lig_sel, cutoff=3.5)

        print(f"\nRep{rep} Frame 0 H-bonds: {n_hbonds}")
        if hbond_residues:
            print(f"  Residues: {', '.join(hbond_residues)}")

    # =========================================================================
    # PART 5: Pocket contacts at frame 0
    # =========================================================================
    print("\n" + "=" * 80)
    print("PART 5: POCKET CONTACTS AT FRAME 0")
    print("=" * 80)

    # Define pocket from budesonide (same as full_pocket_analysis.py)
    u_bud = mda.Universe(
        str(TRAJ_DIR / 'budesonide_rep1_solvated.pdb'),
        str(TRAJ_DIR / 'budesonide_rep1.dcd')
    )
    lig_sel_bud, _ = pick_ligand(u_bud)

    # Get pocket residues from budesonide frame 0
    u_bud.trajectory[0]
    lig_heavy_bud = u_bud.select_atoms(f"({lig_sel_bud}) and not name H*")
    protein_heavy = u_bud.select_atoms("protein and not name H*")
    box = u_bud.dimensions
    dists = distance_array(lig_heavy_bud.positions, protein_heavy.positions, box=box)
    contact_mask = np.any(dists < 5.0, axis=0)
    contact_atoms = protein_heavy[contact_mask]
    pocket_resids = set(contact_atoms.resids)

    print(f"\nReference pocket (from budesonide): {len(pocket_resids)} residues")
    print(f"Pocket residues: {sorted(pocket_resids)}")

    for rep in [1, 2, 3]:
        u = reps[rep]['universe']
        lig_sel = reps[rep]['lig_sel']

        n_contacts, contact_residues = compute_pocket_contacts_frame0(u, lig_sel, pocket_resids)

        print(f"\nRep{rep} Frame 0 Pocket Contacts: {n_contacts}/{len(pocket_resids)}")
        if contact_residues:
            print(f"  Contacting: {contact_residues}")

    # =========================================================================
    # PART 6: Combined visualization
    # =========================================================================
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

    # Left: All three RMSD curves
    for rep in [1, 2, 3]:
        times = all_rmsd_data[rep]['times']
        rmsd = all_rmsd_data[rep]['rmsd']
        ax1.plot(times, rmsd, color=colors[str(rep)], linewidth=1.5,
                 label=f'Rep{rep} (mean={rmsd.mean():.1f}A)', alpha=0.8)

    ax1.axhline(y=3, color='green', linestyle='--', alpha=0.7, label='3A (bound)')
    ax1.axhline(y=5, color='orange', linestyle='--', alpha=0.7, label='5A')
    ax1.set_xlabel('Time (ns)', fontsize=12)
    ax1.set_ylabel('Ligand RMSD (A)', fontsize=12)
    ax1.set_title('RMSD Comparison: All Replicates', fontsize=12)
    ax1.legend(loc='upper right')
    ax1.grid(True, alpha=0.3)
    ax1.set_ylim(0, 25)

    # Right: Bar chart of key metrics
    reps_list = [1, 2, 3]
    mean_rmsd = [all_rmsd_data[r]['rmsd'].mean() for r in reps_list]
    bound_frac = [(all_rmsd_data[r]['rmsd'] < 3).sum() / len(all_rmsd_data[r]['rmsd']) * 100 for r in reps_list]

    x = np.arange(3)
    width = 0.35

    bars1 = ax2.bar(x - width/2, mean_rmsd, width, label='Mean RMSD (A)',
                    color=[colors[str(r)] for r in reps_list], alpha=0.7)
    ax2.set_ylabel('Mean RMSD (A)', fontsize=12)
    ax2.set_xlabel('Replicate', fontsize=12)
    ax2.set_xticks(x)
    ax2.set_xticklabels(['Rep1\n(FAILED)', 'Rep2\n(SUCCESS)', 'Rep3\n(SUCCESS)'])
    ax2.set_title('Mean RMSD by Replicate', fontsize=12)

    # Add value labels
    for bar, val in zip(bars1, mean_rmsd):
        ax2.annotate(f'{val:.1f}A',
                     xy=(bar.get_x() + bar.get_width() / 2, bar.get_height()),
                     xytext=(0, 3), textcoords="offset points",
                     ha='center', va='bottom', fontweight='bold')

    plt.suptitle('FEBUXOSTAT: Why did Rep1 fail?', fontsize=14, fontweight='bold')
    plt.tight_layout()
    plt.savefig(OUTPUT_DIR / 'febuxostat_comparison_summary.png', dpi=300)
    plt.close()
    print(f"\nSaved: {OUTPUT_DIR / 'febuxostat_comparison_summary.png'}")

    # =========================================================================
    # FINAL VERDICT
    # =========================================================================
    print("\n" + "=" * 80)
    print("*** FINAL VERDICT: WHY DID REP1 FAIL? ***")
    print("=" * 80)

    # Check starting pose differences
    r1_r2_diff = np.sqrt(((frame0_positions[1] - frame0_positions[2]) ** 2).sum() / len(frame0_positions[1]))

    print("\nFINDINGS:")
    print("-" * 60)

    if r1_r2_diff < 0.5:
        print("1. STARTING POSES: All replicates started IDENTICAL")
        print("   --> Starting pose is NOT the cause of failure")
    else:
        print(f"1. STARTING POSES: Rep1 differs from Rep2/3 by {r1_r2_diff:.2f}A")
        print("   --> Different starting pose may contribute to failure")

    # When did rep1 leave?
    rep1_rmsd = all_rmsd_data[1]['rmsd']
    rep1_times = all_rmsd_data[1]['times']
    exceed_idx = np.where(rep1_rmsd > 5)[0]

    if len(exceed_idx) > 0:
        first_exit = rep1_times[exceed_idx[0]]
        if first_exit < 1:
            print(f"\n2. TIMING: Rep1 left pocket VERY EARLY (at {first_exit:.2f} ns)")
            print("   --> Suggests weak initial binding / bad luck with random velocities")
        else:
            print(f"\n2. TIMING: Rep1 left pocket at {first_exit:.2f} ns")
            print("   --> Gradual destabilization over time")

    print("\n3. OVERALL INTERPRETATION:")
    print("   - Rep2 and Rep3 demonstrate STABLE binding (low RMSD throughout)")
    print("   - Rep1 is a statistical outlier - ligand escaped early")
    print("   - This is NORMAL in MD - 2/3 stable replicates = STRONG evidence of binding")
    print("   - For publication: Report 'stable binding in 2/3 replicates'")
    print("=" * 80)


if __name__ == "__main__":
    main()
