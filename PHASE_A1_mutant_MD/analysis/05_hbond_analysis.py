#!/usr/bin/env python3
"""
PHASE A1 ANALYSIS - Hydrogen Bond Analysis (FIXED)
Detects and tracks H-bonds between protein and ligand.
Fixed frame counting for accurate occupancy calculation.
"""

import os
import sys
from pathlib import Path
import platform
import warnings
warnings.filterwarnings('ignore')

if 'microsoft' in platform.uname().release.lower() or os.path.exists('/mnt/c'):
    BASE_DIR = Path("/mnt/c/Users/vasud/nod2-screening-data/PHASE_A1_mutant_MD")
else:
    BASE_DIR = Path(r"C:\Users\vasud\nod2-screening-data\PHASE_A1_mutant_MD")

try:
    import MDAnalysis as mda
    from MDAnalysis.analysis.hydrogenbonds import HydrogenBondAnalysis
    import numpy as np
except ImportError:
    print("ERROR: MDAnalysis not installed!")
    sys.exit(1)

# ============================================================
# PARAMETERS (LOCKED)
# ============================================================
HBOND_DISTANCE = 3.5  # Angstroms (donor-acceptor)
HBOND_ANGLE = 150     # Degrees (D-H...A angle minimum)

LIGAND_SELECTION = "resname UNK"
PROTEIN_SELECTION = "protein"

TRAJECTORY_LOCATIONS = {
    'trajectories': BASE_DIR / "trajectories",
    '5080': BASE_DIR / "vast_downloads" / "5080",
    '5090': BASE_DIR / "vast_downloads" / "5090",
    'new_5090': BASE_DIR / "vast_downloads" / "new_5090",
}

SIMULATIONS = [
    ("G908R", "febuxostat", 1, "trajectories"),
    ("G908R", "febuxostat", 2, "5080"),
    ("G908R", "febuxostat", 3, "5080"),
    ("G908R", "natural", 1, "5090"),
    ("G908R", "natural", 2, "new_5090"),
    ("G908R", "natural", 3, "new_5090"),
    ("R702W", "febuxostat", 1, "trajectories"),
    ("R702W", "febuxostat", 2, "trajectories"),
    ("R702W", "febuxostat", 3, "trajectories"),
    ("R702W", "natural", 1, "5080"),
    ("R702W", "natural", 2, "5080"),
    ("R702W", "natural", 3, "new_5090"),
]

def get_file_paths(mutant, ligand, rep, folder_key):
    folder = TRAJECTORY_LOCATIONS[folder_key]
    base_name = f"{mutant}_{ligand}_rep{rep}"
    return folder / f"{base_name}.dcd", folder / f"{base_name}_solvated.pdb"

def analyze_hbonds(pdb_path, dcd_path, stride=10):
    """Analyze hydrogen bonds between protein and ligand."""
    print(f"  Loading: {pdb_path.name}")
    u = mda.Universe(str(pdb_path), str(dcd_path))

    total_frames = len(u.trajectory)
    # Correct frame count calculation
    frames_analyzed = len(range(0, total_frames, stride))
    print(f"  Total frames: {total_frames} | Analyzing: {frames_analyzed} (stride={stride})")

    if frames_analyzed == 0:
        return [], 0, 0

    # H-bond analysis: protein as donors, ligand as acceptors
    all_hbonds = []

    try:
        hbonds_prot_donor = HydrogenBondAnalysis(
            universe=u,
            donors_sel=PROTEIN_SELECTION,
            acceptors_sel=LIGAND_SELECTION,
            d_a_cutoff=HBOND_DISTANCE,
            d_h_a_angle_cutoff=HBOND_ANGLE,
        )
        hbonds_prot_donor.run(step=stride, verbose=False)
        prot_donor_results = hbonds_prot_donor.results.hbonds
    except Exception as e:
        print(f"  Warning (prot->lig): {e}")
        prot_donor_results = np.array([])

    # H-bond analysis: ligand as donors, protein as acceptors
    try:
        hbonds_lig_donor = HydrogenBondAnalysis(
            universe=u,
            donors_sel=LIGAND_SELECTION,
            acceptors_sel=PROTEIN_SELECTION,
            d_a_cutoff=HBOND_DISTANCE,
            d_h_a_angle_cutoff=HBOND_ANGLE,
        )
        hbonds_lig_donor.run(step=stride, verbose=False)
        lig_donor_results = hbonds_lig_donor.results.hbonds
    except Exception as e:
        print(f"  Warning (lig->prot): {e}")
        lig_donor_results = np.array([])

    # Process protein donor H-bonds
    if len(prot_donor_results) > 0:
        for hb in prot_donor_results:
            frame, donor_idx, hydrogen_idx, acceptor_idx, distance, angle = hb
            donor_atom = u.atoms[int(donor_idx)]
            acceptor_atom = u.atoms[int(acceptor_idx)]
            all_hbonds.append({
                'frame': int(frame),
                'donor': f"{donor_atom.resname}{donor_atom.resid}:{donor_atom.name}",
                'acceptor': f"{acceptor_atom.resname}:{acceptor_atom.name}",
                'distance': distance,
                'angle': angle,
                'type': 'protein->ligand'
            })

    # Process ligand donor H-bonds
    if len(lig_donor_results) > 0:
        for hb in lig_donor_results:
            frame, donor_idx, hydrogen_idx, acceptor_idx, distance, angle = hb
            donor_atom = u.atoms[int(donor_idx)]
            acceptor_atom = u.atoms[int(acceptor_idx)]
            all_hbonds.append({
                'frame': int(frame),
                'donor': f"{donor_atom.resname}:{donor_atom.name}",
                'acceptor': f"{acceptor_atom.resname}{acceptor_atom.resid}:{acceptor_atom.name}",
                'distance': distance,
                'angle': angle,
                'type': 'ligand->protein'
            })

    # Calculate occupancies - count UNIQUE FRAMES per pair (not raw counts)
    hbond_pairs = {}
    for hb in all_hbonds:
        key = (hb['donor'], hb['acceptor'])
        if key not in hbond_pairs:
            hbond_pairs[key] = {'frames': set(), 'type': hb['type']}
        hbond_pairs[key]['frames'].add(hb['frame'])

    # Convert to occupancy percentages
    hbond_occupancies = []
    for (donor, acceptor), data in hbond_pairs.items():
        # Occupancy = unique frames with this H-bond / total frames analyzed
        occupancy = (len(data['frames']) / frames_analyzed) * 100
        hbond_occupancies.append({
            'donor': donor,
            'acceptor': acceptor,
            'occupancy': min(occupancy, 100.0),  # Cap at 100%
            'type': data['type'],
            'frame_count': len(data['frames'])
        })

    # Sort by occupancy
    hbond_occupancies.sort(key=lambda x: x['occupancy'], reverse=True)

    return hbond_occupancies, len(all_hbonds), frames_analyzed

def main():
    print("=" * 80)
    print("PHASE A1 ANALYSIS - HYDROGEN BONDS")
    print("=" * 80)
    print(f"H-bond cutoff: {HBOND_DISTANCE} A")
    print(f"H-bond angle: >={HBOND_ANGLE} deg")
    print()

    analysis_dir = BASE_DIR / "analysis"
    analysis_dir.mkdir(exist_ok=True)

    all_results = []
    summary_data = []

    for i, (mutant, ligand, rep, folder_key) in enumerate(SIMULATIONS, 1):
        print(f"[{i}/12] {mutant}_{ligand}_rep{rep}")
        dcd_path, pdb_path = get_file_paths(mutant, ligand, rep, folder_key)

        if not dcd_path.exists() or not pdb_path.exists():
            print(f"  ERROR: Files not found")
            summary_data.append({'mutant': mutant, 'ligand': ligand, 'rep': rep,
                               'total_hbonds': 0, 'unique_pairs': 0, 'status': 'FILE_NOT_FOUND'})
            continue

        try:
            hbond_occupancies, total_hbonds, frames = analyze_hbonds(pdb_path, dcd_path, stride=10)

            unique_pairs = len(hbond_occupancies)
            persistent = len([h for h in hbond_occupancies if h['occupancy'] > 25])

            print(f"  Total H-bond events: {total_hbonds}")
            print(f"  Unique pairs: {unique_pairs} | Persistent (>25%): {persistent}")

            summary_data.append({'mutant': mutant, 'ligand': ligand, 'rep': rep,
                               'total_hbonds': total_hbonds, 'unique_pairs': unique_pairs,
                               'persistent': persistent, 'frames': frames, 'status': 'OK'})

            # Store detailed results
            for hb in hbond_occupancies:
                all_results.append({
                    'mutant': mutant, 'ligand': ligand, 'rep': rep,
                    'donor': hb['donor'], 'acceptor': hb['acceptor'],
                    'occupancy': hb['occupancy'], 'type': hb['type']
                })

            # Print top H-bonds
            if hbond_occupancies:
                print(f"  Top H-bonds:")
                for hb in hbond_occupancies[:5]:
                    print(f"    {hb['donor']} -> {hb['acceptor']}: {hb['occupancy']:.1f}%")

        except Exception as e:
            print(f"  ERROR: {e}")
            summary_data.append({'mutant': mutant, 'ligand': ligand, 'rep': rep,
                               'total_hbonds': 0, 'unique_pairs': 0, 'status': f'ERROR'})

        print()

    # Save detailed results
    csv_path = analysis_dir / "hbond_analysis.csv"
    with open(csv_path, 'w') as f:
        f.write("mutant,ligand,rep,donor,acceptor,occupancy_percent,type\n")
        for r in all_results:
            # Escape commas in donor/acceptor names
            donor = r['donor'].replace(',', ';')
            acceptor = r['acceptor'].replace(',', ';')
            f.write(f"{r['mutant']},{r['ligand']},{r['rep']},{donor},{acceptor},{r['occupancy']:.2f},{r['type']}\n")
    print(f"Detailed H-bonds saved to: {csv_path}")

    # Print summary
    print("\n" + "=" * 80)
    print("H-BOND SUMMARY")
    print("=" * 80)
    print(f"{'Mutant':<8} {'Ligand':<12} {'Rep':<4} {'Total':<8} {'Unique':<8} {'Persistent':<10} {'Status'}")
    print("-" * 70)
    for s in summary_data:
        print(f"{s['mutant']:<8} {s['ligand']:<12} {s['rep']:<4} {s.get('total_hbonds', 'N/A'):<8} {s.get('unique_pairs', 'N/A'):<8} {s.get('persistent', 'N/A'):<10} {s['status']}")

    # Save summary
    summary_path = analysis_dir / "hbond_summary.csv"
    with open(summary_path, 'w') as f:
        f.write("mutant,ligand,rep,total_hbonds,unique_pairs,persistent,frames,status\n")
        for s in summary_data:
            f.write(f"{s['mutant']},{s['ligand']},{s['rep']},{s.get('total_hbonds','')},{s.get('unique_pairs','')},{s.get('persistent','')},{s.get('frames','')},{s['status']}\n")

    print(f"\nSummary saved to: {summary_path}")
    print("\nH-bond analysis complete!")

if __name__ == "__main__":
    main()
