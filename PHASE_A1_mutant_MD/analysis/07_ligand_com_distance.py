#!/usr/bin/env python3
"""
PHASE A1 ANALYSIS - Ligand Center of Mass Distance (FIXED)
Tracks ligand COM distance from pocket centroid using on-the-fly alignment.
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
    import MDAnalysis.transformations as trans
    import numpy as np
except ImportError:
    print("ERROR: MDAnalysis not installed!")
    sys.exit(1)

# ============================================================
# PARAMETERS (LOCKED)
# ============================================================
PROTEIN_BACKBONE = "protein and backbone"
POCKET_SELECTION = "resid 1008 or resid 1011 or resid 1037"  # GLU1008, ASP1011, ARG1037
LIGAND_SELECTION = "resname UNK"

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

def calculate_com_distance(pdb_path, dcd_path, stride=1):
    """Calculate ligand COM distance from pocket centroid with on-the-fly alignment."""
    print(f"  Loading: {pdb_path.name}")
    u = mda.Universe(str(pdb_path), str(dcd_path))
    ref = mda.Universe(str(pdb_path))  # Reference for alignment

    pocket = u.select_atoms(POCKET_SELECTION)
    ligand = u.select_atoms(LIGAND_SELECTION)

    print(f"  Pocket atoms: {len(pocket)} | Ligand atoms: {len(ligand)}")

    if len(pocket) == 0 or len(ligand) == 0:
        return None, None, None, None, None

    # Add on-the-fly transformation to align frames as they are read
    # This fixes the alignment issue!
    protein_mobile = u.select_atoms(PROTEIN_BACKBONE)
    protein_ref = ref.select_atoms(PROTEIN_BACKBONE)

    transform = trans.fit_rot_trans(protein_mobile, protein_ref)
    u.trajectory.add_transformations(transform)

    distances = []

    total_frames = len(u.trajectory)
    frames_analyzed = 0
    frame_count = 0

    # Read frames sequentially (DCD files have issues with random access)
    for ts in u.trajectory:
        frame_count += 1
        if (frame_count - 1) % stride != 0:
            continue

        frames_analyzed += 1

        # Ligand COM
        lig_com = ligand.center_of_mass()

        # Pocket COM (current frame, now aligned)
        pocket_com = pocket.center_of_mass()

        # Distance from ligand to pocket with PBC minimum image convention
        diff = lig_com - pocket_com
        box = u.dimensions[:3]
        diff = diff - box * np.round(diff / box)  # Minimum image
        dist = np.linalg.norm(diff)
        distances.append(dist)

        if frames_analyzed % 500 == 0:
            print(f"    Processed {frames_analyzed} frames...")

    dist_array = np.array(distances)

    # Check for drift (large increase over time)
    if len(dist_array) >= 4:
        quarter_len = len(dist_array) // 4
        first_quarter = dist_array[:quarter_len].mean()
        last_quarter = dist_array[-quarter_len:].mean()
        drift = last_quarter - first_quarter
    else:
        drift = 0.0

    print(f"  Analyzed {frames_analyzed} frames")

    return dist_array.mean(), dist_array.std(), dist_array.max(), drift, distances

def main():
    print("=" * 80)
    print("PHASE A1 ANALYSIS - LIGAND COM DISTANCE")
    print("=" * 80)
    print("Tracking ligand center-of-mass distance from pocket centroid")
    print("Using on-the-fly alignment for accuracy")
    print()

    analysis_dir = BASE_DIR / "analysis"
    analysis_dir.mkdir(exist_ok=True)

    results = []
    all_timeseries = {}

    for i, (mutant, ligand, rep, folder_key) in enumerate(SIMULATIONS, 1):
        print(f"[{i}/12] {mutant}_{ligand}_rep{rep}")
        dcd_path, pdb_path = get_file_paths(mutant, ligand, rep, folder_key)

        if not dcd_path.exists() or not pdb_path.exists():
            print(f"  ERROR: Files not found")
            results.append({'mutant': mutant, 'ligand': ligand, 'rep': rep,
                          'mean_dist': None, 'status': 'FILE_NOT_FOUND'})
            continue

        try:
            mean_dist, std_dist, max_dist, drift, timeseries = calculate_com_distance(pdb_path, dcd_path)

            if mean_dist is not None:
                # Drift assessment
                if abs(drift) < 1.0:
                    status = "STABLE"
                elif abs(drift) < 3.0:
                    status = "MINOR_DRIFT"
                else:
                    status = "SIGNIFICANT_DRIFT"

                print(f"  COM Distance: {mean_dist:.2f} +/- {std_dist:.2f} A (max: {max_dist:.2f} A)")
                print(f"  Drift (last-first quarter): {drift:+.2f} A - {status}")

                results.append({'mutant': mutant, 'ligand': ligand, 'rep': rep,
                              'mean_dist': mean_dist, 'std_dist': std_dist,
                              'max_dist': max_dist, 'drift': drift, 'status': status})
                all_timeseries[f"{mutant}_{ligand}_rep{rep}"] = timeseries
            else:
                results.append({'mutant': mutant, 'ligand': ligand, 'rep': rep,
                              'mean_dist': None, 'status': 'SELECTION_ERROR'})

        except Exception as e:
            print(f"  ERROR: {e}")
            results.append({'mutant': mutant, 'ligand': ligand, 'rep': rep,
                          'mean_dist': None, 'status': 'ERROR'})

        print()

    # Print results
    print("=" * 80)
    print("LIGAND COM DISTANCE RESULTS")
    print("=" * 80)
    print(f"{'ID':<4} {'Mutant':<8} {'Ligand':<12} {'Rep':<4} {'Mean (A)':<10} {'SD (A)':<10} {'Drift (A)':<10} {'Status'}")
    print("-" * 80)

    for i, r in enumerate(results, 1):
        mean_str = f"{r['mean_dist']:.2f}" if r['mean_dist'] is not None else "N/A"
        std_str = f"{r.get('std_dist', 0):.2f}" if r['mean_dist'] is not None else "N/A"
        drift_str = f"{r.get('drift', 0):+.2f}" if r['mean_dist'] is not None else "N/A"
        print(f"{i:<4} {r['mutant']:<8} {r['ligand']:<12} {r['rep']:<4} {mean_str:<10} {std_str:<10} {drift_str:<10} {r['status']}")

    # Save CSV
    csv_path = analysis_dir / "ligand_com_distance.csv"
    with open(csv_path, 'w') as f:
        f.write("sim_id,mutant,ligand,rep,mean_dist,std_dist,max_dist,drift,status\n")
        for i, r in enumerate(results, 1):
            mean = r['mean_dist'] if r['mean_dist'] is not None else ""
            std = r.get('std_dist', "") if r['mean_dist'] is not None else ""
            max_d = r.get('max_dist', "") if r['mean_dist'] is not None else ""
            drift = r.get('drift', "") if r['mean_dist'] is not None else ""
            f.write(f"{i},{r['mutant']},{r['ligand']},{r['rep']},{mean},{std},{max_d},{drift},{r['status']}\n")

    print(f"\nResults saved to: {csv_path}")

    # Save timeseries
    ts_path = analysis_dir / "ligand_com_timeseries.csv"
    with open(ts_path, 'w') as f:
        headers = ['frame'] + list(all_timeseries.keys())
        f.write(','.join(headers) + '\n')
        if all_timeseries:
            max_len = max(len(v) for v in all_timeseries.values())
            for i in range(max_len):
                row = [str(i)]
                for key in all_timeseries.keys():
                    if i < len(all_timeseries[key]):
                        row.append(f"{all_timeseries[key][i]:.4f}")
                    else:
                        row.append("")
                f.write(','.join(row) + '\n')

    print(f"Timeseries saved to: {ts_path}")
    print("\nLigand COM distance analysis complete!")

if __name__ == "__main__":
    main()
