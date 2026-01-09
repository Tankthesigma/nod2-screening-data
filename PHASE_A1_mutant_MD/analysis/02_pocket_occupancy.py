#!/usr/bin/env python3
"""
PHASE A1 ANALYSIS - CHECKPOINT 2: Pocket Occupancy
Calculates ligand pocket occupancy for all 12 simulations.

Pocket Definition:
- Residues: GLU1008, ASP1011, ARG1037 (heavy atoms only)
- Cutoff: 5.0 Angstroms
- Occupied = ANY ligand heavy atom within cutoff of ANY pocket heavy atom
"""

import os
import sys
from pathlib import Path
import platform
import warnings
warnings.filterwarnings('ignore')

# Detect if running in WSL or Windows
if 'microsoft' in platform.uname().release.lower() or os.path.exists('/mnt/c'):
    BASE_DIR = Path("/mnt/c/Users/vasud/nod2-screening-data/PHASE_A1_mutant_MD")
else:
    BASE_DIR = Path(r"C:\Users\vasud\nod2-screening-data\PHASE_A1_mutant_MD")

# Check for MDAnalysis
try:
    import MDAnalysis as mda
    from MDAnalysis.analysis import distances
    import numpy as np
except ImportError:
    print("ERROR: MDAnalysis not installed!")
    print("Run: pip install MDAnalysis")
    sys.exit(1)

# ============================================================
# ANALYSIS PARAMETERS (LOCKED - DO NOT CHANGE)
# ============================================================
POCKET_SELECTION = "resid 1008 or resid 1011 or resid 1037"  # GLU1008, ASP1011, ARG1037
LIGAND_SELECTION = "resname UNK"                              # Ligand
CUTOFF = 5.0  # Angstroms

# All trajectory locations
TRAJECTORY_LOCATIONS = {
    'trajectories': BASE_DIR / "trajectories",
    '5080': BASE_DIR / "vast_downloads" / "5080",
    '5090': BASE_DIR / "vast_downloads" / "5090",
    'new_5090': BASE_DIR / "vast_downloads" / "new_5090",
}

# ============================================================
# FILE MAPPING (from checkpoint 1)
# ============================================================
SIMULATIONS = [
    # (mutant, ligand, rep, folder_key)
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
    """Get DCD and PDB paths for a simulation."""
    folder = TRAJECTORY_LOCATIONS[folder_key]
    base_name = f"{mutant}_{ligand}_rep{rep}"
    dcd_path = folder / f"{base_name}.dcd"
    pdb_path = folder / f"{base_name}_solvated.pdb"
    return dcd_path, pdb_path

def calculate_pocket_occupancy(pdb_path, dcd_path, stride=1):
    """
    Calculate pocket occupancy for a single trajectory.

    Returns:
        occupancy_percent: float (0-100)
        occupied_frames: int
        total_frames: int
        min_distance_per_frame: list of floats
    """
    print(f"  Loading topology: {pdb_path.name}")
    u = mda.Universe(str(pdb_path), str(dcd_path))

    # Select atoms
    pocket = u.select_atoms(POCKET_SELECTION)
    ligand = u.select_atoms(LIGAND_SELECTION)

    print(f"  Pocket atoms: {len(pocket)} | Ligand atoms: {len(ligand)}")

    if len(pocket) == 0:
        print(f"  ERROR: No pocket atoms found with selection: {POCKET_SELECTION}")
        return None, None, None, None

    if len(ligand) == 0:
        print(f"  ERROR: No ligand atoms found with selection: {LIGAND_SELECTION}")
        return None, None, None, None

    # Analyze trajectory
    total_frames = len(u.trajectory)
    expected_frames = (total_frames + stride - 1) // stride

    print(f"  Total frames: {total_frames} | Analyzing: ~{expected_frames} (stride={stride})")

    occupied_count = 0
    min_distances = []
    frame_count = 0
    analyzed_count = 0

    # Read frames sequentially (DCD files have issues with random access)
    for ts in u.trajectory:
        frame_count += 1
        if (frame_count - 1) % stride != 0:
            continue

        analyzed_count += 1

        # Calculate all pairwise distances between ligand and pocket
        # Use box=u.dimensions for PBC minimum image convention
        dist_matrix = distances.distance_array(ligand.positions, pocket.positions, box=u.dimensions)

        # Minimum distance in this frame
        min_dist = dist_matrix.min()
        min_distances.append(min_dist)

        # Frame is occupied if any distance < cutoff
        if min_dist < CUTOFF:
            occupied_count += 1

        # Progress update every 500 analyzed frames
        if analyzed_count % 500 == 0:
            print(f"    Processed {analyzed_count} frames...")

    occupancy_percent = (occupied_count / analyzed_count) * 100

    return occupancy_percent, occupied_count, analyzed_count, min_distances

def main():
    print("=" * 80)
    print("PHASE A1 ANALYSIS - CHECKPOINT 2: POCKET OCCUPANCY")
    print("=" * 80)
    print()
    print("Parameters (LOCKED):")
    print(f"  Pocket residues: GLU1008, ASP1011, ARG1037")
    print(f"  Pocket selection: {POCKET_SELECTION}")
    print(f"  Ligand selection: {LIGAND_SELECTION}")
    print(f"  Cutoff: {CUTOFF} Å")
    print()

    # Create analysis directory
    analysis_dir = BASE_DIR / "analysis"
    analysis_dir.mkdir(exist_ok=True)

    # Results storage
    results = []

    # Process each simulation
    for i, (mutant, ligand, rep, folder_key) in enumerate(SIMULATIONS, 1):
        print("-" * 60)
        print(f"[{i}/12] {mutant}_{ligand}_rep{rep}")
        print("-" * 60)

        dcd_path, pdb_path = get_file_paths(mutant, ligand, rep, folder_key)

        if not dcd_path.exists():
            print(f"  ERROR: DCD not found: {dcd_path}")
            results.append({
                'mutant': mutant, 'ligand': ligand, 'rep': rep,
                'occupancy': None, 'status': 'FILE_NOT_FOUND'
            })
            continue

        if not pdb_path.exists():
            print(f"  ERROR: PDB not found: {pdb_path}")
            results.append({
                'mutant': mutant, 'ligand': ligand, 'rep': rep,
                'occupancy': None, 'status': 'FILE_NOT_FOUND'
            })
            continue

        try:
            # Use stride=1 for full analysis (2000 frames is manageable)
            occupancy, occupied, total, min_dists = calculate_pocket_occupancy(
                pdb_path, dcd_path, stride=1
            )

            if occupancy is not None:
                # Determine status
                if occupancy >= 95:
                    status = "EXCELLENT"
                elif occupancy >= 80:
                    status = "GOOD"
                elif occupancy >= 50:
                    status = "PARTIAL"
                else:
                    status = "POOR"

                avg_min_dist = np.mean(min_dists)
                max_min_dist = np.max(min_dists)

                print(f"  RESULT: {occupancy:.1f}% occupancy ({occupied}/{total} frames)")
                print(f"  Avg min distance: {avg_min_dist:.2f} Å | Max: {max_min_dist:.2f} Å")
                print(f"  Status: {status}")

                results.append({
                    'mutant': mutant, 'ligand': ligand, 'rep': rep,
                    'occupancy': occupancy, 'occupied_frames': occupied,
                    'total_frames': total, 'avg_min_dist': avg_min_dist,
                    'max_min_dist': max_min_dist, 'status': status
                })
            else:
                results.append({
                    'mutant': mutant, 'ligand': ligand, 'rep': rep,
                    'occupancy': None, 'status': 'SELECTION_ERROR'
                })

        except Exception as e:
            print(f"  ERROR: {str(e)}")
            results.append({
                'mutant': mutant, 'ligand': ligand, 'rep': rep,
                'occupancy': None, 'status': f'ERROR: {str(e)[:50]}'
            })

        print()

    # ============================================================
    # RESULTS SUMMARY
    # ============================================================
    print("=" * 80)
    print("POCKET OCCUPANCY RESULTS")
    print("=" * 80)
    print()

    # Table 1: Per-simulation results
    print("Table 1: Individual Simulation Results")
    print("-" * 90)
    print(f"{'ID':<4} {'Mutant':<8} {'Ligand':<12} {'Rep':<4} {'Occupancy':<12} {'Avg Dist (Å)':<14} {'Status':<12}")
    print("-" * 90)

    for i, r in enumerate(results, 1):
        occ_str = f"{r['occupancy']:.1f}%" if r['occupancy'] is not None else "N/A"
        dist_str = f"{r.get('avg_min_dist', 0):.2f}" if r['occupancy'] is not None else "N/A"
        print(f"{i:<4} {r['mutant']:<8} {r['ligand']:<12} {r['rep']:<4} {occ_str:<12} {dist_str:<14} {r['status']:<12}")

    print("-" * 90)
    print()

    # Table 2: Aggregated by mutant/ligand
    print("Table 2: Aggregated Results (Mean ± SD)")
    print("-" * 100)
    print(f"{'Mutant':<8} {'Ligand':<12} {'Rep1':<10} {'Rep2':<10} {'Rep3':<10} {'Avg ± SD':<15} {'WT Baseline':<12} {'Verdict':<10}")
    print("-" * 100)

    # Group results
    groups = {}
    for r in results:
        key = (r['mutant'], r['ligand'])
        if key not in groups:
            groups[key] = {1: None, 2: None, 3: None}
        groups[key][r['rep']] = r['occupancy']

    for (mutant, ligand), reps in sorted(groups.items()):
        rep1 = f"{reps[1]:.1f}%" if reps[1] is not None else "N/A"
        rep2 = f"{reps[2]:.1f}%" if reps[2] is not None else "N/A"
        rep3 = f"{reps[3]:.1f}%" if reps[3] is not None else "N/A"

        # Calculate mean and SD
        valid_vals = [v for v in reps.values() if v is not None]
        if len(valid_vals) >= 2:
            mean_occ = np.mean(valid_vals)
            std_occ = np.std(valid_vals)
            avg_str = f"{mean_occ:.1f} ± {std_occ:.1f}%"
        elif len(valid_vals) == 1:
            avg_str = f"{valid_vals[0]:.1f}%"
        else:
            avg_str = "N/A"

        # WT baseline reference
        if ligand == "febuxostat":
            wt_baseline = "99-100%"
        else:
            wt_baseline = "100%"

        # Verdict
        if len(valid_vals) >= 2:
            mean_val = np.mean(valid_vals)
            if mean_val >= 95:
                verdict = "PRESERVED"
            elif mean_val >= 80:
                verdict = "PARTIAL"
            else:
                verdict = "IMPAIRED"
        else:
            verdict = "?"

        print(f"{mutant:<8} {ligand:<12} {rep1:<10} {rep2:<10} {rep3:<10} {avg_str:<15} {wt_baseline:<12} {verdict:<10}")

    print("-" * 100)
    print()

    # Success criteria reminder
    print("Success Criteria:")
    print("  >=95% = Binding preserved")
    print("  80-95% = Partial binding")
    print("  <80% = Binding compromised")
    print()

    # Save results to CSV
    csv_path = analysis_dir / "pocket_occupancy.csv"
    with open(csv_path, 'w') as f:
        f.write("sim_id,mutant,ligand,rep,occupancy_percent,occupied_frames,total_frames,avg_min_dist,max_min_dist,status\n")
        for i, r in enumerate(results, 1):
            occ = r['occupancy'] if r['occupancy'] is not None else ""
            occ_frames = r.get('occupied_frames', "")
            tot_frames = r.get('total_frames', "")
            avg_dist = r.get('avg_min_dist', "")
            max_dist = r.get('max_min_dist', "")
            f.write(f"{i},{r['mutant']},{r['ligand']},{r['rep']},{occ},{occ_frames},{tot_frames},{avg_dist},{max_dist},{r['status']}\n")

    print(f"Results saved to: {csv_path}")
    print()
    print("=" * 80)
    print("CHECKPOINT 2 COMPLETE - WAITING FOR APPROVAL")
    print("=" * 80)

    return results

if __name__ == "__main__":
    results = main()
