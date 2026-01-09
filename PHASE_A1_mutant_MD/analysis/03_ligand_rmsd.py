#!/usr/bin/env python3
"""
PHASE A1 ANALYSIS - Ligand RMSD (FIXED)
Calculates ligand RMSD vs first frame after backbone alignment.
Uses rms.RMSD class for proper alignment.
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
    from MDAnalysis.analysis import rms
    import numpy as np
except ImportError:
    print("ERROR: MDAnalysis not installed!")
    sys.exit(1)

# ============================================================
# PARAMETERS (LOCKED)
# ============================================================
PROTEIN_SELECTION = "protein and backbone"  # For alignment
LIGAND_SELECTION = "resname UNK"  # Ligand (LOCKED - no H* filter here, RMSD class handles it)

# File mapping
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

def calculate_ligand_rmsd(pdb_path, dcd_path, stride=1):
    """Calculate ligand RMSD using rms.RMSD class (proper alignment)."""
    print(f"  Loading: {pdb_path.name}")
    u = mda.Universe(str(pdb_path), str(dcd_path))
    ref = mda.Universe(str(pdb_path))  # Just topology for reference (frame 0)

    # Check ligand exists
    ligand = u.select_atoms(LIGAND_SELECTION)
    print(f"  Ligand atoms: {len(ligand)}")

    if len(ligand) == 0:
        return None, None, None, None

    # Use rms.RMSD class - this properly aligns and calculates RMSD
    # Column 0: frame, Column 1: time, Column 2: backbone RMSD, Column 3: ligand RMSD
    R = rms.RMSD(u, ref,
                 select=PROTEIN_SELECTION,  # Align on backbone
                 groupselections=[LIGAND_SELECTION])  # Calculate RMSD for ligand
    R.run(step=stride, verbose=False)

    # Extract ligand RMSD (column 3)
    rmsd_values = R.results.rmsd[:, 3]

    print(f"  Analyzed {len(rmsd_values)} frames")

    return rmsd_values.mean(), rmsd_values.std(), rmsd_values.max(), rmsd_values.tolist()

def main():
    print("=" * 80)
    print("PHASE A1 ANALYSIS - LIGAND RMSD")
    print("=" * 80)
    print(f"Alignment: {PROTEIN_SELECTION}")
    print(f"RMSD of: {LIGAND_SELECTION}")
    print()

    analysis_dir = BASE_DIR / "analysis"
    plots_dir = analysis_dir / "plots"
    plots_dir.mkdir(parents=True, exist_ok=True)

    results = []
    all_timeseries = {}

    for i, (mutant, ligand, rep, folder_key) in enumerate(SIMULATIONS, 1):
        print(f"[{i}/12] {mutant}_{ligand}_rep{rep}")
        dcd_path, pdb_path = get_file_paths(mutant, ligand, rep, folder_key)

        if not dcd_path.exists() or not pdb_path.exists():
            print(f"  ERROR: Files not found")
            results.append({'mutant': mutant, 'ligand': ligand, 'rep': rep,
                          'mean_rmsd': None, 'std_rmsd': None, 'max_rmsd': None, 'status': 'FILE_NOT_FOUND'})
            continue

        try:
            mean_rmsd, std_rmsd, max_rmsd, timeseries = calculate_ligand_rmsd(pdb_path, dcd_path)

            if mean_rmsd is not None:
                print(f"  RMSD: {mean_rmsd:.2f} ± {std_rmsd:.2f} Å (max: {max_rmsd:.2f} Å)")

                # Stability assessment
                if mean_rmsd < 2.0:
                    status = "STABLE"
                elif mean_rmsd < 4.0:
                    status = "MODERATE"
                else:
                    status = "UNSTABLE"

                results.append({'mutant': mutant, 'ligand': ligand, 'rep': rep,
                              'mean_rmsd': mean_rmsd, 'std_rmsd': std_rmsd, 'max_rmsd': max_rmsd, 'status': status})
                all_timeseries[f"{mutant}_{ligand}_rep{rep}"] = timeseries
            else:
                results.append({'mutant': mutant, 'ligand': ligand, 'rep': rep,
                              'mean_rmsd': None, 'std_rmsd': None, 'max_rmsd': None, 'status': 'SELECTION_ERROR'})
        except Exception as e:
            print(f"  ERROR: {e}")
            results.append({'mutant': mutant, 'ligand': ligand, 'rep': rep,
                          'mean_rmsd': None, 'std_rmsd': None, 'max_rmsd': None, 'status': f'ERROR'})
        print()

    # Print results table
    print("=" * 80)
    print("LIGAND RMSD RESULTS")
    print("=" * 80)
    print(f"{'ID':<4} {'Mutant':<8} {'Ligand':<12} {'Rep':<4} {'Mean (Å)':<10} {'SD (Å)':<10} {'Max (Å)':<10} {'Status':<10}")
    print("-" * 80)
    for i, r in enumerate(results, 1):
        mean_str = f"{r['mean_rmsd']:.2f}" if r['mean_rmsd'] is not None else "N/A"
        std_str = f"{r['std_rmsd']:.2f}" if r['std_rmsd'] is not None else "N/A"
        max_str = f"{r['max_rmsd']:.2f}" if r['max_rmsd'] is not None else "N/A"
        print(f"{i:<4} {r['mutant']:<8} {r['ligand']:<12} {r['rep']:<4} {mean_str:<10} {std_str:<10} {max_str:<10} {r['status']:<10}")

    # Save CSV
    csv_path = analysis_dir / "ligand_rmsd.csv"
    with open(csv_path, 'w') as f:
        f.write("sim_id,mutant,ligand,rep,mean_rmsd,std_rmsd,max_rmsd,status\n")
        for i, r in enumerate(results, 1):
            mean = r['mean_rmsd'] if r['mean_rmsd'] is not None else ""
            std = r['std_rmsd'] if r['std_rmsd'] is not None else ""
            max_r = r['max_rmsd'] if r['max_rmsd'] is not None else ""
            f.write(f"{i},{r['mutant']},{r['ligand']},{r['rep']},{mean},{std},{max_r},{r['status']}\n")

    print(f"\nResults saved to: {csv_path}")

    # Save timeseries for plotting
    ts_path = analysis_dir / "ligand_rmsd_timeseries.csv"
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
    print("\nLigand RMSD analysis complete!")

    return results

if __name__ == "__main__":
    main()
