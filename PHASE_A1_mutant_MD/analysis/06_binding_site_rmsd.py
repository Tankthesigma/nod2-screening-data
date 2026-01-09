#!/usr/bin/env python3
"""
PHASE A1 ANALYSIS - Binding Site RMSD (FIXED)
Calculates RMSD of pocket residues using rms.RMSD class for proper alignment.
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
PROTEIN_BACKBONE = "protein and backbone"  # For alignment
POCKET_SELECTION = "resid 1008 or resid 1011 or resid 1037"  # GLU1008, ASP1011, ARG1037

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

def calculate_binding_site_rmsd(pdb_path, dcd_path, stride=1):
    """Calculate RMSD of binding site using rms.RMSD class (proper alignment)."""
    print(f"  Loading: {pdb_path.name}")
    u = mda.Universe(str(pdb_path), str(dcd_path))
    ref = mda.Universe(str(pdb_path))  # Just topology for reference

    pocket = u.select_atoms(POCKET_SELECTION)
    print(f"  Pocket atoms: {len(pocket)}")

    if len(pocket) == 0:
        return None, None, None, None

    # Use rms.RMSD class - properly aligns and calculates RMSD
    # Column 0: frame, Column 1: time, Column 2: backbone RMSD, Column 3: pocket RMSD
    R = rms.RMSD(u, ref,
                 select=PROTEIN_BACKBONE,  # Align on backbone
                 groupselections=[POCKET_SELECTION])  # Calculate RMSD for pocket
    R.run(step=stride, verbose=False)

    # Extract pocket RMSD (column 3)
    rmsd_values = R.results.rmsd[:, 3]

    print(f"  Analyzed {len(rmsd_values)} frames")

    return rmsd_values.mean(), rmsd_values.std(), rmsd_values.max(), rmsd_values.tolist()

def main():
    print("=" * 80)
    print("PHASE A1 ANALYSIS - BINDING SITE RMSD")
    print("=" * 80)
    print(f"Pocket residues: GLU1008, ASP1011, ARG1037")
    print(f"Selection: {POCKET_SELECTION}")
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
                          'mean_rmsd': None, 'status': 'FILE_NOT_FOUND'})
            continue

        try:
            mean_rmsd, std_rmsd, max_rmsd, timeseries = calculate_binding_site_rmsd(pdb_path, dcd_path)

            if mean_rmsd is not None:
                # Stability assessment
                if mean_rmsd < 1.5:
                    status = "VERY_STABLE"
                elif mean_rmsd < 2.5:
                    status = "STABLE"
                elif mean_rmsd < 4.0:
                    status = "MODERATE"
                else:
                    status = "UNSTABLE"

                print(f"  Pocket RMSD: {mean_rmsd:.2f} ± {std_rmsd:.2f} Å (max: {max_rmsd:.2f} Å) - {status}")

                results.append({'mutant': mutant, 'ligand': ligand, 'rep': rep,
                              'mean_rmsd': mean_rmsd, 'std_rmsd': std_rmsd,
                              'max_rmsd': max_rmsd, 'status': status})
                all_timeseries[f"{mutant}_{ligand}_rep{rep}"] = timeseries
            else:
                results.append({'mutant': mutant, 'ligand': ligand, 'rep': rep,
                              'mean_rmsd': None, 'status': 'SELECTION_ERROR'})

        except Exception as e:
            print(f"  ERROR: {e}")
            results.append({'mutant': mutant, 'ligand': ligand, 'rep': rep,
                          'mean_rmsd': None, 'status': 'ERROR'})

        print()

    # Print results
    print("=" * 80)
    print("BINDING SITE RMSD RESULTS")
    print("=" * 80)
    print(f"{'ID':<4} {'Mutant':<8} {'Ligand':<12} {'Rep':<4} {'Mean (Å)':<10} {'SD (Å)':<10} {'Max (Å)':<10} {'Status'}")
    print("-" * 80)

    for i, r in enumerate(results, 1):
        mean_str = f"{r['mean_rmsd']:.2f}" if r['mean_rmsd'] is not None else "N/A"
        std_str = f"{r.get('std_rmsd', 0):.2f}" if r['mean_rmsd'] is not None else "N/A"
        max_str = f"{r.get('max_rmsd', 0):.2f}" if r['mean_rmsd'] is not None else "N/A"
        print(f"{i:<4} {r['mutant']:<8} {r['ligand']:<12} {r['rep']:<4} {mean_str:<10} {std_str:<10} {max_str:<10} {r['status']}")

    # Save CSV
    csv_path = analysis_dir / "binding_site_rmsd.csv"
    with open(csv_path, 'w') as f:
        f.write("sim_id,mutant,ligand,rep,mean_rmsd,std_rmsd,max_rmsd,status\n")
        for i, r in enumerate(results, 1):
            mean = r['mean_rmsd'] if r['mean_rmsd'] is not None else ""
            std = r.get('std_rmsd', "") if r['mean_rmsd'] is not None else ""
            max_r = r.get('max_rmsd', "") if r['mean_rmsd'] is not None else ""
            f.write(f"{i},{r['mutant']},{r['ligand']},{r['rep']},{mean},{std},{max_r},{r['status']}\n")

    print(f"\nResults saved to: {csv_path}")
    print("\nBinding site RMSD analysis complete!")

if __name__ == "__main__":
    main()
