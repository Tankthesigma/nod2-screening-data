#!/usr/bin/env python3
"""
PHASE A1 ANALYSIS - Master Runner
Runs all analysis scripts in sequence.

Usage: python run_all_analysis.py [--skip-discovery]
"""

import subprocess
import sys
import os
from pathlib import Path
import platform

if 'microsoft' in platform.uname().release.lower() or os.path.exists('/mnt/c'):
    BASE_DIR = Path("/mnt/c/Users/vasud/nod2-screening-data/PHASE_A1_mutant_MD")
else:
    BASE_DIR = Path(r"C:\Users\vasud\nod2-screening-data\PHASE_A1_mutant_MD")

ANALYSIS_DIR = BASE_DIR / "analysis"

# Scripts to run in order
SCRIPTS = [
    ("01_file_discovery.py", "File Discovery"),
    ("02_pocket_occupancy.py", "Pocket Occupancy"),
    ("03_ligand_rmsd.py", "Ligand RMSD"),
    ("04_contact_analysis.py", "Contact Analysis"),
    ("05_hbond_analysis.py", "H-bond Analysis"),
    ("06_binding_site_rmsd.py", "Binding Site RMSD"),
    ("07_ligand_com_distance.py", "Ligand COM Distance"),
    ("08_energy_analysis.py", "Energy Analysis"),
    ("09_generate_report.py", "Report Generation"),
]

def run_script(script_name, description):
    """Run a single analysis script."""
    script_path = ANALYSIS_DIR / script_name

    if not script_path.exists():
        print(f"  ERROR: Script not found: {script_path}")
        return False

    print(f"\n{'='*60}")
    print(f"Running: {description}")
    print(f"Script: {script_name}")
    print('='*60)

    try:
        result = subprocess.run(
            [sys.executable, str(script_path)],
            capture_output=False,
            text=True,
            cwd=str(BASE_DIR)
        )
        return result.returncode == 0
    except Exception as e:
        print(f"  ERROR: {e}")
        return False

def main():
    print("=" * 70)
    print("PHASE A1 ANALYSIS - MASTER RUNNER")
    print("=" * 70)
    print()
    print(f"Base directory: {BASE_DIR}")
    print(f"Analysis directory: {ANALYSIS_DIR}")
    print()

    # Check for skip flag
    skip_discovery = "--skip-discovery" in sys.argv

    # Create plots directory
    plots_dir = ANALYSIS_DIR / "plots"
    plots_dir.mkdir(parents=True, exist_ok=True)

    results = []

    for script_name, description in SCRIPTS:
        if skip_discovery and script_name == "01_file_discovery.py":
            print(f"\nSkipping: {description}")
            results.append((description, "SKIPPED"))
            continue

        success = run_script(script_name, description)
        results.append((description, "OK" if success else "FAILED"))

        if not success:
            print(f"\nWARNING: {description} failed!")
            # Continue anyway

    # Print summary
    print("\n" + "=" * 70)
    print("ANALYSIS SUMMARY")
    print("=" * 70)

    for desc, status in results:
        status_symbol = "[OK]" if status == "OK" else "[FAIL]" if status == "FAILED" else "[--]"
        print(f"  {status_symbol} {desc}: {status}")

    print()
    print("Output files in:", ANALYSIS_DIR)
    print()

    # List output files
    if ANALYSIS_DIR.exists():
        csv_files = list(ANALYSIS_DIR.glob("*.csv"))
        md_files = list(ANALYSIS_DIR.glob("*.md"))
        txt_files = list(ANALYSIS_DIR.glob("*.txt"))

        print("CSV files:")
        for f in csv_files:
            print(f"  - {f.name}")

        print("\nReport files:")
        for f in md_files + txt_files:
            print(f"  - {f.name}")

    print("\n" + "=" * 70)
    print("ANALYSIS COMPLETE")
    print("=" * 70)

if __name__ == "__main__":
    main()
