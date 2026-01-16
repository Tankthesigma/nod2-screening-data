#!/usr/bin/env python
"""
Launch all 60 FEP windows for Natural Product CID_10120.
Run this after system setup is complete.
"""
import os
import subprocess
import sys
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed

BASE_DIR = Path("C:\Users\vasud\nod2-screening-data\fep_pmx_natural")
SYSTEMS = ['wt_complex', 'mut_complex', 'solvent']
N_WINDOWS = 20

def run_window(sys_name, window_idx):
    """Run a single FEP window."""
    window_dir = BASE_DIR / sys_name / f"window_{window_idx:02d}"
    script = window_dir / "run_window.py"
    log_file = window_dir / "run.log"

    print(f"[START] {sys_name}/window_{window_idx:02d}")

    with open(log_file, 'w') as log:
        result = subprocess.run(
            [sys.executable, str(script)],
            cwd=str(window_dir),
            stdout=log,
            stderr=subprocess.STDOUT,
        )

    if result.returncode == 0:
        print(f"[DONE] {sys_name}/window_{window_idx:02d}")
    else:
        print(f"[FAIL] {sys_name}/window_{window_idx:02d}")

    return (sys_name, window_idx, result.returncode)

def main():
    print("="*70)
    print("LAUNCHING ALL 60 FEP WINDOWS")
    print("="*70)

    # Build list of all windows
    windows = []
    for sys_name in SYSTEMS:
        for i in range(N_WINDOWS):
            windows.append((sys_name, i))

    print(f"Total windows to run: {len(windows)}")

    # Run sequentially (can be parallelized with GPU resources)
    results = []
    for sys_name, window_idx in windows:
        result = run_window(sys_name, window_idx)
        results.append(result)

    # Summary
    print()
    print("="*70)
    print("SUMMARY")
    print("="*70)

    success = sum(1 for r in results if r[2] == 0)
    failed = sum(1 for r in results if r[2] != 0)

    print(f"Success: {success}/{len(results)}")
    print(f"Failed: {failed}/{len(results)}")

    if failed > 0:
        print("\nFailed windows:")
        for sys_name, window_idx, rc in results:
            if rc != 0:
                print(f"  {sys_name}/window_{window_idx:02d}")

if __name__ == "__main__":
    main()
