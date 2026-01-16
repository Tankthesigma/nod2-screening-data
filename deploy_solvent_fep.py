#!/usr/bin/env python
"""
Deploy and run solvent leg FEP simulations on Vast.ai GPU instance.

Windows to run:
- Febuxostat solvent: 17, 18 (2 windows)
- CID_10120 (Bufadienolide) solvent: 0-19 (20 windows)
  NOTE: Previously mislabeled as CID_10592/Dihydrocortisol - correct ID is CID_10120

Total: 22 windows
"""
import subprocess
import sys
import os
from pathlib import Path
import time

# Configuration
N_GPUS = 8  # 8x 4070S Ti
BASE_DIR = Path("/workspace/fep_data")

# Window assignments
JOBS = [
    # (compound, window_idx, gpu_id)
    # Febuxostat solvent (2 windows)
    ("fep_pmx", "solvent", 17, 0),
    ("fep_pmx", "solvent", 18, 1),
    # CID_10120 solvent (20 windows) - distribute across GPUs
    ("fep_pmx_natural", "solvent", 0, 2),
    ("fep_pmx_natural", "solvent", 1, 3),
    ("fep_pmx_natural", "solvent", 2, 4),
    ("fep_pmx_natural", "solvent", 3, 5),
    ("fep_pmx_natural", "solvent", 4, 6),
    ("fep_pmx_natural", "solvent", 5, 7),
    ("fep_pmx_natural", "solvent", 6, 0),  # Reuse GPU 0 after feb completes
    ("fep_pmx_natural", "solvent", 7, 1),  # Reuse GPU 1 after feb completes
    ("fep_pmx_natural", "solvent", 8, 2),
    ("fep_pmx_natural", "solvent", 9, 3),
    ("fep_pmx_natural", "solvent", 10, 4),
    ("fep_pmx_natural", "solvent", 11, 5),
    ("fep_pmx_natural", "solvent", 12, 6),
    ("fep_pmx_natural", "solvent", 13, 7),
    ("fep_pmx_natural", "solvent", 14, 0),
    ("fep_pmx_natural", "solvent", 15, 1),
    ("fep_pmx_natural", "solvent", 16, 2),
    ("fep_pmx_natural", "solvent", 17, 3),
    ("fep_pmx_natural", "solvent", 18, 4),
    ("fep_pmx_natural", "solvent", 19, 5),
]


def run_window(compound, system, window_idx, gpu_id):
    """Run a single FEP window."""
    cmd = [
        "python", "run_fep_gpu.py",
        str(BASE_DIR / compound),
        system,
        str(window_idx),
        str(window_idx),
        str(gpu_id)
    ]
    print(f"[GPU {gpu_id}] {compound}/{system} window {window_idx}")
    return subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)


def run_sequential_on_gpu(jobs_for_gpu):
    """Run jobs sequentially on a single GPU."""
    for compound, system, window_idx, gpu_id in jobs_for_gpu:
        print(f"\n{'='*60}")
        print(f"[GPU {gpu_id}] Starting {compound}/{system} window {window_idx}")
        print(f"{'='*60}")

        result = subprocess.run(
            ["python", "run_fep_gpu.py",
             str(BASE_DIR / compound), system,
             str(window_idx), str(window_idx), str(gpu_id)],
            capture_output=True, text=True
        )

        if result.returncode != 0:
            print(f"ERROR on {compound}/{system} window {window_idx}:")
            print(result.stderr[-1000:] if result.stderr else "No stderr")
        else:
            print(f"[GPU {gpu_id}] Completed window {window_idx}")


def main():
    """Run all solvent FEP windows."""
    print("=" * 70)
    print("SOLVENT LEG FEP DEPLOYMENT")
    print("=" * 70)
    print(f"\nTotal windows to run: {len(JOBS)}")
    print(f"  - Febuxostat: 2 (windows 17-18)")
    print(f"  - CID_10120: 20 (windows 0-19)")
    print(f"\nUsing {N_GPUS} GPUs")

    # Group jobs by GPU for sequential execution
    gpu_jobs = {i: [] for i in range(N_GPUS)}
    for job in JOBS:
        compound, system, window_idx, gpu_id = job
        gpu_jobs[gpu_id].append(job)

    print("\nGPU assignments:")
    for gpu_id, jobs in gpu_jobs.items():
        windows = [f"{j[0].split('_')[-1]}:{j[2]}" for j in jobs]
        print(f"  GPU {gpu_id}: {windows}")

    # Run all GPUs in parallel, each running its jobs sequentially
    import multiprocessing
    processes = []

    for gpu_id in range(N_GPUS):
        if gpu_jobs[gpu_id]:
            p = multiprocessing.Process(
                target=run_sequential_on_gpu,
                args=(gpu_jobs[gpu_id],)
            )
            p.start()
            processes.append(p)

    # Wait for all to complete
    for p in processes:
        p.join()

    print("\n" + "=" * 70)
    print("ALL SOLVENT FEP SIMULATIONS COMPLETE")
    print("=" * 70)

    # Check for completion
    print("\nVerifying results...")
    missing = []

    # Febuxostat
    for w in [17, 18]:
        unk = BASE_DIR / "fep_pmx" / "solvent" / f"window_{w:02d}" / "u_nk.npy"
        if not unk.exists():
            missing.append(f"fep_pmx/solvent/window_{w:02d}")

    # Natural
    for w in range(20):
        unk = BASE_DIR / "fep_pmx_natural" / "solvent" / f"window_{w:02d}" / "u_nk.npy"
        if not unk.exists():
            missing.append(f"fep_pmx_natural/solvent/window_{w:02d}")

    if missing:
        print(f"\nMISSING u_nk.npy files: {len(missing)}")
        for m in missing[:5]:
            print(f"  - {m}")
    else:
        print("\n[OK] All 22 windows completed successfully!")

    return 0 if not missing else 1


if __name__ == "__main__":
    sys.exit(main())
