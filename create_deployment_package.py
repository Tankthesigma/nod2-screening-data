#!/usr/bin/env python
"""Create deployment package for Vast.ai solvent FEP simulations."""
import zipfile
import os
from pathlib import Path

BASE = Path("/mnt/c/Users/vasud/nod2-screening-data")
OUTPUT = BASE / "solvent_fep_deployment.zip"

# Files to include
FILES_TO_INCLUDE = [
    # Runner scripts
    ("fep_pmx/run_fep_gpu.py", "fep_pmx/run_fep_gpu.py"),
    ("fep_pmx_natural/run_fep_gpu.py", "fep_pmx_natural/run_fep_gpu.py"),
    ("deploy_solvent_fep.py", "deploy_solvent_fep.py"),

    # Febuxostat solvent (only system files, no u_nk.npy)
    ("fep_pmx/solvent/alchemical_system.xml", "fep_pmx/solvent/alchemical_system.xml"),
    ("fep_pmx/solvent/topology.pdb", "fep_pmx/solvent/topology.pdb"),
    ("fep_pmx/solvent/positions.npy", "fep_pmx/solvent/positions.npy"),
    ("fep_pmx/solvent/lambda_schedule.npy", "fep_pmx/solvent/lambda_schedule.npy"),

    # CID_10120 (Bufadienolide) solvent - NOTE: Previously mislabeled as CID_10592
    ("fep_pmx_natural/solvent/alchemical_system.xml", "fep_pmx_natural/solvent/alchemical_system.xml"),
    ("fep_pmx_natural/solvent/topology.pdb", "fep_pmx_natural/solvent/topology.pdb"),
    ("fep_pmx_natural/solvent/positions.npy", "fep_pmx_natural/solvent/positions.npy"),
    ("fep_pmx_natural/solvent/lambda_schedule.npy", "fep_pmx_natural/solvent/lambda_schedule.npy"),
    ("fep_pmx_natural/solvent/ligand_indices.npy", "fep_pmx_natural/solvent/ligand_indices.npy"),
]


def main():
    print("Creating deployment package...")
    print(f"Output: {OUTPUT}")

    with zipfile.ZipFile(OUTPUT, 'w', zipfile.ZIP_DEFLATED) as zf:
        for src, dst in FILES_TO_INCLUDE:
            src_path = BASE / src
            if src_path.exists():
                print(f"  + {dst}")
                zf.write(src_path, dst)
            else:
                print(f"  ! MISSING: {src}")

    # Check size
    size_mb = OUTPUT.stat().st_size / (1024 * 1024)
    print(f"\nPackage size: {size_mb:.1f} MB")
    print(f"\nReady for upload: {OUTPUT}")


if __name__ == "__main__":
    main()
