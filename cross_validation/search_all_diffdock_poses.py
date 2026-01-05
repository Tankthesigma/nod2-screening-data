#!/usr/bin/env python3
"""Search all 100 DiffDock poses for any that agree with GNINA pocket."""

import numpy as np
import os
import re

DIFFDOCK_LRR_DIR = r"C:\Users\vasud\nod2-screening-data\cross_validation\diffdock_lrr"
POCKET_CENTER = np.array([-38.10, -15.89, 22.69])  # From analysis


def load_sdf_coords(sdf_file):
    """Load ligand heavy atom coordinates from SDF file."""
    coords = []
    with open(sdf_file, 'r') as f:
        lines = f.readlines()

    counts_idx = None
    for i, line in enumerate(lines):
        if 'V2000' in line or 'V3000' in line:
            counts_idx = i
            break

    if counts_idx is None:
        return None

    try:
        n_atoms = int(lines[counts_idx][:3].strip())
    except ValueError:
        return None

    for i in range(counts_idx + 1, counts_idx + 1 + n_atoms):
        if i >= len(lines):
            break
        line = lines[i]
        try:
            x = float(line[0:10].strip())
            y = float(line[10:20].strip())
            z = float(line[20:30].strip())
            atom_symbol = line[31:34].strip()
            if atom_symbol != 'H':
                coords.append([x, y, z])
        except (ValueError, IndexError):
            continue

    return np.array(coords) if coords else None


def main():
    print("Searching all 100 DiffDock poses for pocket agreement...")
    print(f"Pocket center: {POCKET_CENTER}")
    print(f"GNINA distance to pocket: 4.64 A")
    print()

    results = []

    for rank in range(1, 101):
        sdf_files = [f for f in os.listdir(DIFFDOCK_LRR_DIR)
                     if f.startswith(f'rank{rank}_') and f.endswith('.sdf')]
        if sdf_files:
            sdf_path = os.path.join(DIFFDOCK_LRR_DIR, sdf_files[0])
            coords = load_sdf_coords(sdf_path)
            if coords is not None:
                centroid = np.mean(coords, axis=0)
                dist = np.linalg.norm(centroid - POCKET_CENTER)
                conf = sdf_files[0].split('confidence')[1].replace('.sdf', '')
                results.append((rank, dist, conf))

    # Sort by distance
    results.sort(key=lambda x: x[1])

    print("Top 20 closest poses to GNINA pocket:")
    print("-" * 50)
    for rank, dist, conf in results[:20]:
        marker = " <-- CLOSE!" if dist < 10 else ""
        print(f"  Rank {rank:3d} (conf={conf:8s}): {dist:.2f} A{marker}")

    print()
    print("Summary:")
    print(f"  Closest pose: Rank {results[0][0]} at {results[0][1]:.2f} A")
    print(f"  Poses within 10A of pocket: {sum(1 for r in results if r[1] < 10)}")
    print(f"  Poses within 20A of pocket: {sum(1 for r in results if r[1] < 20)}")
    print(f"  Poses within 30A of pocket: {sum(1 for r in results if r[1] < 30)}")

    if results[0][1] > 20:
        print("\n  CONCLUSION: DiffDock found a DIFFERENT binding site than GNINA!")
    elif results[0][1] < 10:
        print(f"\n  CONCLUSION: DiffDock Rank {results[0][0]} agrees with GNINA!")


if __name__ == "__main__":
    main()
