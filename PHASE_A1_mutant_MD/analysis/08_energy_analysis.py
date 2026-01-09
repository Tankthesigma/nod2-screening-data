#!/usr/bin/env python3
"""
PHASE A1 ANALYSIS - Energy Analysis (FIXED)
Parses log files for energy and temperature stability.
Fixed log file paths to match actual project structure.
"""

import os
import sys
from pathlib import Path
import platform

if 'microsoft' in platform.uname().release.lower() or os.path.exists('/mnt/c'):
    BASE_DIR = Path("/mnt/c/Users/vasud/nod2-screening-data/PHASE_A1_mutant_MD")
else:
    BASE_DIR = Path(r"C:\Users\vasud\nod2-screening-data\PHASE_A1_mutant_MD")

try:
    import numpy as np
except ImportError:
    print("ERROR: numpy not installed!")
    sys.exit(1)

# Log file locations - FIXED to match actual structure
# Local logs are in BASE_DIR/logs, Vast logs are in same folder as DCDs
def get_log_path(mutant, ligand, rep, folder_key):
    """Get log file path based on folder location."""
    base_name = f"{mutant}_{ligand}_rep{rep}"

    if folder_key == 'trajectories':
        # Local simulations - logs in BASE_DIR/logs
        return BASE_DIR / "logs" / f"{base_name}.log"
    else:
        # Vast downloads - logs in same folder as DCDs
        folder_map = {
            '5080': BASE_DIR / "vast_downloads" / "5080",
            '5090': BASE_DIR / "vast_downloads" / "5090",
            'new_5090': BASE_DIR / "vast_downloads" / "new_5090",
        }
        return folder_map[folder_key] / f"{base_name}.log"

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

def parse_log_file(log_path):
    """
    Parse OpenMM log file.
    Expected format: Progress%,Step,Time(ps),PotentialEnergy,KineticEnergy,Temperature,Speed,ETA
    """
    if not log_path.exists():
        return None

    potential_energies = []
    kinetic_energies = []
    temperatures = []
    times = []

    with open(log_path, 'r') as f:
        for line in f:
            # Skip header
            if line.startswith('#') or 'Progress' in line or line.strip() == '':
                continue

            parts = line.strip().split(',')
            if len(parts) >= 6:
                try:
                    # Format: progress,step,time,potential,kinetic,temperature,speed,eta
                    time_ps = float(parts[2])
                    potential = float(parts[3])
                    kinetic = float(parts[4])
                    temp = float(parts[5])

                    times.append(time_ps)
                    potential_energies.append(potential)
                    kinetic_energies.append(kinetic)
                    temperatures.append(temp)
                except (ValueError, IndexError):
                    continue

    if not temperatures:
        return None

    return {
        'times': np.array(times),
        'potential': np.array(potential_energies),
        'kinetic': np.array(kinetic_energies),
        'temperature': np.array(temperatures),
    }

def main():
    print("=" * 80)
    print("PHASE A1 ANALYSIS - ENERGY & TEMPERATURE")
    print("=" * 80)
    print()

    analysis_dir = BASE_DIR / "analysis"
    analysis_dir.mkdir(exist_ok=True)

    results = []

    for i, (mutant, ligand, rep, folder_key) in enumerate(SIMULATIONS, 1):
        print(f"[{i}/12] {mutant}_{ligand}_rep{rep}")
        log_path = get_log_path(mutant, ligand, rep, folder_key)

        if not log_path.exists():
            print(f"  Log file not found: {log_path}")
            results.append({
                'mutant': mutant, 'ligand': ligand, 'rep': rep,
                'temp_mean': None, 'status': 'LOG_NOT_FOUND'
            })
            continue

        data = parse_log_file(log_path)

        if data is None:
            print(f"  Could not parse log file")
            results.append({
                'mutant': mutant, 'ligand': ligand, 'rep': rep,
                'temp_mean': None, 'status': 'PARSE_ERROR'
            })
            continue

        # Temperature analysis
        temp_mean = data['temperature'].mean()
        temp_std = data['temperature'].std()
        temp_min = data['temperature'].min()
        temp_max = data['temperature'].max()

        # Potential energy analysis
        pe_mean = data['potential'].mean()
        pe_std = data['potential'].std()

        # Total energy
        total_energy = data['potential'] + data['kinetic']
        te_mean = total_energy.mean()
        te_std = total_energy.std()

        # Stability check
        if temp_std < 5.0 and 300 < temp_mean < 320:
            temp_status = "STABLE"
        elif temp_std < 10.0:
            temp_status = "ACCEPTABLE"
        else:
            temp_status = "UNSTABLE"

        print(f"  Temperature: {temp_mean:.1f} +/- {temp_std:.1f} K (range: {temp_min:.1f}-{temp_max:.1f} K) - {temp_status}")
        print(f"  Potential Energy: {pe_mean/1000:.1f} +/- {pe_std/1000:.1f} MJ/mol")

        results.append({
            'mutant': mutant, 'ligand': ligand, 'rep': rep,
            'temp_mean': temp_mean, 'temp_std': temp_std,
            'temp_min': temp_min, 'temp_max': temp_max,
            'pe_mean': pe_mean, 'pe_std': pe_std,
            'te_mean': te_mean, 'te_std': te_std,
            'n_points': len(data['temperature']),
            'status': temp_status
        })

        print()

    # Print summary
    print("=" * 80)
    print("ENERGY ANALYSIS SUMMARY")
    print("=" * 80)
    print(f"{'ID':<4} {'Mutant':<8} {'Ligand':<12} {'Rep':<4} {'Temp (K)':<15} {'PE (MJ/mol)':<15} {'Status'}")
    print("-" * 80)

    for i, r in enumerate(results, 1):
        if r['temp_mean'] is not None:
            temp_str = f"{r['temp_mean']:.1f}+/-{r['temp_std']:.1f}"
            pe_str = f"{r['pe_mean']/1000:.1f}+/-{r['pe_std']/1000:.1f}"
        else:
            temp_str = "N/A"
            pe_str = "N/A"
        print(f"{i:<4} {r['mutant']:<8} {r['ligand']:<12} {r['rep']:<4} {temp_str:<15} {pe_str:<15} {r['status']}")

    # Save CSV
    csv_path = analysis_dir / "energy_analysis.csv"
    with open(csv_path, 'w') as f:
        f.write("sim_id,mutant,ligand,rep,temp_mean,temp_std,temp_min,temp_max,pe_mean,pe_std,te_mean,te_std,n_points,status\n")
        for i, r in enumerate(results, 1):
            temp_mean = r['temp_mean'] if r['temp_mean'] is not None else ""
            temp_std = r.get('temp_std', "") if r['temp_mean'] is not None else ""
            temp_min = r.get('temp_min', "") if r['temp_mean'] is not None else ""
            temp_max = r.get('temp_max', "") if r['temp_mean'] is not None else ""
            pe_mean = r.get('pe_mean', "") if r['temp_mean'] is not None else ""
            pe_std = r.get('pe_std', "") if r['temp_mean'] is not None else ""
            te_mean = r.get('te_mean', "") if r['temp_mean'] is not None else ""
            te_std = r.get('te_std', "") if r['temp_mean'] is not None else ""
            n_points = r.get('n_points', "") if r['temp_mean'] is not None else ""
            f.write(f"{i},{r['mutant']},{r['ligand']},{r['rep']},{temp_mean},{temp_std},{temp_min},{temp_max},{pe_mean},{pe_std},{te_mean},{te_std},{n_points},{r['status']}\n")

    print(f"\nResults saved to: {csv_path}")
    print("\nEnergy analysis complete!")

if __name__ == "__main__":
    main()
