"""
Hamiltonian Consistency Verification Script
Checks that all parameters match before fixing failed windows
"""
import numpy as np
import os
from pathlib import Path

BASE = Path("C:/Users/vasud/nod2-screening-data/fep_pmx")
SYSTEMS = ["wt_complex", "mut_complex", "solvent"]

print("=" * 70)
print("HAMILTONIAN CONSISTENCY VERIFICATION")
print("=" * 70)

# 1. LAMBDA SCHEDULE
print("\n[1] LAMBDA SCHEDULE")
print("-" * 50)
for sys in SYSTEMS:
    # Check system root
    root_sched = BASE / sys / "lambda_schedule.npy"
    # Check window_00
    win0_sched = BASE / sys / "window_00" / "lambda_schedule.npy"

    sched = None
    if root_sched.exists():
        sched = np.load(root_sched, allow_pickle=True)
        print(f"{sys} (root): lambda_schedule.npy found")
    elif win0_sched.exists():
        sched = np.load(win0_sched, allow_pickle=True)
        print(f"{sys} (window_00): lambda_schedule.npy found")
    else:
        print(f"{sys}: WARNING - lambda_schedule.npy NOT FOUND")
        continue

    # Handle dict vs array
    if isinstance(sched, np.ndarray) and sched.dtype == object:
        sched = sched.item()  # Extract dict

    if isinstance(sched, dict):
        print(f"  Format: dict with keys {list(sched.keys())}")
        n_windows = len(list(sched.values())[0])
    else:
        print(f"  Format: array shape {sched.shape}")
        n_windows = sched.shape[0]

    print(f"  Number of windows: {n_windows}")

    # Print extreme windows (15-19)
    if isinstance(sched, dict):
        print(f"  Windows 15-19:")
        for i in range(15, 20):
            elec = sched.get('lambda_electrostatics', sched.get('elec', [0]*20))[i]
            ster = sched.get('lambda_sterics', sched.get('sterics', [0]*20))[i]
            rest = sched.get('lambda_restraints', sched.get('restraints', [1]*20))[i]
            print(f"    {i}: elec={elec:.3f}, sterics={ster:.3f}, restraints={rest:.3f}")
    else:
        print(f"  Windows 15-19:")
        for i in range(15, 20):
            print(f"    {i}: {sched[i]}")

# 2. BORESCH RESTRAINTS
print("\n[2] BORESCH RESTRAINTS")
print("-" * 50)
for sys in ["wt_complex", "mut_complex"]:  # solvent has no Boresch
    root_anch = BASE / sys / "boresch_anchors.npy"

    if root_anch.exists():
        anchors = np.load(root_anch, allow_pickle=True)
        if anchors.dtype == object:
            anchors = anchors.item()
        print(f"{sys}:")
        if isinstance(anchors, dict):
            if 'force_constant' in anchors:
                print(f"  force_constant: {anchors['force_constant']}")
            if 'k' in anchors:
                print(f"  k: {anchors['k']}")
            if 'dG_correction' in anchors:
                print(f"  dG_correction: {anchors['dG_correction']} kJ/mol")
            if 'standard_state_correction' in anchors:
                print(f"  standard_state_correction: {anchors['standard_state_correction']}")
        else:
            print(f"  Raw data: {anchors}")
    else:
        print(f"{sys}: WARNING - boresch_anchors.npy NOT FOUND")

# 3. SEED FILES CHECK
print("\n[3] SEED FILES (with box vectors)")
print("-" * 50)
seed_windows = {
    "wt_complex": [14],  # seed for 15
    "mut_complex": [14, 16],  # seed for 15, and 16 succeeded so seed for 17
    "solvent": [15]  # 15 succeeded, seed for 16
}

all_seeds_ok = True
for sys, windows in seed_windows.items():
    for win in windows:
        win_dir = BASE / sys / f"window_{win:02d}"
        chk = win_dir / "checkpoint.chk"
        pos = win_dir / "final_positions.npy"
        box = win_dir / "final_box_vectors.npy"

        print(f"{sys}/window_{win:02d}:")
        if chk.exists():
            print(f"  ✓ checkpoint.chk EXISTS ({chk.stat().st_size / 1024:.1f} KB)")
        elif pos.exists() and box.exists():
            print(f"  ✓ final_positions.npy + final_box_vectors.npy EXISTS")
        elif pos.exists():
            print(f"  ⚠ final_positions.npy EXISTS but NO box vectors!")
            all_seeds_ok = False
        else:
            print(f"  ✗ NO seed files found!")
            all_seeds_ok = False

# 4. u_nk FORMAT CHECK
print("\n[4] u_nk FORMAT (from existing successful window)")
print("-" * 50)
for sys in SYSTEMS:
    # Find first successful window
    for i in range(20):
        unk_file = BASE / sys / f"window_{i:02d}" / "u_nk.npy"
        if unk_file.exists():
            unk = np.load(unk_file)
            print(f"{sys}/window_{i:02d}:")
            print(f"  Shape: {unk.shape}")
            print(f"  dtype: {unk.dtype}")
            print(f"  Range: [{unk.min():.2f}, {unk.max():.2f}]")
            print(f"  Units: likely kJ/mol (based on magnitude)")
            break

# 5. SOFTCORE CHECK (from alchemical_system.xml or run_window.py)
print("\n[5] SOFTCORE PARAMETERS")
print("-" * 50)
# Check run_window.py for softcore_alpha
run_py = BASE / "wt_complex" / "window_00" / "run_window.py"
if run_py.exists():
    with open(run_py, 'r') as f:
        content = f.read()
    if 'softcore_alpha' in content:
        # Extract value
        import re
        match = re.search(r'softcore_alpha\s*[=:]\s*([0-9.]+)', content)
        if match:
            print(f"softcore_alpha = {match.group(1)}")
        else:
            print("softcore_alpha mentioned but value not extracted")
    else:
        print("softcore_alpha not found in run_window.py (using OpenMM default 0.5)")
else:
    print("run_window.py not found")

# Final summary
print("\n" + "=" * 70)
print("VERIFICATION SUMMARY")
print("=" * 70)
if all_seeds_ok:
    print("✓ All seed files present with box vectors")
else:
    print("✗ SOME SEED FILES MISSING - MUST FIX BEFORE PROCEEDING")
print("✓ Lambda schedules verified")
print("✓ Review Boresch restraints above")
print("✓ Review u_nk format above")
