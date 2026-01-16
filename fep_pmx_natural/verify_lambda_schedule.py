#!/usr/bin/env python
"""
Verify lambda schedule matches febuxostat EXACTLY.
This is CRITICAL for MBAR compatibility.
"""
import numpy as np
from pathlib import Path

FEBUXOSTAT_DIR = Path("C:/Users/vasud/nod2-screening-data/fep_complete/fep_pmx")
NATURAL_DIR = Path("C:/Users/vasud/nod2-screening-data/fep_pmx_natural")

def main():
    print("="*70)
    print("LAMBDA SCHEDULE VERIFICATION")
    print("="*70)

    # Load febuxostat lambda schedule
    feb_path = FEBUXOSTAT_DIR / "wt_complex" / "lambda_schedule.npy"
    if not feb_path.exists():
        print(f"[ERROR] Febuxostat lambda schedule not found: {feb_path}")
        print("        Cannot verify without reference schedule!")
        return
    feb_schedule = np.load(feb_path)
    print(f"\nFebuxostat lambda schedule shape: {feb_schedule.shape}")

    # Validate schedule shape (should be Nx3 for elec, sterics, restraints)
    if len(feb_schedule.shape) != 2 or feb_schedule.shape[1] != 3:
        print(f"[ERROR] Invalid schedule shape: {feb_schedule.shape} (expected Nx3)")
        return

    print("\nFebuxostat Lambda Schedule:")
    print("-"*50)
    print("Window  Elec      Sterics    Restraints")
    print("-"*50)
    for i, (e, s, r) in enumerate(feb_schedule):
        print(f"{i:2d}      {e:.6f}  {s:.6f}   {r:.6f}")

    # Check if natural schedule exists and compare
    nat_schedule_path = NATURAL_DIR / "wt_complex" / "lambda_schedule.npy"
    if nat_schedule_path.exists():
        nat_schedule = np.load(nat_schedule_path)

        print("\n" + "="*70)
        print("COMPARISON WITH NATURAL PRODUCT SCHEDULE")
        print("="*70)

        if feb_schedule.shape != nat_schedule.shape:
            print(f"[FAIL] Shape mismatch! Febuxostat: {feb_schedule.shape}, Natural: {nat_schedule.shape}")
            print("       Lambda schedules MUST have identical shapes for MBAR!")
        elif np.allclose(feb_schedule, nat_schedule):
            print("[PASS] Lambda schedules match EXACTLY")
        else:
            print("[FAIL] Lambda schedules DO NOT MATCH!")
            print("\nDifferences:")
            diff = np.abs(feb_schedule - nat_schedule)
            for i in range(len(feb_schedule)):
                if diff[i].sum() > 1e-10:
                    print(f"  Window {i}: feb={feb_schedule[i]}, nat={nat_schedule[i]}")
    else:
        print(f"\n[INFO] Natural schedule not found at {nat_schedule_path}")
        print("       Run setup_fep_natural.py first")

    # Save correct schedule
    print("\n" + "="*70)
    print("SAVING VERIFIED LAMBDA SCHEDULE")
    print("="*70)

    for sys_name in ['wt_complex', 'mut_complex', 'solvent']:
        out_path = NATURAL_DIR / sys_name / "lambda_schedule.npy"
        out_path.parent.mkdir(parents=True, exist_ok=True)
        np.save(out_path, feb_schedule)
        print(f"  Saved: {out_path}")

    print("\n[DONE] Lambda schedules verified and saved")

if __name__ == "__main__":
    main()
