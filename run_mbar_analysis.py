#!/usr/bin/env python
"""
MBAR Analysis for NOD2-Crohn's FEP Calculations

Analyzes all 60 FEP windows across 3 systems:
- solvent: ligand decoupling in water
- wt_complex: ligand decoupling from WT NOD2
- mut_complex: ligand decoupling from R702W mutant NOD2

Calculates:
- ΔG_decouple for each system
- ΔG_bind for WT and Mut
- ΔΔG_bind = ΔG_bind(Mut) - ΔG_bind(WT)
"""

import os
import json
import numpy as np

# Try to import pymbar
try:
    import pymbar
    from pymbar import MBAR
    print(f"Using pymbar version: {pymbar.__version__}")
except ImportError:
    print("ERROR: pymbar not installed. Run: pip install pymbar")
    exit(1)

# Constants
KB = 0.008314462618  # kJ/mol/K (Boltzmann constant)
TEMPERATURE = 310.0  # K
BETA = 1.0 / (KB * TEMPERATURE)  # 1/(kJ/mol)
KJ_TO_KCAL = 1.0 / 4.184

# Base path
BASE_PATH = "C:/Users/vasud/nod2-screening-data/fep_complete/fep_pmx"

def load_lambda_schedule(system_dir):
    """Load and print lambda schedule to verify ordering."""
    schedule_path = os.path.join(system_dir, "lambda_schedule.npy")
    schedule = np.load(schedule_path)

    print(f"\n  Lambda schedule shape: {schedule.shape}")
    print(f"  Window 0 (first):  elec={schedule[0,0]:.3f}, sterics={schedule[0,1]:.3f}, restraints={schedule[0,2]:.3f}")
    print(f"  Window 19 (last):  elec={schedule[19,0]:.3f}, sterics={schedule[19,1]:.3f}, restraints={schedule[19,2]:.3f}")

    # Determine direction
    if schedule[0, 1] > schedule[19, 1]:  # sterics decreases
        print("  Direction: λ=1 (coupled) → λ=0 (decoupled)")
        coupled_idx = 0
        decoupled_idx = 19
    else:
        print("  Direction: λ=0 (decoupled) → λ=1 (coupled)")
        coupled_idx = 19
        decoupled_idx = 0

    return schedule, coupled_idx, decoupled_idx


def load_u_nk_data(system_dir, system_name):
    """Load all u_nk.npy files for a system."""
    n_windows = 20
    u_nk_list = []
    N_k = []

    print(f"\n  Loading {system_name} windows:")
    for k in range(n_windows):
        window_dir = os.path.join(system_dir, f"window_{k:02d}")
        u_nk_path = os.path.join(window_dir, "u_nk.npy")

        if not os.path.exists(u_nk_path):
            print(f"    ERROR: Missing {u_nk_path}")
            return None, None

        # Load potential energies (kJ/mol)
        u_nk = np.load(u_nk_path)
        n_samples = u_nk.shape[0]

        # Convert to reduced potentials: u* = beta * U
        u_nk_reduced = u_nk * BETA

        u_nk_list.append(u_nk_reduced)
        N_k.append(n_samples)

        if k in [0, 9, 19]:
            print(f"    Window {k:02d}: {n_samples} samples, shape {u_nk.shape}")

    N_k = np.array(N_k)
    print(f"  Total samples: {N_k.sum()}")
    print(f"  N_k: {N_k}")

    return u_nk_list, N_k


def build_u_kn_matrix(u_nk_list, N_k):
    """
    Build the u_kn matrix for MBAR.

    u_kn[k, n] = reduced potential of sample n evaluated at state k

    For each window i, we have samples that need to be evaluated at all states k.
    The u_nk data already contains energies at all states for each sample.
    """
    K = len(N_k)  # number of states
    N_total = N_k.sum()  # total samples

    # Initialize u_kn matrix
    u_kn = np.zeros((K, N_total))

    # Fill in the matrix
    sample_idx = 0
    for i, u_nk in enumerate(u_nk_list):
        n_i = N_k[i]  # samples from state i
        # u_nk[n, k] has shape (n_samples, K)
        # We need u_kn[k, samples_from_state_i]
        for k in range(K):
            u_kn[k, sample_idx:sample_idx + n_i] = u_nk[:, k]
        sample_idx += n_i

    return u_kn


def compute_overlap_matrix(mbar_obj):
    """Compute and return the overlap matrix."""
    try:
        overlap = mbar_obj.compute_overlap()
        if hasattr(overlap, 'matrix'):
            return overlap.matrix
        else:
            return overlap['matrix']
    except Exception as e:
        print(f"  Warning: Could not compute overlap matrix: {e}")
        return None


def check_overlap(overlap_matrix, threshold=0.03):
    """Check for poor overlap between neighboring states."""
    issues = []
    K = overlap_matrix.shape[0]

    for i in range(K - 1):
        overlap_val = overlap_matrix[i, i + 1]
        if overlap_val < threshold:
            issues.append((i, i + 1, overlap_val))

    return issues


def run_mbar_for_system(system_name, system_dir):
    """Run MBAR analysis for a single system."""
    print(f"\n{'='*60}")
    print(f"MBAR ANALYSIS: {system_name}")
    print('='*60)

    # Load lambda schedule
    schedule, coupled_idx, decoupled_idx = load_lambda_schedule(system_dir)

    # Load u_nk data
    u_nk_list, N_k = load_u_nk_data(system_dir, system_name)
    if u_nk_list is None:
        return None, None, None

    # Build u_kn matrix
    print("\n  Building u_kn matrix...")
    u_kn = build_u_kn_matrix(u_nk_list, N_k)
    print(f"  u_kn shape: {u_kn.shape}")

    # Check for NaN/Inf
    if not np.all(np.isfinite(u_kn)):
        n_nan = np.sum(~np.isfinite(u_kn))
        print(f"  WARNING: {n_nan} non-finite values in u_kn!")

    # Run MBAR
    print("\n  Running MBAR...")
    try:
        mbar = MBAR(u_kn, N_k, verbose=False)
        print("  MBAR converged successfully!")
    except Exception as e:
        print(f"  ERROR: MBAR failed: {e}")
        return None, None, None

    # Get free energy differences
    try:
        results = mbar.compute_free_energy_differences()
        if hasattr(results, 'df'):
            delta_f = results.df
            delta_f_err = results.ddf
        else:
            delta_f = results['Delta_f']
            delta_f_err = results['dDelta_f']
    except Exception as e:
        print(f"  ERROR getting free energies: {e}")
        return None, None, None

    # ΔG_decouple = G(decoupled) - G(coupled) = f[decoupled_idx] - f[coupled_idx]
    # In reduced units (kT)
    dG_reduced = delta_f[coupled_idx, decoupled_idx]
    dG_err_reduced = delta_f_err[coupled_idx, decoupled_idx]

    # Convert to kJ/mol: ΔG = dG_reduced / beta
    dG_kj = dG_reduced / BETA
    dG_err_kj = dG_err_reduced / BETA

    # Convert to kcal/mol
    dG_kcal = dG_kj * KJ_TO_KCAL
    dG_err_kcal = dG_err_kj * KJ_TO_KCAL

    print(f"\n  ΔG_decouple (coupled → decoupled):")
    print(f"    {dG_kj:.2f} ± {dG_err_kj:.2f} kJ/mol")
    print(f"    {dG_kcal:.2f} ± {dG_err_kcal:.2f} kcal/mol")

    # Compute overlap matrix
    print("\n  Computing overlap matrix...")
    overlap = compute_overlap_matrix(mbar)

    if overlap is not None:
        # Check for poor overlap
        issues = check_overlap(overlap)
        if issues:
            print("  WARNING: Poor overlap detected:")
            for i, j, val in issues:
                print(f"    States {i}↔{j}: {val:.4f}")
        else:
            print("  Overlap OK (all neighboring states > 0.03)")

        # Print condensed overlap for key transitions
        print("\n  Key overlaps (diagonal neighbors):")
        for i in [0, 5, 10, 14, 15, 16, 17, 18]:
            if i < 19:
                print(f"    {i:2d}↔{i+1:2d}: {overlap[i, i+1]:.4f}")

    return dG_kcal, dG_err_kcal, overlap


def load_boresch_correction(system_dir, system_name):
    """
    Load Boresch restraint correction.

    The analytical correction for Boresch restraints is:
    ΔG_Boresch = kT * ln(V° / V_site)  (POSITIVE)

    where V_site is the configurational volume of the restrained ligand
    and V° is the standard state volume (1660 Å³ = 1 L/mol).

    Since V° >> V_site, this correction is typically +2 to +8 kcal/mol.
    """
    # Try JSON first
    json_path = os.path.join(system_dir, "boresch_restraints.json")
    if os.path.exists(json_path):
        with open(json_path, 'r') as f:
            data = json.load(f)
            if 'dG_correction_kcal' in data:
                return data['dG_correction_kcal']
            elif 'dG_correction_kj' in data:
                return data['dG_correction_kj'] * KJ_TO_KCAL

    # Try anchors file - this contains the pre-computed correction
    anchors_path = os.path.join(system_dir, "boresch_anchors.npy")
    if os.path.exists(anchors_path):
        try:
            anchors = np.load(anchors_path, allow_pickle=True).item()
            if 'standard_state_correction_kcal_mol' in anchors:
                correction = anchors['standard_state_correction_kcal_mol']
                print(f"  Loaded from boresch_anchors.npy: {correction:.3f} kcal/mol")
                return correction
        except Exception as e:
            print(f"  Warning: Could not read boresch_anchors.npy: {e}")

    # Fallback estimate (should not be reached if setup was done correctly)
    # Note: Boresch correction is POSITIVE (releasing restraints to standard state)
    print(f"  WARNING: Using estimated Boresch correction for {system_name}")
    return 7.6  # kcal/mol (typical positive value)


def main():
    print("="*60)
    print("NOD2-Crohn's FEP MBAR Analysis")
    print("="*60)
    print(f"\nTemperature: {TEMPERATURE} K")
    print(f"Beta: {BETA:.6f} mol/kJ")
    print(f"kB: {KB} kJ/mol/K")

    results = {}
    overlaps = {}

    # Analyze each system
    systems = ['solvent', 'wt_complex', 'mut_complex']

    for system_name in systems:
        system_dir = os.path.join(BASE_PATH, system_name)
        dG, dG_err, overlap = run_mbar_for_system(system_name, system_dir)

        if dG is not None:
            results[system_name] = (dG, dG_err)
            overlaps[system_name] = overlap
        else:
            print(f"\nERROR: Failed to analyze {system_name}")
            return

    # Load Boresch corrections
    print("\n" + "="*60)
    print("BORESCH RESTRAINT CORRECTIONS")
    print("="*60)

    dG_boresch_wt = load_boresch_correction(
        os.path.join(BASE_PATH, "wt_complex"), "wt_complex")
    dG_boresch_mut = load_boresch_correction(
        os.path.join(BASE_PATH, "mut_complex"), "mut_complex")

    print(f"\n  WT Boresch correction: {dG_boresch_wt:.2f} kcal/mol")
    print(f"  Mut Boresch correction: {dG_boresch_mut:.2f} kcal/mol")

    # Calculate binding free energies
    print("\n" + "="*60)
    print("BINDING FREE ENERGY CALCULATION")
    print("="*60)

    dG_solvent, dG_solvent_err = results['solvent']
    dG_wt_complex, dG_wt_complex_err = results['wt_complex']
    dG_mut_complex, dG_mut_complex_err = results['mut_complex']

    print(f"\nDecoupling free energies:")
    print(f"  ΔG_solvent:     {dG_solvent:8.2f} ± {dG_solvent_err:.2f} kcal/mol")
    print(f"  ΔG_wt_complex:  {dG_wt_complex:8.2f} ± {dG_wt_complex_err:.2f} kcal/mol")
    print(f"  ΔG_mut_complex: {dG_mut_complex:8.2f} ± {dG_mut_complex_err:.2f} kcal/mol")

    # ΔG_bind = ΔG_solvent - ΔG_complex + ΔG_Boresch
    # (This is the standard ABFE thermodynamic cycle)

    dG_bind_wt = dG_solvent - dG_wt_complex + dG_boresch_wt
    dG_bind_wt_err = np.sqrt(dG_solvent_err**2 + dG_wt_complex_err**2)

    dG_bind_mut = dG_solvent - dG_mut_complex + dG_boresch_mut
    dG_bind_mut_err = np.sqrt(dG_solvent_err**2 + dG_mut_complex_err**2)

    print(f"\nBinding free energies:")
    print(f"  ΔG_bind(WT):  {dG_bind_wt:8.2f} ± {dG_bind_wt_err:.2f} kcal/mol")
    print(f"  ΔG_bind(Mut): {dG_bind_mut:8.2f} ± {dG_bind_mut_err:.2f} kcal/mol")

    # ΔΔG_bind = ΔG_bind(Mut) - ΔG_bind(WT)
    # Note: ΔG_solvent cancels out in ΔΔG, so error only depends on complex terms
    # ΔΔG = (ΔG_solvent - ΔG_mut + Boresch_mut) - (ΔG_solvent - ΔG_wt + Boresch_wt)
    #     = ΔG_wt - ΔG_mut + (Boresch_mut - Boresch_wt)
    ddG_bind = dG_bind_mut - dG_bind_wt
    ddG_bind_err = np.sqrt(dG_wt_complex_err**2 + dG_mut_complex_err**2)

    print(f"\n" + "="*60)
    print("FINAL RESULT")
    print("="*60)
    print(f"\n  ΔΔG_bind = ΔG_bind(Mut) - ΔG_bind(WT)")
    print(f"  ΔΔG_bind = {ddG_bind:.2f} ± {ddG_bind_err:.2f} kcal/mol")

    # Interpretation
    print(f"\nInterpretation:")
    if ddG_bind > 0:
        print(f"  Positive ΔΔG: Mutation WEAKENS binding by {ddG_bind:.2f} kcal/mol")
    else:
        print(f"  Negative ΔΔG: Mutation STRENGTHENS binding by {-ddG_bind:.2f} kcal/mol")

    # Compare to MM-GBSA
    mmgbsa_result = -2.7  # kcal/mol from previous calculation
    print(f"\n  MM-GBSA result: {mmgbsa_result:.1f} kcal/mol")
    print(f"  FEP result:     {ddG_bind:.2f} ± {ddG_bind_err:.2f} kcal/mol")

    if np.sign(ddG_bind) == np.sign(mmgbsa_result):
        print(f"  Agreement: Both methods predict mutation {'weakens' if ddG_bind > 0 else 'strengthens'} binding")
    else:
        print(f"  Disagreement: Methods predict opposite effects!")

    # Summary table
    print("\n" + "="*60)
    print("SUMMARY TABLE")
    print("="*60)
    print(f"{'Quantity':<25} {'Value (kcal/mol)':<20} {'Error':<10}")
    print("-"*55)
    print(f"{'ΔG_decouple(solvent)':<25} {dG_solvent:>15.2f}     ±{dG_solvent_err:.2f}")
    print(f"{'ΔG_decouple(wt_complex)':<25} {dG_wt_complex:>15.2f}     ±{dG_wt_complex_err:.2f}")
    print(f"{'ΔG_decouple(mut_complex)':<25} {dG_mut_complex:>15.2f}     ±{dG_mut_complex_err:.2f}")
    print(f"{'ΔG_Boresch(WT)':<25} {dG_boresch_wt:>15.2f}")
    print(f"{'ΔG_Boresch(Mut)':<25} {dG_boresch_mut:>15.2f}")
    print("-"*55)
    print(f"{'ΔG_bind(WT)':<25} {dG_bind_wt:>15.2f}     ±{dG_bind_wt_err:.2f}")
    print(f"{'ΔG_bind(Mut)':<25} {dG_bind_mut:>15.2f}     ±{dG_bind_mut_err:.2f}")
    print("-"*55)
    print(f"{'ΔΔG_bind':<25} {ddG_bind:>15.2f}     ±{ddG_bind_err:.2f}")
    print("="*55)

    # Note about uncertainties
    print("\nNotes:")
    print("- Uncertainties are from MBAR statistical error only")
    print("- Correlation/subsampling not applied (may underestimate error)")
    print("- Boresch corrections are estimates if not pre-computed")

    return results, ddG_bind, ddG_bind_err


if __name__ == "__main__":
    main()
