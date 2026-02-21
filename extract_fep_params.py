#!/usr/bin/env python
"""Extract FEP parameters from febuxostat setup for reuse with natural product."""
import numpy as np
import os
import json

BASE = "C:/Users/vasud/nod2-screening-data/fep_complete/fep_pmx"

print("="*70)
print("EXTRACTING FEBUXOSTAT FEP PARAMETERS")
print("="*70)

# Lambda schedule
sched = np.load(f"{BASE}/wt_complex/lambda_schedule.npy")
print(f"\nLambda schedule (K={len(sched)} windows):")
print("Window  Elec    Sterics  Restraints")
print("-"*40)
for i in range(len(sched)):
    print(f"{i:2d}      {sched[i,0]:.4f}  {sched[i,1]:.4f}   {sched[i,2]:.4f}")

# Save lambda arrays
lambda_elec = sched[:, 0].tolist()
lambda_sterics = sched[:, 1].tolist()
lambda_restraints = sched[:, 2].tolist()

# Boresch anchors
print("\n" + "="*70)
print("BORESCH ANCHORS")
print("="*70)
for sys_name in ['wt_complex', 'mut_complex']:
    anchors_file = f"{BASE}/{sys_name}/boresch_anchors.npy"
    if os.path.exists(anchors_file):
        anchors = np.load(anchors_file, allow_pickle=True).item()
        print(f"\n{sys_name}:")
        for k, v in anchors.items():
            print(f"  {k}: {v}")

# Check alchemical system for softcore parameters
print("\n" + "="*70)
print("ALCHEMICAL SYSTEM PARAMETERS")
print("="*70)
xml_file = f"{BASE}/wt_complex/alchemical_system.xml"
with open(xml_file, 'r') as f:
    xml_content = f.read()

# Extract key parameters
if 'alpha' in xml_content.lower():
    # Look for softcore alpha
    import re
    alpha_match = re.search(r'softcore.*?alpha["\s=:]+([0-9.]+)', xml_content, re.IGNORECASE)
    if alpha_match:
        print(f"Softcore alpha: {alpha_match.group(1)}")

# Check for CustomBondForce (Boresch restraints)
if 'CustomBondForce' in xml_content:
    print("CustomBondForce present (Boresch restraints)")
if 'CustomAngleForce' in xml_content:
    print("CustomAngleForce present (Boresch restraints)")
if 'CustomTorsionForce' in xml_content:
    print("CustomTorsionForce present (Boresch restraints)")

# Count atoms in topology
print("\n" + "="*70)
print("SYSTEM SIZES")
print("="*70)
for sys_name in ['wt_complex', 'mut_complex', 'solvent']:
    pdb_file = f"{BASE}/{sys_name}/topology.pdb"
    if os.path.exists(pdb_file):
        with open(pdb_file, 'r') as f:
            atom_count = sum(1 for line in f if line.startswith('ATOM') or line.startswith('HETATM'))
        print(f"{sys_name}: {atom_count} atoms")

# Check ion counts
print("\n" + "="*70)
print("ION COUNTS")
print("="*70)
for sys_name in ['wt_complex', 'mut_complex', 'solvent']:
    pdb_file = f"{BASE}/{sys_name}/topology.pdb"
    if os.path.exists(pdb_file):
        with open(pdb_file, 'r') as f:
            content = f.read()
        na_count = content.count(' NA ')
        cl_count = content.count(' CL ')
        print(f"{sys_name}: Na+={na_count}, Cl-={cl_count}")

# Save parameters to JSON
params = {
    'lambda_electrostatics': lambda_elec,
    'lambda_sterics': lambda_sterics,
    'lambda_restraints': lambda_restraints,
    'n_windows': len(sched),
    'temperature_K': 310.0,
    'friction_per_ps': 1.0,
    'timestep_fs': 2.0,
    'equil_steps': 50000,
    'prod_steps': 500000,
    'energy_interval': 500,
}

with open('C:/Users/vasud/nod2-screening-data/febuxostat_fep_params.json', 'w') as f:
    json.dump(params, f, indent=2)
print("\n[SAVED] febuxostat_fep_params.json")

print("\n" + "="*70)
print("EXTRACTION COMPLETE")
print("="*70)
