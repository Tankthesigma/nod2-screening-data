#!/usr/bin/env python
"""MM-GBSA File Discovery Script - Dry Run"""
import MDAnalysis as mda
import os
import warnings
warnings.filterwarnings('ignore')

base = r'C:\Users\vasud\nod2-screening-data'

# Define all 25 simulations
simulations = [
    # Phase 6 - Wild Type (13 sims)
    # Febuxostat
    {'num': 1, 'phase': 'Phase6', 'system': 'WT', 'ligand': 'Febuxostat', 'rep': 1,
     'dcd': r'PHASE_6\trajectories\febuxostat_rep1.dcd',
     'pdb': r'PHASE_6\trajectories\febuxostat_rep1_solvated.pdb',
     'sdf': r'PHASE_6\structures\febuxostat_docked.sdf',
     'xml': r'PHASE_6\checkpoints\febuxostat_rep1_final.xml'},
    {'num': 2, 'phase': 'Phase6', 'system': 'WT', 'ligand': 'Febuxostat', 'rep': 2,
     'dcd': r'PHASE_6\trajectories\febuxostat_rep2.dcd',
     'pdb': r'PHASE_6\trajectories\febuxostat_rep2_solvated.pdb',
     'sdf': r'PHASE_6\structures\febuxostat_docked.sdf',
     'xml': r'PHASE_6\checkpoints\febuxostat_rep2_final.xml'},
    {'num': 3, 'phase': 'Phase6', 'system': 'WT', 'ligand': 'Febuxostat', 'rep': 3,
     'dcd': r'PHASE_6\trajectories\febuxostat_rep3.dcd',
     'pdb': r'PHASE_6\trajectories\febuxostat_rep1_solvated.pdb',  # Substitute rep1 PDB
     'sdf': r'PHASE_6\structures\febuxostat_docked.sdf',
     'xml': r'PHASE_6\checkpoints\febuxostat_rep3_final.xml',
     'note': 'Using rep1 PDB'},
    # Natural
    {'num': 4, 'phase': 'Phase6', 'system': 'WT', 'ligand': 'Natural', 'rep': 1,
     'dcd': r'PHASE_6\trajectories\natural_top_rep1.dcd',
     'pdb': r'PHASE_6\trajectories\natural_top_rep1_solvated.pdb',
     'sdf': r'PHASE_6\structures\natural_top_docked.sdf',
     'xml': r'PHASE_6\checkpoints\natural_top_rep1_final.xml'},
    {'num': 5, 'phase': 'Phase6', 'system': 'WT', 'ligand': 'Natural', 'rep': 2,
     'dcd': r'PHASE_6\trajectories\natural_top_rep2.dcd',
     'pdb': r'PHASE_6\trajectories\natural_top_rep2_solvated_new.pdb',
     'sdf': r'PHASE_6\structures\natural_top_docked.sdf',
     'xml': r'PHASE_6\checkpoints\natural_top_rep2_final.xml',
     'note': 'Using _new PDB'},
    {'num': 6, 'phase': 'Phase6', 'system': 'WT', 'ligand': 'Natural', 'rep': 3,
     'dcd': r'PHASE_6\trajectories\natural_top_rep3.dcd',
     'pdb': r'PHASE_6\trajectories\natural_top_rep3_solvated.pdb',
     'sdf': r'PHASE_6\structures\natural_top_docked.sdf',
     'xml': r'PHASE_6\checkpoints\natural_top_rep3_final.xml'},
    # Budesonide (positive control)
    {'num': 7, 'phase': 'Phase6', 'system': 'WT', 'ligand': 'Budesonide', 'rep': 1,
     'dcd': r'PHASE_6\trajectories\budesonide_rep1.dcd',
     'pdb': r'PHASE_6\trajectories\budesonide_rep1_solvated.pdb',
     'sdf': r'PHASE_6\structures\budesonide_docked.sdf',
     'xml': r'PHASE_6\checkpoints\budesonide_rep1_final.xml'},
    {'num': 8, 'phase': 'Phase6', 'system': 'WT', 'ligand': 'Budesonide', 'rep': 2,
     'dcd': r'PHASE_6\trajectories\budesonide_rep2.dcd',
     'pdb': r'PHASE_6\trajectories\budesonide_rep2_solvated.pdb',
     'sdf': r'PHASE_6\structures\budesonide_docked.sdf',
     'xml': r'PHASE_6\checkpoints\budesonide_rep2_final.xml'},
    {'num': 9, 'phase': 'Phase6', 'system': 'WT', 'ligand': 'Budesonide', 'rep': 3,
     'dcd': r'PHASE_6\trajectories\budesonide_rep3.dcd',
     'pdb': r'PHASE_6\trajectories\budesonide_rep3_solvated.pdb',
     'sdf': r'PHASE_6\structures\budesonide_docked.sdf',
     'xml': r'PHASE_6\checkpoints\budesonide_rep3_final.xml'},
    # Ursodiol
    {'num': 10, 'phase': 'Phase6', 'system': 'WT', 'ligand': 'Ursodiol', 'rep': 1,
     'dcd': r'PHASE_6\trajectories\ursodiol_rep1.dcd',
     'pdb': r'PHASE_6\trajectories\ursodiol_rep3_solvated.pdb',
     'sdf': r'PHASE_6\structures\ursodiol_docked.sdf',
     'xml': r'PHASE_6\checkpoints\ursodiol_rep1_final.xml',
     'note': 'Using rep3 PDB'},
    {'num': 11, 'phase': 'Phase6', 'system': 'WT', 'ligand': 'Ursodiol', 'rep': 2,
     'dcd': r'PHASE_6\trajectories\ursodiol_rep2.dcd',
     'pdb': r'PHASE_6\trajectories\ursodiol_rep3_solvated.pdb',
     'sdf': r'PHASE_6\structures\ursodiol_docked.sdf',
     'xml': r'PHASE_6\checkpoints\ursodiol_rep2_final.xml',
     'note': 'Using rep3 PDB'},
    {'num': 12, 'phase': 'Phase6', 'system': 'WT', 'ligand': 'Ursodiol', 'rep': 3,
     'dcd': r'PHASE_6\trajectories\ursodiol_rep3.dcd',
     'pdb': r'PHASE_6\trajectories\ursodiol_rep3_solvated.pdb',
     'sdf': r'PHASE_6\structures\ursodiol_docked.sdf',
     'xml': r'PHASE_6\checkpoints\ursodiol_rep3_final.xml'},
    # Decoy (negative control)
    {'num': 13, 'phase': 'Phase6', 'system': 'WT', 'ligand': 'Decoy', 'rep': 1,
     'dcd': r'PHASE_6\trajectories\decoy_rep1.dcd',
     'pdb': r'PHASE_6\trajectories\decoy_rep1_solvated.pdb',
     'sdf': r'PHASE_6\structures\decoy_docked.sdf',
     'xml': r'PHASE_6\checkpoints\decoy_rep1_final.xml',
     'note': 'NEGATIVE CONTROL'},

    # Phase A1 - Mutants (12 sims)
    # R702W Febuxostat
    {'num': 14, 'phase': 'PhaseA1', 'system': 'R702W', 'ligand': 'Febuxostat', 'rep': 1,
     'dcd': r'PHASE_A1_mutant_MD\trajectories\R702W_febuxostat_rep1.dcd',
     'pdb': r'PHASE_A1_mutant_MD\trajectories\R702W_febuxostat_rep1_solvated.pdb',
     'sdf': r'PHASE_A1_mutant_MD\ligands\febuxostat_docked.sdf',
     'xml': r'PHASE_A1_mutant_MD\checkpoints\R702W_febuxostat_rep1_final.xml'},
    {'num': 15, 'phase': 'PhaseA1', 'system': 'R702W', 'ligand': 'Febuxostat', 'rep': 2,
     'dcd': r'PHASE_A1_mutant_MD\trajectories\R702W_febuxostat_rep2.dcd',
     'pdb': r'PHASE_A1_mutant_MD\trajectories\R702W_febuxostat_rep2_solvated.pdb',
     'sdf': r'PHASE_A1_mutant_MD\ligands\febuxostat_docked.sdf',
     'xml': r'PHASE_A1_mutant_MD\checkpoints\R702W_febuxostat_rep2_final.xml'},
    {'num': 16, 'phase': 'PhaseA1', 'system': 'R702W', 'ligand': 'Febuxostat', 'rep': 3,
     'dcd': r'PHASE_A1_mutant_MD\trajectories\R702W_febuxostat_rep3.dcd',
     'pdb': r'PHASE_A1_mutant_MD\trajectories\R702W_febuxostat_rep3_solvated.pdb',
     'sdf': r'PHASE_A1_mutant_MD\ligands\febuxostat_docked.sdf',
     'xml': r'PHASE_A1_mutant_MD\checkpoints\R702W_febuxostat_rep3_final.xml'},
    # R702W Natural
    {'num': 17, 'phase': 'PhaseA1', 'system': 'R702W', 'ligand': 'Natural', 'rep': 1,
     'dcd': r'PHASE_A1_mutant_MD\vast_downloads\5080\R702W_natural_rep1.dcd',
     'pdb': r'PHASE_A1_mutant_MD\vast_downloads\5080\R702W_natural_rep1_solvated.pdb',
     'sdf': r'PHASE_A1_mutant_MD\ligands\natural_top_docked.sdf',
     'xml': None},
    {'num': 18, 'phase': 'PhaseA1', 'system': 'R702W', 'ligand': 'Natural', 'rep': 2,
     'dcd': r'PHASE_A1_mutant_MD\vast_downloads\5080\R702W_natural_rep2.dcd',
     'pdb': r'PHASE_A1_mutant_MD\vast_downloads\5080\R702W_natural_rep2_solvated.pdb',
     'sdf': r'PHASE_A1_mutant_MD\ligands\natural_top_docked.sdf',
     'xml': None},
    {'num': 19, 'phase': 'PhaseA1', 'system': 'R702W', 'ligand': 'Natural', 'rep': 3,
     'dcd': r'PHASE_A1_mutant_MD\vast_downloads\new_5090\R702W_natural_rep3.dcd',
     'pdb': r'PHASE_A1_mutant_MD\vast_downloads\new_5090\R702W_natural_rep3_solvated.pdb',
     'sdf': r'PHASE_A1_mutant_MD\ligands\natural_top_docked.sdf',
     'xml': None},
    # G908R Febuxostat
    {'num': 20, 'phase': 'PhaseA1', 'system': 'G908R', 'ligand': 'Febuxostat', 'rep': 1,
     'dcd': r'PHASE_A1_mutant_MD\trajectories\G908R_febuxostat_rep1.dcd',
     'pdb': r'PHASE_A1_mutant_MD\trajectories\G908R_febuxostat_rep1_solvated.pdb',
     'sdf': r'PHASE_A1_mutant_MD\ligands\febuxostat_docked.sdf',
     'xml': r'PHASE_A1_mutant_MD\checkpoints\G908R_febuxostat_rep1_final.xml',
     'note': 'OUTLIER 45.8% occ'},
    {'num': 21, 'phase': 'PhaseA1', 'system': 'G908R', 'ligand': 'Febuxostat', 'rep': 2,
     'dcd': r'PHASE_A1_mutant_MD\vast_downloads\5080\G908R_febuxostat_rep2.dcd',
     'pdb': r'PHASE_A1_mutant_MD\vast_downloads\5080\G908R_febuxostat_rep2_solvated.pdb',
     'sdf': r'PHASE_A1_mutant_MD\ligands\febuxostat_docked.sdf',
     'xml': None},
    {'num': 22, 'phase': 'PhaseA1', 'system': 'G908R', 'ligand': 'Febuxostat', 'rep': 3,
     'dcd': r'PHASE_A1_mutant_MD\vast_downloads\5080\G908R_febuxostat_rep3.dcd',
     'pdb': r'PHASE_A1_mutant_MD\vast_downloads\5080\G908R_febuxostat_rep3_solvated.pdb',
     'sdf': r'PHASE_A1_mutant_MD\ligands\febuxostat_docked.sdf',
     'xml': None},
    # G908R Natural
    {'num': 23, 'phase': 'PhaseA1', 'system': 'G908R', 'ligand': 'Natural', 'rep': 1,
     'dcd': r'PHASE_A1_mutant_MD\vast_downloads\5090\G908R_natural_rep1.dcd',
     'pdb': r'PHASE_A1_mutant_MD\vast_downloads\5090\G908R_natural_rep1_solvated.pdb',
     'sdf': r'PHASE_A1_mutant_MD\ligands\natural_top_docked.sdf',
     'xml': None},
    {'num': 24, 'phase': 'PhaseA1', 'system': 'G908R', 'ligand': 'Natural', 'rep': 2,
     'dcd': r'PHASE_A1_mutant_MD\vast_downloads\new_5090\G908R_natural_rep2.dcd',
     'pdb': r'PHASE_A1_mutant_MD\vast_downloads\new_5090\G908R_natural_rep2_solvated.pdb',
     'sdf': r'PHASE_A1_mutant_MD\ligands\natural_top_docked.sdf',
     'xml': None},
    {'num': 25, 'phase': 'PhaseA1', 'system': 'G908R', 'ligand': 'Natural', 'rep': 3,
     'dcd': r'PHASE_A1_mutant_MD\vast_downloads\new_5090\G908R_natural_rep3.dcd',
     'pdb': r'PHASE_A1_mutant_MD\vast_downloads\new_5090\G908R_natural_rep3_solvated.pdb',
     'sdf': r'PHASE_A1_mutant_MD\ligands\natural_top_docked.sdf',
     'xml': None},
]

print("=" * 130)
print("MM-GBSA FILE DISCOVERY - DRY RUN")
print("=" * 130)
print(f"{'#':<3} {'Phase':<8} {'System':<7} {'Ligand':<12} {'Rep':<4} {'DCD':<5} {'PDB':<5} {'SDF':<5} {'XML':<5} {'Frames':<7} {'Status':<25}")
print("-" * 130)

ok_count = 0
fail_count = 0
issues = []

for sim in simulations:
    dcd_path = os.path.join(base, sim['dcd']) if sim['dcd'] else None
    pdb_path = os.path.join(base, sim['pdb']) if sim['pdb'] else None
    sdf_path = os.path.join(base, sim['sdf']) if sim['sdf'] else None
    xml_path = os.path.join(base, sim['xml']) if sim.get('xml') else None

    dcd_ok = os.path.exists(dcd_path) if dcd_path else False
    pdb_ok = os.path.exists(pdb_path) if pdb_path else False
    sdf_ok = os.path.exists(sdf_path) if sdf_path else False
    xml_ok = os.path.exists(xml_path) if xml_path else False

    frames = '-'
    if dcd_ok and pdb_ok:
        try:
            u = mda.Universe(pdb_path, dcd_path)
            frames = str(len(u.trajectory))
        except Exception as e:
            frames = 'ERR'
            issues.append(f"#{sim['num']}: Load error - {str(e)[:50]}")

    status = 'OK'
    if not dcd_ok:
        status = 'MISSING DCD'
        issues.append(f"#{sim['num']}: Missing DCD - {sim['dcd']}")
    elif not pdb_ok:
        status = 'MISSING PDB'
        issues.append(f"#{sim['num']}: Missing PDB - {sim['pdb']}")
    elif not sdf_ok:
        status = 'MISSING SDF'
        issues.append(f"#{sim['num']}: Missing SDF - {sim['sdf']}")

    note = sim.get('note', '')
    if note and status == 'OK':
        status = note

    if dcd_ok and pdb_ok and sdf_ok:
        ok_count += 1
    else:
        fail_count += 1

    dcd_mark = 'Y' if dcd_ok else 'N'
    pdb_mark = 'Y' if pdb_ok else 'N'
    sdf_mark = 'Y' if sdf_ok else 'N'
    xml_mark = 'Y' if xml_ok else '-'

    print(f"{sim['num']:<3} {sim['phase']:<8} {sim['system']:<7} {sim['ligand']:<12} {sim['rep']:<4} "
          f"{dcd_mark:<5} {pdb_mark:<5} {sdf_mark:<5} {xml_mark:<5} {frames:<7} {status:<25}")

print("-" * 130)
print(f"\nSUMMARY: {ok_count} OK, {fail_count} FAILED out of 25 simulations")

if issues:
    print(f"\nISSUES ({len(issues)}):")
    for issue in issues:
        print(f"  - {issue}")

print("\n" + "=" * 130)
print("AWAITING APPROVAL TO PROCEED TO PILOT TEST")
print("=" * 130)
