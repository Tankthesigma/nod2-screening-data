#!/usr/bin/env python3
"""
PHASE A1 ANALYSIS - Report Generator
Compiles all analysis results into a final markdown report.
"""

import os
import sys
from pathlib import Path
import platform
from datetime import datetime

if 'microsoft' in platform.uname().release.lower() or os.path.exists('/mnt/c'):
    BASE_DIR = Path("/mnt/c/Users/vasud/nod2-screening-data/PHASE_A1_mutant_MD")
else:
    BASE_DIR = Path(r"C:\Users\vasud\nod2-screening-data\PHASE_A1_mutant_MD")

def read_csv(filepath):
    """Read CSV file and return list of dicts."""
    if not filepath.exists():
        return None

    with open(filepath, 'r') as f:
        lines = f.readlines()

    if len(lines) < 2:
        return None

    headers = lines[0].strip().split(',')
    data = []
    for line in lines[1:]:
        values = line.strip().split(',')
        row = dict(zip(headers, values))
        data.append(row)

    return data

def generate_report():
    analysis_dir = BASE_DIR / "analysis"

    # Read all CSV files
    pocket_data = read_csv(analysis_dir / "pocket_occupancy.csv")
    rmsd_data = read_csv(analysis_dir / "ligand_rmsd.csv")
    contact_data = read_csv(analysis_dir / "key_residue_contacts.csv")
    hbond_data = read_csv(analysis_dir / "hbond_summary.csv")
    binding_site_data = read_csv(analysis_dir / "binding_site_rmsd.csv")
    com_data = read_csv(analysis_dir / "ligand_com_distance.csv")
    energy_data = read_csv(analysis_dir / "energy_analysis.csv")

    report = []
    report.append("# PHASE A1 ANALYSIS REPORT")
    report.append(f"## NOD2 Mutant MD Simulations")
    report.append(f"**Generated:** {datetime.now().strftime('%Y-%m-%d %H:%M')}")
    report.append("")

    # Executive Summary
    report.append("---")
    report.append("## Executive Summary")
    report.append("")

    if pocket_data:
        # Calculate averages by mutant/ligand
        groups = {}
        for row in pocket_data:
            key = (row['mutant'], row['ligand'])
            if key not in groups:
                groups[key] = []
            try:
                groups[key].append(float(row['occupancy_percent']))
            except:
                pass

        report.append("### Pocket Occupancy Results")
        report.append("")
        report.append("| Mutant | Ligand | Avg Occupancy | Verdict |")
        report.append("|--------|--------|---------------|---------|")

        for (mutant, ligand), values in sorted(groups.items()):
            if values:
                avg = sum(values) / len(values)
                if avg >= 95:
                    verdict = "BINDING PRESERVED"
                elif avg >= 80:
                    verdict = "~ PARTIAL BINDING"
                else:
                    verdict = "BINDING IMPAIRED"
                report.append(f"| {mutant} | {ligand} | {avg:.1f}% | {verdict} |")

        report.append("")

    # Methods Section
    report.append("---")
    report.append("## Methods")
    report.append("")
    report.append("### Simulation Parameters")
    report.append("- **Mutations:** R702W, G908R (Crohn's disease-associated)")
    report.append("- **Ligands:** Febuxostat (drug), Natural (MDP control)")
    report.append("- **Replicates:** 3 per condition")
    report.append("- **Duration:** 20 ns each (240 ns total)")
    report.append("- **Force field:** AMBER14 (protein) + OpenFF-2.1.0 (ligand)")
    report.append("- **Temperature:** 310.15 K")
    report.append("")
    report.append("### Analysis Definitions")
    report.append("- **Pocket residues:** GLU1008, ASP1011, ARG1037")
    report.append("- **Pocket occupancy cutoff:** 5.0 Å")
    report.append("- **Contact cutoff:** 4.0 Å")
    report.append("- **H-bond criteria:** 3.5 A, >=150 deg")
    report.append("")

    # Detailed Results
    report.append("---")
    report.append("## Detailed Results")
    report.append("")

    # Table 1: Pocket Occupancy
    if pocket_data:
        report.append("### Table 1: Pocket Occupancy")
        report.append("")
        report.append("| ID | Mutant | Ligand | Rep | Occupancy | Status |")
        report.append("|----|--------|--------|-----|-----------|--------|")
        for i, row in enumerate(pocket_data, 1):
            occ = f"{float(row['occupancy_percent']):.1f}%" if row['occupancy_percent'] != '' else "N/A"
            report.append(f"| {i} | {row['mutant']} | {row['ligand']} | {row['rep']} | {occ} | {row['status']} |")
        report.append("")

    # Table 2: Ligand RMSD
    if rmsd_data:
        report.append("### Table 2: Ligand RMSD")
        report.append("")
        report.append("| ID | Mutant | Ligand | Rep | Mean (Å) | SD (Å) | Max (Å) | Status |")
        report.append("|----|--------|--------|-----|----------|--------|---------|--------|")
        for i, row in enumerate(rmsd_data, 1):
            mean = f"{float(row['mean_rmsd']):.2f}" if row['mean_rmsd'] != '' else "N/A"
            std = f"{float(row['std_rmsd']):.2f}" if row.get('std_rmsd', '') != '' else "N/A"
            max_r = f"{float(row['max_rmsd']):.2f}" if row.get('max_rmsd', '') != '' else "N/A"
            report.append(f"| {i} | {row['mutant']} | {row['ligand']} | {row['rep']} | {mean} | {std} | {max_r} | {row['status']} |")
        report.append("")

    # Table 3: Key Residue Contacts
    if contact_data:
        report.append("### Table 3: Key Residue Contact Frequencies")
        report.append("")
        report.append("| Mutant | Ligand | Rep | GLU1008 | ASP1011 | ARG1037 |")
        report.append("|--------|--------|-----|---------|---------|---------|")

        # Group by sim
        sims = {}
        for row in contact_data:
            key = (row['mutant'], row['ligand'], row['rep'])
            if key not in sims:
                sims[key] = {}
            sims[key][int(row['resid'])] = float(row['frequency_percent'])

        for (mutant, ligand, rep), residues in sorted(sims.items()):
            g1008 = f"{residues.get(1008, 0):.1f}%"
            a1011 = f"{residues.get(1011, 0):.1f}%"
            r1037 = f"{residues.get(1037, 0):.1f}%"
            report.append(f"| {mutant} | {ligand} | {rep} | {g1008} | {a1011} | {r1037} |")
        report.append("")

    # Table 4: Binding Site RMSD
    if binding_site_data:
        report.append("### Table 4: Binding Site RMSD")
        report.append("")
        report.append("| ID | Mutant | Ligand | Rep | Mean (Å) | Status |")
        report.append("|----|--------|--------|-----|----------|--------|")
        for i, row in enumerate(binding_site_data, 1):
            mean = f"{float(row['mean_rmsd']):.2f}" if row['mean_rmsd'] != '' else "N/A"
            report.append(f"| {i} | {row['mutant']} | {row['ligand']} | {row['rep']} | {mean} | {row['status']} |")
        report.append("")

    # Table 5: Hydrogen Bond Summary
    if hbond_data:
        report.append("### Table 5: Hydrogen Bond Summary")
        report.append("")
        report.append("| ID | Mutant | Ligand | Rep | Total Events | Unique Pairs | Persistent (>25%) | Status |")
        report.append("|----|--------|--------|-----|--------------|--------------|-------------------|--------|")
        for i, row in enumerate(hbond_data, 1):
            total = row.get('total_hbonds', 'N/A')
            unique = row.get('unique_pairs', 'N/A')
            persistent = row.get('persistent', 'N/A')
            report.append(f"| {i} | {row['mutant']} | {row['ligand']} | {row['rep']} | {total} | {unique} | {persistent} | {row['status']} |")
        report.append("")

    # Table 6: Ligand COM Distance
    if com_data:
        report.append("### Table 6: Ligand COM Distance from Pocket")
        report.append("")
        report.append("| ID | Mutant | Ligand | Rep | Mean (Å) | SD (Å) | Drift (Å) | Status |")
        report.append("|----|--------|--------|-----|----------|--------|-----------|--------|")
        for i, row in enumerate(com_data, 1):
            mean = f"{float(row['mean_dist']):.2f}" if row['mean_dist'] != '' else "N/A"
            std = f"{float(row['std_dist']):.2f}" if row.get('std_dist', '') != '' else "N/A"
            drift = f"{float(row['drift']):+.2f}" if row.get('drift', '') != '' else "N/A"
            report.append(f"| {i} | {row['mutant']} | {row['ligand']} | {row['rep']} | {mean} | {std} | {drift} | {row['status']} |")
        report.append("")

    # Table 7: Energy Analysis
    if energy_data:
        report.append("### Table 7: Temperature Stability")
        report.append("")
        report.append("| ID | Mutant | Ligand | Rep | Temp (K) | Status |")
        report.append("|----|--------|--------|-----|----------|--------|")
        for i, row in enumerate(energy_data, 1):
            if row['temp_mean'] != '':
                temp = f"{float(row['temp_mean']):.1f}±{float(row['temp_std']):.1f}"
            else:
                temp = "N/A"
            report.append(f"| {i} | {row['mutant']} | {row['ligand']} | {row['rep']} | {temp} | {row['status']} |")
        report.append("")

    # Conclusions
    report.append("---")
    report.append("## Conclusions")
    report.append("")
    report.append("### Key Findings")
    report.append("")

    if pocket_data:
        # Analyze results
        feb_r702w = [float(r['occupancy_percent']) for r in pocket_data if r['mutant'] == 'R702W' and r['ligand'] == 'febuxostat' and r['occupancy_percent'] != '']
        feb_g908r = [float(r['occupancy_percent']) for r in pocket_data if r['mutant'] == 'G908R' and r['ligand'] == 'febuxostat' and r['occupancy_percent'] != '']
        nat_r702w = [float(r['occupancy_percent']) for r in pocket_data if r['mutant'] == 'R702W' and r['ligand'] == 'natural' and r['occupancy_percent'] != '']
        nat_g908r = [float(r['occupancy_percent']) for r in pocket_data if r['mutant'] == 'G908R' and r['ligand'] == 'natural' and r['occupancy_percent'] != '']

        report.append("1. **R702W Mutant:**")
        if feb_r702w:
            avg = sum(feb_r702w) / len(feb_r702w)
            report.append(f"   - Febuxostat binding: {avg:.1f}% occupancy")
        if nat_r702w:
            avg = sum(nat_r702w) / len(nat_r702w)
            report.append(f"   - Natural ligand binding: {avg:.1f}% occupancy")
        report.append("")

        report.append("2. **G908R Mutant:**")
        if feb_g908r:
            avg = sum(feb_g908r) / len(feb_g908r)
            report.append(f"   - Febuxostat binding: {avg:.1f}% occupancy")
        if nat_g908r:
            avg = sum(nat_g908r) / len(nat_g908r)
            report.append(f"   - Natural ligand binding: {avg:.1f}% occupancy")
        report.append("")

    report.append("### Precision Medicine Implications")
    report.append("")
    report.append("Based on pocket occupancy results:")
    report.append("- **Occupancy >=95%:** Febuxostat binding preserved, patient likely responsive")
    report.append("- **Occupancy 80-95%:** Partial binding, monitor treatment response")
    report.append("- **Occupancy <80%:** Binding impaired, consider alternative therapy")
    report.append("")

    report.append("---")
    report.append("## References")
    report.append("")
    report.append("- Wild-type NOD2 baseline: 99-100% pocket occupancy (Phase 6/7)")
    report.append("- R702W mutation: ~10% of Crohn's patients")
    report.append("- G908R mutation: ~4.6% of Crohn's patients")
    report.append("")

    # Write report
    report_path = analysis_dir / "PHASE_A1_ANALYSIS_REPORT.md"
    with open(report_path, 'w') as f:
        f.write('\n'.join(report))

    print(f"Report saved to: {report_path}")

    # Also create summary tables file
    summary_path = analysis_dir / "SUMMARY_TABLES.md"
    with open(summary_path, 'w') as f:
        f.write("# PHASE A1 SUMMARY TABLES\n\n")
        f.write("See PHASE_A1_ANALYSIS_REPORT.md for full report.\n")

    print(f"Summary saved to: {summary_path}")

    return report_path

def main():
    print("=" * 80)
    print("PHASE A1 ANALYSIS - REPORT GENERATOR")
    print("=" * 80)
    print()

    report_path = generate_report()

    print()
    print("=" * 80)
    print("REPORT GENERATION COMPLETE")
    print("=" * 80)

if __name__ == "__main__":
    main()
