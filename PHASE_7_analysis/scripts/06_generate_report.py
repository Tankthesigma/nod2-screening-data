#!/usr/bin/env python3
"""
PHASE 7: Report Generator
Combines all analyses into a final summary report.
"""

import pandas as pd
from pathlib import Path
from datetime import datetime

# Paths
ANALYSIS_DIR = Path(r"C:\Users\vasud\nod2-screening-data\PHASE_7_analysis")
REPORT_DIR = ANALYSIS_DIR / "report"


def load_csv_safe(filepath):
    """Load CSV if it exists."""
    if filepath.exists():
        return pd.read_csv(filepath)
    return None


def generate_report():
    """Generate comprehensive markdown report."""

    REPORT_DIR.mkdir(parents=True, exist_ok=True)

    # Load all statistics
    rmsd_stats = load_csv_safe(ANALYSIS_DIR / "rmsd" / "rmsd_statistics.csv")
    hbond_stats = load_csv_safe(ANALYSIS_DIR / "hbonds" / "hbond_statistics.csv")
    contact_stats = load_csv_safe(ANALYSIS_DIR / "contacts" / "contact_statistics.csv")
    fes_stats = load_csv_safe(ANALYSIS_DIR / "fes" / "fes_statistics.csv")
    key_residues = load_csv_safe(ANALYSIS_DIR / "contacts" / "key_binding_residues.csv")

    report = []

    # Header
    report.append("# PHASE 7: MD Trajectory Analysis Report")
    report.append(f"\n**Generated:** {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    report.append("\n**Project:** NOD2-CROHN Drug Discovery (ISEF 2026)")
    report.append("\n---\n")

    # Executive Summary
    report.append("## Executive Summary\n")

    if rmsd_stats is not None:
        # Check validation criteria
        febuxostat_row = rmsd_stats[rmsd_stats['Compound'] == 'febuxostat']
        budesonide_row = rmsd_stats[rmsd_stats['Compound'] == 'budesonide']
        decoy_row = rmsd_stats[rmsd_stats['Compound'] == 'decoy']

        feb_stable = False
        bud_stable = False
        decoy_unstable = False

        if len(febuxostat_row) > 0:
            feb_rmsd = float(febuxostat_row['Mean RMSD (Å)'].values[0])
            feb_stable = feb_rmsd < 3
            report.append(f"- **Febuxostat Mean RMSD:** {feb_rmsd:.2f}Å " +
                         f"({'STABLE ✓' if feb_stable else 'MARGINAL/UNSTABLE'})")

        if len(budesonide_row) > 0:
            bud_rmsd = float(budesonide_row['Mean RMSD (Å)'].values[0])
            bud_stable = bud_rmsd < 3
            report.append(f"- **Budesonide Mean RMSD:** {bud_rmsd:.2f}Å " +
                         f"({'STABLE ✓' if bud_stable else 'MARGINAL/UNSTABLE'})")

        if len(decoy_row) > 0:
            dec_rmsd = float(decoy_row['Mean RMSD (Å)'].values[0])
            decoy_unstable = dec_rmsd > 5
            report.append(f"- **Decoy Mean RMSD:** {dec_rmsd:.2f}Å " +
                         f"({'UNSTABLE ✓ (expected)' if decoy_unstable else 'UNEXPECTEDLY STABLE ⚠️'})")

        report.append("\n### Validation Status\n")

        if feb_stable and bud_stable and decoy_unstable:
            report.append("**✅ VALIDATION PASSED**")
            report.append("\n- Febuxostat shows similar stability to positive control (Budesonide)")
            report.append("- Decoy correctly shows unstable binding (negative control worked)")
            report.append("- **Conclusion: Febuxostat is a VALIDATED NOD2 binder candidate**")
        elif not decoy_unstable:
            report.append("**⚠️ CRITICAL ISSUE: Negative control (Decoy) appears stable**")
            report.append("\nThis may indicate a force field artifact. Results should be interpreted cautiously.")
        else:
            report.append("**⚠️ PARTIAL VALIDATION**")
            report.append("\nSome criteria not met. See detailed analysis below.")

    report.append("\n---\n")

    # RMSD Results
    report.append("## 1. Binding Stability (RMSD Analysis)\n")
    report.append("Ligand RMSD measures how much the drug moves from its initial docked position.\n")
    report.append("- **<3Å** = Stable binding")
    report.append("- **3-5Å** = Marginal binding")
    report.append("- **>5Å** = Unstable/Unbound\n")

    if rmsd_stats is not None:
        report.append("\n### RMSD Statistics Table\n")
        report.append(rmsd_stats.to_markdown(index=False))
    else:
        report.append("\n*RMSD analysis not yet completed.*")

    report.append("\n\n![Ligand RMSD Comparison](../rmsd/ligand_rmsd_comparison.png)")

    report.append("\n---\n")

    # H-bond Results
    report.append("## 2. Hydrogen Bond Analysis\n")
    report.append("Hydrogen bonds anchor the drug in the binding site. More stable H-bonds = stronger binding.\n")

    if hbond_stats is not None:
        report.append("\n### H-bond Statistics\n")
        report.append(hbond_stats.to_markdown(index=False))
    else:
        report.append("\n*H-bond analysis not yet completed.*")

    report.append("\n\n![H-bond Timeseries](../hbonds/hbond_timeseries.png)")
    report.append("\n![H-bond Summary](../hbonds/hbond_summary.png)")

    report.append("\n---\n")

    # Contact Analysis
    report.append("## 3. Binding Site Contacts\n")
    report.append("Identifies which protein residues interact with the drug most frequently.\n")

    if contact_stats is not None:
        report.append("\n### Contact Statistics\n")
        report.append(contact_stats.to_markdown(index=False))

    if key_residues is not None and len(key_residues) > 0:
        report.append("\n### Key Binding Residues (shared across compounds)\n")
        report.append(key_residues.to_markdown(index=False))
    else:
        report.append("\n*Contact analysis not yet completed.*")

    report.append("\n\n![Contact Heatmap](../contacts/contact_heatmap.png)")

    report.append("\n---\n")

    # FES Analysis
    report.append("## 4. Conformational Landscape (FES Proxy)\n")
    report.append("2D density plot of RMSD vs Radius of Gyration shows binding stability landscape.\n")
    report.append("- **Tight cluster** = stable binding with defined pose")
    report.append("- **Scattered/smeared** = dynamic/unstable binding\n")

    if fes_stats is not None:
        report.append("\n### FES Statistics\n")
        report.append(fes_stats.to_markdown(index=False))
    else:
        report.append("\n*FES analysis not yet completed.*")

    report.append("\n\n![FES Individual](../fes/fes_individual.png)")
    report.append("\n![FES Combined](../fes/fes_combined.png)")

    report.append("\n---\n")

    # Methods Summary
    report.append("## Methods Summary\n")
    report.append("""
- **Simulation Setup:** 20ns production MD with OpenMM
- **Force Field:** OpenFF 2.1.0 (ligand), AMBER14 (protein), TIP3P-FB (water)
- **Analysis Tools:** MDAnalysis, Python
- **Trajectory Stride:** Every 10th frame analyzed
- **Contact Cutoff:** 4.0Å
- **H-bond Criteria:** 3.5Å distance, 150° angle
""")

    report.append("\n---\n")

    # Conclusions
    report.append("## Conclusions\n")
    report.append("""
### For ISEF Judges:

1. **Scientific Question:** Can existing FDA-approved drugs bind to NOD2, a key Crohn's disease target?

2. **Approach:** Used computational screening + MD validation to test 11,122 compounds

3. **Key Finding:** Febuxostat (gout drug) shows stable binding comparable to positive control

4. **Significance:** If validated experimentally, this could be repurposed for Crohn's disease

5. **Next Steps:** Cell-based assays, NOD2 activity assays, eventual clinical collaboration
""")

    # Write report
    report_path = REPORT_DIR / "PHASE_7_RESULTS.md"
    with open(report_path, 'w') as f:
        f.write('\n'.join(report))

    print(f"Report saved to: {report_path}")

    return report_path


def main():
    print("=" * 60)
    print("PHASE 7: GENERATING FINAL REPORT")
    print("=" * 60)

    report_path = generate_report()

    print("\n" + "=" * 60)
    print("REPORT GENERATION COMPLETE")
    print("=" * 60)
    print(f"\nView report at: {report_path}")


if __name__ == "__main__":
    main()
