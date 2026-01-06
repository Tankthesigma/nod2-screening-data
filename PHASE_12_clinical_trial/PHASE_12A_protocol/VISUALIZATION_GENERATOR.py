#!/usr/bin/env python3
"""
PHASE 12: Clinical Trial Visualization Generator

Generates:
1. PATIENT_FLOWCHART.png - Screening to randomization flow
2. TRIAL_TIMELINE.png - 27-month Gantt chart
"""

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
import math
import os

# Output directory
OUTPUT_DIR = os.path.dirname(os.path.abspath(__file__))

# Calculate sample size dynamically (same logic as SAMPLE_SIZE_CALCULATION.py)
def get_sample_size_per_arm():
    """Calculate N per arm using same parameters as main script."""
    from scipy import stats
    alpha = 0.05
    power = 0.80
    p1 = 0.25
    p2 = 0.50
    dropout = 0.15

    z_alpha = stats.norm.ppf(1 - alpha/2)
    z_beta = stats.norm.ppf(power)
    p_bar = (p1 + p2) / 2
    effect = p2 - p1

    n_raw = (2 * ((z_alpha + z_beta) ** 2) * p_bar * (1 - p_bar)) / (effect ** 2)
    n_adjusted = math.ceil(math.ceil(n_raw) / (1 - dropout))
    return n_adjusted

N_PER_ARM = get_sample_size_per_arm()


def create_patient_flowchart():
    """Create patient flow diagram from screening to randomization."""

    fig, ax = plt.subplots(figsize=(12, 14))
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 14)
    ax.axis('off')

    # Title
    ax.text(5, 13.5, "Patient Flowchart: NOD2-Stratified Febuxostat Trial",
            ha='center', va='center', fontsize=14, fontweight='bold')

    # Box style
    box_style = dict(boxstyle="round,pad=0.3", facecolor='lightblue', edgecolor='navy', linewidth=2)
    exclude_style = dict(boxstyle="round,pad=0.3", facecolor='lightcoral', edgecolor='darkred', linewidth=2)
    arm_style = dict(boxstyle="round,pad=0.3", facecolor='lightgreen', edgecolor='darkgreen', linewidth=2)

    # Boxes (top to bottom)
    # Screening
    ax.text(5, 12.5, "Screened Patients\n(n=~300)", ha='center', va='center',
            fontsize=11, bbox=box_style)

    # Arrow
    ax.annotate('', xy=(5, 11.5), xytext=(5, 12),
                arrowprops=dict(arrowstyle='->', color='navy', lw=2))

    # Genotyping
    ax.text(5, 11, "NOD2 Genotyping\n(CLIA-certified lab)", ha='center', va='center',
            fontsize=11, bbox=box_style)

    # Split arrows
    ax.annotate('', xy=(2.5, 9.8), xytext=(5, 10.5),
                arrowprops=dict(arrowstyle='->', color='darkred', lw=2))
    ax.annotate('', xy=(5, 9.8), xytext=(5, 10.5),
                arrowprops=dict(arrowstyle='->', color='navy', lw=2))

    # Excluded - L1007fs
    ax.text(2.5, 9.3, "EXCLUDED\nL1007fs carriers\n(n=~30-45)", ha='center', va='center',
            fontsize=10, bbox=exclude_style)
    ax.text(2.5, 8.3, "Pocket deleted:\nGLU1008, ASP1011,\nARG1037 absent", ha='center', va='center',
            fontsize=8, fontstyle='italic', color='darkred')

    # Eligible genotypes
    ax.text(5, 9.3, "Eligible Genotypes\nWT, R702W, G908R,\nN852S, M863V", ha='center', va='center',
            fontsize=10, bbox=box_style)

    # Arrow
    ax.annotate('', xy=(5, 8), xytext=(5, 8.7),
                arrowprops=dict(arrowstyle='->', color='navy', lw=2))

    # Other exclusions
    ax.text(7.5, 8.3, "Other Exclusions:\n- CDAI <220 or >450\n- Severe hepatic/renal\n- Active infection\n(n=~50-70)",
            ha='center', va='center', fontsize=9, bbox=exclude_style)
    ax.annotate('', xy=(7.5, 8.3), xytext=(6, 8.3),
                arrowprops=dict(arrowstyle='->', color='darkred', lw=1.5))

    # Eligible for randomization
    ax.text(5, 7.5, f"Eligible for Randomization\n(n={N_PER_ARM * 3})", ha='center', va='center',
            fontsize=11, bbox=box_style)

    # Arrow to stratification
    ax.annotate('', xy=(5, 6.5), xytext=(5, 7),
                arrowprops=dict(arrowstyle='->', color='navy', lw=2))

    # Stratified randomization box
    ax.text(5, 6, "Stratified Randomization\n1:1:1", ha='center', va='center',
            fontsize=11, bbox=dict(boxstyle="round,pad=0.3", facecolor='lightyellow',
                                   edgecolor='orange', linewidth=2))

    # Stratification factors
    ax.text(8, 6, "Stratification:\n- Genotype (WT vs missense)\n- CDAI (220-300 vs 301-450)\n- Prior biologic (Y/N)",
            ha='left', va='center', fontsize=8, fontstyle='italic')

    # Three arrows to arms
    ax.annotate('', xy=(2.5, 4.5), xytext=(4, 5.5),
                arrowprops=dict(arrowstyle='->', color='darkgreen', lw=2))
    ax.annotate('', xy=(5, 4.5), xytext=(5, 5.5),
                arrowprops=dict(arrowstyle='->', color='darkgreen', lw=2))
    ax.annotate('', xy=(7.5, 4.5), xytext=(6, 5.5),
                arrowprops=dict(arrowstyle='->', color='darkgreen', lw=2))

    # Three arms
    ax.text(2.5, 4, f"Arm 1: Placebo\n(n={N_PER_ARM})", ha='center', va='center',
            fontsize=11, bbox=arm_style)
    ax.text(5, 4, f"Arm 2: Febuxostat\n40mg QD\n(n={N_PER_ARM})\n[Exploratory]", ha='center', va='center',
            fontsize=10, bbox=arm_style)
    ax.text(7.5, 4, f"Arm 3: Febuxostat\n80mg QD\n(n={N_PER_ARM})\n[Primary]", ha='center', va='center',
            fontsize=10, bbox=arm_style)

    # Treatment period
    ax.annotate('', xy=(5, 2.5), xytext=(5, 3.2),
                arrowprops=dict(arrowstyle='->', color='navy', lw=2))
    ax.text(2.5, 2.5, "12-Week Treatment", ha='center', va='center', fontsize=9)
    ax.text(5, 2.5, "12-Week Treatment", ha='center', va='center', fontsize=9)
    ax.text(7.5, 2.5, "12-Week Treatment", ha='center', va='center', fontsize=9)

    # Endpoints
    ax.annotate('', xy=(5, 1.5), xytext=(5, 2),
                arrowprops=dict(arrowstyle='->', color='navy', lw=2))
    ax.text(5, 1, "Week 12: Primary Endpoint\nCDAI-100 Response\n(+4-week follow-up)", ha='center', va='center',
            fontsize=11, bbox=dict(boxstyle="round,pad=0.3", facecolor='gold',
                                   edgecolor='darkorange', linewidth=2))

    # Save
    output_file = os.path.join(OUTPUT_DIR, "PATIENT_FLOWCHART.png")
    plt.tight_layout()
    plt.savefig(output_file, dpi=150, bbox_inches='tight', facecolor='white')
    plt.close()

    print(f"  Created: {output_file}")
    return output_file


def create_trial_timeline():
    """Create Gantt chart for 27-month trial timeline."""

    fig, ax = plt.subplots(figsize=(14, 8))

    # Timeline phases
    phases = [
        ("Site Selection & IRB", 0, 3, 'lightblue'),
        ("Drug Supply & Logistics", 1, 3, 'lightyellow'),
        ("Patient Enrollment", 3, 12, 'lightgreen'),
        ("Treatment (Rolling)", 4, 16, 'lightcoral'),
        ("Follow-up (Rolling)", 5, 18, 'plum'),
        ("Database Lock", 18, 2, 'lightskyblue'),
        ("Statistical Analysis", 19, 3, 'wheat'),
        ("Manuscript Prep", 21, 3, 'lightpink'),
        ("Regulatory Submission", 23, 4, 'lavender'),
    ]

    # Draw bars
    for i, (name, start, duration, color) in enumerate(phases):
        y = len(phases) - i - 1
        ax.barh(y, duration, left=start, height=0.6, color=color, edgecolor='black', linewidth=1)
        ax.text(start + duration/2, y, name, ha='center', va='center', fontsize=9, fontweight='bold')

    # Milestones
    milestones = [
        (3, "First Patient In", 'green'),
        (15, "Last Patient In", 'blue'),
        (18, "Last Patient Out", 'orange'),
        (20, "Database Lock", 'red'),
        (27, "Results", 'purple'),
    ]

    for month, label, color in milestones:
        ax.axvline(x=month, color=color, linestyle='--', linewidth=1.5, alpha=0.7)
        ax.text(month, len(phases) + 0.3, label, ha='center', va='bottom', fontsize=8,
                color=color, rotation=45)

    # Formatting
    ax.set_xlim(-0.5, 28)
    ax.set_ylim(-1.2, len(phases) + 1)
    ax.set_xlabel("Month", fontsize=12)
    ax.set_ylabel("")
    ax.set_title("Phase II Trial Timeline: Febuxostat in Crohn's Disease\n(27-Month Plan)",
                 fontsize=14, fontweight='bold')

    # X-axis ticks
    ax.set_xticks(range(0, 28, 3))
    ax.set_xticklabels([f"M{i}" for i in range(0, 28, 3)])

    # Remove y-axis ticks
    ax.set_yticks([])

    # Grid
    ax.grid(axis='x', alpha=0.3)

    # Legend for milestones
    ax.text(0.5, -0.8, "Milestones: ", fontsize=9, fontweight='bold', transform=ax.transData)
    legend_text = "FPI (M3) | LPI (M15) | LPO (M18) | DBL (M20) | Results (M27)"
    ax.text(3, -0.8, legend_text, fontsize=8, transform=ax.transData)

    # Save
    output_file = os.path.join(OUTPUT_DIR, "TRIAL_TIMELINE.png")
    plt.tight_layout()
    plt.savefig(output_file, dpi=150, bbox_inches='tight', facecolor='white')
    plt.close()

    print(f"  Created: {output_file}")
    return output_file


def main():
    print("=" * 60)
    print("PHASE 12: Clinical Trial Visualization Generator")
    print("=" * 60)

    print("\nGenerating visualizations...")

    # Create flowchart
    print("\n1. Patient Flowchart")
    flowchart_file = create_patient_flowchart()

    # Create timeline
    print("\n2. Trial Timeline (Gantt Chart)")
    timeline_file = create_trial_timeline()

    print("\n" + "=" * 60)
    print("COMPLETE")
    print("=" * 60)
    print(f"\nOutput files:")
    print(f"  - {flowchart_file}")
    print(f"  - {timeline_file}")


if __name__ == "__main__":
    main()
