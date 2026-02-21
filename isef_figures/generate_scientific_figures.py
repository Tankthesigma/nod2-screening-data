"""
ISEF 2026 - NOD2 Project Scientific Figures
Author: Tanmay Vasudeva

REAL scientific publication style - simple, clean, no fancy effects.
"""

import matplotlib.pyplot as plt
import numpy as np

# Standard scientific style
plt.style.use('default')
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = 11
plt.rcParams['axes.linewidth'] = 1
plt.rcParams['axes.spines.top'] = False
plt.rcParams['axes.spines.right'] = False
plt.rcParams['figure.facecolor'] = 'white'
plt.rcParams['axes.facecolor'] = 'white'
plt.rcParams['savefig.facecolor'] = 'white'

OUTPUT_DIR = "C:/Users/vasud/nod2-screening-data/isef_figures/"


def figure1_pipeline():
    """Simple pipeline - just boxes and arrows, no fancy effects."""
    fig, ax = plt.subplots(figsize=(12, 3), dpi=300)
    ax.set_xlim(0, 12)
    ax.set_ylim(0, 3)
    ax.axis('off')

    stages = [
        ('9,566\nCompounds', 'FDA + Natural'),
        ('GNINA\nDocking', 'AlphaFold NOD2'),
        ('ML Scoring', 'AUC 0.85-0.93'),
        ('MD Simulation', '520 ns'),
        ('FEP', 'Gold Standard'),
    ]

    x_positions = [1, 3.2, 5.4, 7.6, 9.8]

    for i, (x, (main, sub)) in enumerate(zip(x_positions, stages)):
        # Simple rectangle
        rect = plt.Rectangle((x-0.8, 0.8), 1.6, 1.4,
                            fill=False, edgecolor='black', linewidth=1.5)
        ax.add_patch(rect)

        # Text
        ax.text(x, 1.7, main, ha='center', va='center', fontsize=10, fontweight='bold')
        ax.text(x, 1.2, sub, ha='center', va='center', fontsize=8, color='gray')

        # Arrow to next
        if i < len(stages) - 1:
            ax.annotate('', xy=(x_positions[i+1]-0.9, 1.5), xytext=(x+0.9, 1.5),
                       arrowprops=dict(arrowstyle='->', color='black', lw=1.5))

    ax.text(6, 2.7, 'Computational Pipeline', ha='center', fontsize=12, fontweight='bold')

    plt.tight_layout()
    plt.savefig(f'{OUTPUT_DIR}figure1_pipeline_scientific.png', dpi=300, bbox_inches='tight')
    plt.close()
    print("Generated: figure1_pipeline_scientific.png")


def figure2_rmsd():
    """Standard bar chart - RMSD comparison."""
    fig, ax = plt.subplots(figsize=(4, 4), dpi=300)

    conditions = ['WT', 'R702W']
    means = [4.05, 3.27]
    sems = [0.15, 0.12]  # Standard error

    x = np.arange(len(conditions))
    bars = ax.bar(x, means, yerr=sems, width=0.6,
                  color=['white', 'white'],
                  edgecolor='black', linewidth=1.5,
                  capsize=5, error_kw={'linewidth': 1.5})

    # Hatch pattern to distinguish
    bars[1].set_hatch('///')

    ax.set_xticks(x)
    ax.set_xticklabels(conditions)
    ax.set_ylabel('Mean RMSD (Å)')
    ax.set_ylim(0, 5.5)

    # Add values above bars
    for i, (m, s) in enumerate(zip(means, sems)):
        ax.text(i, m + s + 0.2, f'{m:.2f}', ha='center', fontsize=10)

    # Sample sizes
    ax.text(0, -0.4, 'n=4', ha='center', fontsize=9, color='gray')
    ax.text(1, -0.4, 'n=3', ha='center', fontsize=9, color='gray')

    ax.set_title('Backbone RMSD (Apo)', fontweight='bold', pad=10)

    plt.tight_layout()
    plt.savefig(f'{OUTPUT_DIR}figure2_rmsd_scientific.png', dpi=300, bbox_inches='tight')
    plt.close()
    print("Generated: figure2_rmsd_scientific.png")


def figure3_heterogeneity():
    """Bar chart - % frames >5A."""
    fig, ax = plt.subplots(figsize=(4, 4), dpi=300)

    conditions = ['WT', 'R702W']
    values = [12.6, 0.0]

    x = np.arange(len(conditions))
    bars = ax.bar(x, values, width=0.6,
                  color=['white', 'white'],
                  edgecolor='black', linewidth=1.5)
    bars[1].set_hatch('///')

    ax.set_xticks(x)
    ax.set_xticklabels(conditions)
    ax.set_ylabel('% Frames with RMSD > 5 Å')
    ax.set_ylim(0, 18)

    # Add values
    ax.text(0, 12.6 + 0.8, '12.6%', ha='center', fontsize=10)
    ax.text(1, 0.8, '0%', ha='center', fontsize=10)

    ax.set_title('Conformational Heterogeneity', fontweight='bold', pad=10)

    plt.tight_layout()
    plt.savefig(f'{OUTPUT_DIR}figure3_heterogeneity_scientific.png', dpi=300, bbox_inches='tight')
    plt.close()
    print("Generated: figure3_heterogeneity_scientific.png")


def figure4_fep():
    """FEP ddG - clean scientific style, NO overlapping text."""
    fig, ax = plt.subplots(figsize=(5, 4), dpi=300)

    compounds = ['Febuxostat', 'Bufadienolide']
    ddg = [2.34, -0.44]
    errors = [0.26, 0.37]

    x = np.arange(len(compounds))
    bars = ax.bar(x, ddg, yerr=errors, width=0.5,
                  color=['white', 'white'],
                  edgecolor='black', linewidth=1.5,
                  capsize=5, error_kw={'linewidth': 1.5})

    bars[0].set_hatch('///')  # Febuxostat hatched

    # Zero line
    ax.axhline(y=0, color='black', linestyle='-', linewidth=0.8)

    ax.set_xticks(x)
    ax.set_xticklabels(compounds)
    ax.set_ylabel('ΔΔG (kcal/mol)')
    ax.set_ylim(-1.5, 3.5)

    # Significance - placed clearly above error bars
    ax.text(0, 2.34 + 0.26 + 0.3, '***', ha='center', fontsize=14, fontweight='bold')
    ax.text(1, -0.44 - 0.37 - 0.2, 'n.s.', ha='center', fontsize=10)

    # Values next to bars (not inside)
    ax.text(0.35, 2.34/2, '+2.34', ha='left', va='center', fontsize=9)
    ax.text(1.35, -0.44/2, '-0.44', ha='left', va='center', fontsize=9)

    ax.set_title('FEP Binding: ΔΔG (R702W − WT)', fontweight='bold', pad=10)

    plt.tight_layout()
    plt.savefig(f'{OUTPUT_DIR}figure4_fep_scientific.png', dpi=300, bbox_inches='tight')
    plt.close()
    print("Generated: figure4_fep_scientific.png")


def figure5_heterogeneity_combined():
    """Combined panel showing WT vs R702W dynamics - simpler."""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7, 3.5), dpi=300)

    # Left: RMSD
    conditions = ['WT', 'R702W']
    means = [4.05, 3.27]
    sems = [0.15, 0.12]

    x = np.arange(2)
    bars1 = ax1.bar(x, means, yerr=sems, width=0.5,
                   color=['white', 'gray'],
                   edgecolor='black', linewidth=1.5,
                   capsize=4)
    ax1.set_xticks(x)
    ax1.set_xticklabels(conditions)
    ax1.set_ylabel('RMSD (Å)')
    ax1.set_ylim(0, 5.5)
    ax1.set_title('A. Mean RMSD', fontweight='bold', loc='left')

    for i, m in enumerate(means):
        ax1.text(i, m + 0.3, f'{m:.2f}', ha='center', fontsize=9)

    # Right: % >5A
    values = [12.6, 0.0]
    bars2 = ax2.bar(x, values, width=0.5,
                   color=['white', 'gray'],
                   edgecolor='black', linewidth=1.5)
    ax2.set_xticks(x)
    ax2.set_xticklabels(conditions)
    ax2.set_ylabel('% Frames > 5 Å')
    ax2.set_ylim(0, 18)
    ax2.set_title('B. Conformational Heterogeneity', fontweight='bold', loc='left')

    ax2.text(0, 12.6 + 0.8, '12.6%', ha='center', fontsize=9)
    ax2.text(1, 0.8, '0%', ha='center', fontsize=9)

    plt.tight_layout()
    plt.savefig(f'{OUTPUT_DIR}figure5_dynamics_scientific.png', dpi=300, bbox_inches='tight')
    plt.close()
    print("Generated: figure5_dynamics_scientific.png")


def figure8_occupancy():
    """Pocket occupancy bar chart."""
    fig, ax = plt.subplots(figsize=(5, 4), dpi=300)

    compounds = ['Febuxostat', 'Bufadienolide', 'Decoy']
    occupancy = [70, 80, 0]

    x = np.arange(len(compounds))
    bars = ax.bar(x, occupancy, width=0.5,
                  color=['white', 'white', 'white'],
                  edgecolor='black', linewidth=1.5)
    bars[0].set_hatch('///')
    bars[1].set_hatch('xxx')

    ax.set_xticks(x)
    ax.set_xticklabels(compounds, fontsize=9)
    ax.set_ylabel('Pocket Occupancy (%)')
    ax.set_ylim(0, 100)

    for i, v in enumerate(occupancy):
        ax.text(i, v + 3, f'{v}%', ha='center', fontsize=10)

    ax.set_title('MD Binding Validation', fontweight='bold', pad=10)

    plt.tight_layout()
    plt.savefig(f'{OUTPUT_DIR}figure8_occupancy_scientific.png', dpi=300, bbox_inches='tight')
    plt.close()
    print("Generated: figure8_occupancy_scientific.png")


if __name__ == "__main__":
    print("Generating scientific-style figures...")
    print("=" * 40)

    figure1_pipeline()
    figure2_rmsd()
    figure3_heterogeneity()
    figure4_fep()
    figure5_heterogeneity_combined()
    figure8_occupancy()

    print("=" * 40)
    print("Done! Figures saved to:", OUTPUT_DIR)
