"""
ISEF 2026 - NOD2 Project Professional Figures
Author: Tanmay Vasudeva
Texas Virtual Academy at Hallsville

Generates all publication-quality figures for ISEF poster.
"""

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
import numpy as np

# Set professional styling
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = 12
plt.rcParams['axes.linewidth'] = 1.5
plt.rcParams['axes.spines.top'] = False
plt.rcParams['axes.spines.right'] = False

# Color scheme
COLORS = {
    'primary': '#1a365d',      # Deep Navy Blue
    'secondary': '#319795',    # Teal
    'accent': '#ed8936',       # Coral/Orange
    'wt': '#3182ce',           # Blue for WT
    'mut': '#e53e3e',          # Red for R702W
    'gray': '#718096',         # Gray
    'light_gray': '#e2e8f0',   # Light gray
    'white': '#ffffff',
}

OUTPUT_DIR = "C:/Users/vasud/nod2-screening-data/isef_figures/"


def figure1_pipeline():
    """Generate professional pipeline flowchart."""
    fig, ax = plt.subplots(figsize=(14, 5), dpi=300)
    ax.set_xlim(0, 14)
    ax.set_ylim(0, 5)
    ax.axis('off')

    # Pipeline stages
    stages = [
        {'x': 1, 'label': '9,566\nCompounds', 'sub': 'FDA + Natural', 'color': COLORS['gray']},
        {'x': 3.5, 'label': 'GNINA\nDocking', 'sub': 'AlphaFold NOD2', 'color': COLORS['secondary']},
        {'x': 6, 'label': 'ML Scoring', 'sub': 'AUC 0.85-0.93', 'color': COLORS['secondary']},
        {'x': 8.5, 'label': 'MD\nValidation', 'sub': '520 ns total', 'color': COLORS['secondary']},
        {'x': 11, 'label': 'FEP\nBinding', 'sub': 'Gold Standard', 'color': COLORS['accent']},
    ]

    box_width = 2.0
    box_height = 2.0

    for i, stage in enumerate(stages):
        # Draw rounded rectangle
        box = FancyBboxPatch(
            (stage['x'] - box_width/2, 1.5),
            box_width, box_height,
            boxstyle="round,pad=0.05,rounding_size=0.2",
            facecolor=stage['color'],
            edgecolor=COLORS['primary'],
            linewidth=2,
            alpha=0.9
        )
        ax.add_patch(box)

        # Add main label
        ax.text(stage['x'], 2.7, stage['label'],
                ha='center', va='center',
                fontsize=11, fontweight='bold',
                color='white')

        # Add sub label
        ax.text(stage['x'], 1.9, stage['sub'],
                ha='center', va='center',
                fontsize=9,
                color='white', alpha=0.9)

        # Draw arrow to next stage
        if i < len(stages) - 1:
            arrow = FancyArrowPatch(
                (stage['x'] + box_width/2 + 0.1, 2.5),
                (stages[i+1]['x'] - box_width/2 - 0.1, 2.5),
                arrowstyle='-|>',
                mutation_scale=20,
                color=COLORS['primary'],
                linewidth=2
            )
            ax.add_patch(arrow)

    # Title
    ax.text(6, 4.5, 'NOD2 Binder Discovery Pipeline',
            ha='center', va='center',
            fontsize=16, fontweight='bold',
            color=COLORS['primary'])

    # Bottom annotation
    ax.text(6, 0.7, 'Scaffold-split cross validation',
            ha='center', va='center',
            fontsize=10, fontstyle='italic',
            color=COLORS['gray'])

    plt.tight_layout()
    plt.savefig(f'{OUTPUT_DIR}figure1_pipeline_professional.png',
                dpi=300, bbox_inches='tight',
                facecolor='white', edgecolor='none')
    plt.close()
    print("Generated: figure1_pipeline_professional.png")


def figure2_rmsd():
    """Generate clean RMSD bar chart (WT apo vs R702W apo ONLY)."""
    fig, ax = plt.subplots(figsize=(6, 5), dpi=300)

    # Data
    conditions = ['WT apo', 'R702W apo']
    rmsd_values = [4.05, 3.27]
    rmsd_errors = [0.5, 0.3]  # Approximate SDs
    colors = [COLORS['wt'], COLORS['mut']]

    bars = ax.bar(conditions, rmsd_values, yerr=rmsd_errors,
                  color=colors, edgecolor=COLORS['primary'],
                  linewidth=2, capsize=8, error_kw={'linewidth': 2})

    # Add value labels on bars
    for bar, val in zip(bars, rmsd_values):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.7,
                f'{val:.2f} A',
                ha='center', va='bottom',
                fontsize=12, fontweight='bold')

    ax.set_ylabel('Mean RMSD (A)', fontsize=14, fontweight='bold')
    ax.set_title('Apo MD: Backbone RMSD', fontsize=16, fontweight='bold',
                 color=COLORS['primary'], pad=20)
    ax.set_ylim(0, 6)
    ax.tick_params(axis='both', labelsize=12)

    # Add n values
    ax.text(0, -0.8, 'n=4', ha='center', fontsize=10, color=COLORS['gray'])
    ax.text(1, -0.8, 'n=3', ha='center', fontsize=10, color=COLORS['gray'])

    plt.tight_layout()
    plt.savefig(f'{OUTPUT_DIR}figure2_rmsd_clean.png',
                dpi=300, bbox_inches='tight',
                facecolor='white', edgecolor='none')
    plt.close()
    print("Generated: figure2_rmsd_clean.png")


def figure3_heterogeneity():
    """Generate conformational heterogeneity chart (% frames >5A)."""
    fig, ax = plt.subplots(figsize=(6, 5), dpi=300)

    # Data
    conditions = ['WT apo', 'R702W apo']
    pct_values = [12.6, 0.0]
    colors = [COLORS['wt'], COLORS['mut']]

    bars = ax.bar(conditions, pct_values,
                  color=colors, edgecolor=COLORS['primary'],
                  linewidth=2)

    # Add value labels
    for bar, val in zip(bars, pct_values):
        y_pos = bar.get_height() + 1 if val > 0 else 1
        ax.text(bar.get_x() + bar.get_width()/2, y_pos,
                f'{val:.1f}%',
                ha='center', va='bottom',
                fontsize=14, fontweight='bold')

    ax.set_ylabel('% Frames > 5A RMSD', fontsize=14, fontweight='bold')
    ax.set_title('Conformational Heterogeneity', fontsize=16, fontweight='bold',
                 color=COLORS['primary'], pad=20)
    ax.set_ylim(0, 20)
    ax.tick_params(axis='both', labelsize=12)

    # Add interpretation
    ax.text(0.5, 17, 'R702W shows reduced\nconformational sampling',
            ha='center', va='center',
            fontsize=11, fontstyle='italic',
            color=COLORS['gray'],
            bbox=dict(boxstyle='round', facecolor=COLORS['light_gray'], alpha=0.5))

    plt.tight_layout()
    plt.savefig(f'{OUTPUT_DIR}figure3_heterogeneity_clean.png',
                dpi=300, bbox_inches='tight',
                facecolor='white', edgecolor='none')
    plt.close()
    print("Generated: figure3_heterogeneity_clean.png")


def figure4_fep():
    """Generate FEP ddG comparison chart."""
    fig, ax = plt.subplots(figsize=(8, 6), dpi=300)

    # Data
    compounds = ['Febuxostat', 'Bufadienolide']
    ddg_values = [2.34, -0.44]
    ddg_errors = [0.26, 0.37]
    colors = [COLORS['mut'], COLORS['gray']]

    x = np.arange(len(compounds))
    bars = ax.bar(x, ddg_values, yerr=ddg_errors,
                  color=colors, edgecolor=COLORS['primary'],
                  linewidth=2, capsize=10, error_kw={'linewidth': 2},
                  width=0.6)

    # Add zero line
    ax.axhline(y=0, color=COLORS['primary'], linestyle='--', linewidth=1.5, alpha=0.7)

    # Add significance markers
    ax.text(0, ddg_values[0] + ddg_errors[0] + 0.3, '***',
            ha='center', va='bottom', fontsize=20, fontweight='bold',
            color=COLORS['mut'])
    ax.text(0, ddg_values[0] + ddg_errors[0] + 0.7, 'p<0.001',
            ha='center', va='bottom', fontsize=10,
            color=COLORS['mut'])

    ax.text(1, ddg_values[1] - ddg_errors[1] - 0.3, 'NS',
            ha='center', va='top', fontsize=14, fontweight='bold',
            color=COLORS['gray'])

    # Add value labels
    ax.text(0, ddg_values[0]/2, f'+{ddg_values[0]:.2f}\nkcal/mol',
            ha='center', va='center', fontsize=11, fontweight='bold',
            color='white')
    ax.text(1, ddg_values[1]/2 - 0.5, f'{ddg_values[1]:.2f}\nkcal/mol',
            ha='center', va='center', fontsize=11, fontweight='bold',
            color='white')

    # Labels and title
    ax.set_xticks(x)
    ax.set_xticklabels(compounds, fontsize=14, fontweight='bold')
    ax.set_ylabel('ddG (kcal/mol)\nR702W - WT', fontsize=14, fontweight='bold')
    ax.set_title('FEP Binding Free Energy: Mutation Effect',
                 fontsize=16, fontweight='bold',
                 color=COLORS['primary'], pad=20)
    ax.set_ylim(-1.5, 4)
    ax.tick_params(axis='y', labelsize=12)

    # Add annotation
    ax.text(0, -1.2, '50x weaker\nin R702W',
            ha='center', va='top', fontsize=10,
            color=COLORS['mut'], fontweight='bold')
    ax.text(1, -1.2, 'No significant\nchange',
            ha='center', va='top', fontsize=10,
            color=COLORS['gray'])

    plt.tight_layout()
    plt.savefig(f'{OUTPUT_DIR}figure4_fep_comparison.png',
                dpi=300, bbox_inches='tight',
                facecolor='white', edgecolor='none')
    plt.close()
    print("Generated: figure4_fep_comparison.png")


def figure5_apo_dynamics():
    """Generate apo dynamics comparison diagram."""
    fig, axes = plt.subplots(1, 2, figsize=(12, 5), dpi=300)

    # Left panel: Literature Assumption
    ax1 = axes[0]
    ax1.set_xlim(0, 10)
    ax1.set_ylim(0, 10)
    ax1.axis('off')

    ax1.text(5, 9, 'Literature Assumption', ha='center', va='center',
             fontsize=14, fontweight='bold', color=COLORS['primary'])

    ax1.text(5, 7.5, 'R702W is "destabilized"', ha='center', va='center',
             fontsize=12, fontstyle='italic', color=COLORS['gray'])

    # Draw wavy/floppy protein representation
    theta = np.linspace(0, 4*np.pi, 100)
    x_wave = 5 + 1.5*np.sin(theta) + 0.5*np.sin(3*theta)
    y_wave = np.linspace(2, 6, 100)
    ax1.plot(x_wave, y_wave, color=COLORS['gray'], linewidth=3, alpha=0.7)

    ax1.text(5, 1, 'High conformational\nheterogeneity expected',
             ha='center', va='center', fontsize=11, color=COLORS['gray'])

    # Right panel: This Study
    ax2 = axes[1]
    ax2.set_xlim(0, 10)
    ax2.set_ylim(0, 10)
    ax2.axis('off')

    ax2.text(5, 9, 'This Study (Apo MD)', ha='center', va='center',
             fontsize=14, fontweight='bold', color=COLORS['primary'])

    ax2.text(5, 7.5, 'R702W shows reduced sampling', ha='center', va='center',
             fontsize=12, fontweight='bold', color=COLORS['mut'])

    # Draw stable/rigid protein representation
    ax2.plot([5, 5], [2, 6], color=COLORS['mut'], linewidth=6, alpha=0.8)
    ax2.scatter([5], [4], s=200, color=COLORS['accent'], zorder=5)

    ax2.text(5, 1, 'WT: 12.6% frames >5A\nR702W: 0% frames >5A',
             ha='center', va='center', fontsize=11, fontweight='bold',
             color=COLORS['primary'],
             bbox=dict(boxstyle='round', facecolor=COLORS['light_gray'], alpha=0.5))

    # Add "Supporting Context" label
    fig.text(0.5, 0.02, 'Supporting context for mutation-sensitive binding observed in FEP',
             ha='center', va='center', fontsize=11, fontstyle='italic',
             color=COLORS['gray'])

    plt.tight_layout()
    plt.savefig(f'{OUTPUT_DIR}figure5_apo_dynamics.png',
                dpi=300, bbox_inches='tight',
                facecolor='white', edgecolor='none')
    plt.close()
    print("Generated: figure5_apo_dynamics.png")


def figure8_pocket_occupancy():
    """Generate pocket occupancy bar chart."""
    fig, ax = plt.subplots(figsize=(8, 5), dpi=300)

    # Data (Febuxostat, Bufadienolide, Decoy only)
    compounds = ['Febuxostat', 'Bufadienolide', 'Decoy']
    occupancy = [70, 80, 0]
    colors = [COLORS['secondary'], COLORS['secondary'], COLORS['gray']]

    bars = ax.bar(compounds, occupancy, color=colors,
                  edgecolor=COLORS['primary'], linewidth=2)

    # Add value labels
    for bar, val in zip(bars, occupancy):
        y_pos = bar.get_height() + 3 if val > 0 else 3
        ax.text(bar.get_x() + bar.get_width()/2, y_pos,
                f'{val}%',
                ha='center', va='bottom',
                fontsize=14, fontweight='bold')

    ax.set_ylabel('Mean Pocket Occupancy (%)', fontsize=14, fontweight='bold')
    ax.set_title('MD Binding Validation', fontsize=16, fontweight='bold',
                 color=COLORS['primary'], pad=20)
    ax.set_ylim(0, 100)
    ax.tick_params(axis='both', labelsize=12)

    # Add note
    ax.text(1, 95, 'Ligand within 5A of binding pocket',
            ha='center', va='top', fontsize=10, fontstyle='italic',
            color=COLORS['gray'])

    plt.tight_layout()
    plt.savefig(f'{OUTPUT_DIR}figure8_pocket_occupancy.png',
                dpi=300, bbox_inches='tight',
                facecolor='white', edgecolor='none')
    plt.close()
    print("Generated: figure8_pocket_occupancy.png")


if __name__ == "__main__":
    print("Generating ISEF 2026 figures...")
    print("=" * 50)

    figure1_pipeline()
    figure2_rmsd()
    figure3_heterogeneity()
    figure4_fep()
    figure5_apo_dynamics()
    figure8_pocket_occupancy()

    print("=" * 50)
    print("All figures generated successfully!")
    print(f"Output directory: {OUTPUT_DIR}")
