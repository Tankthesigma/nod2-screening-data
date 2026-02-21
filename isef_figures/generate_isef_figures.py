#!/usr/bin/env python3
"""
ISEF Figure Generator for NOD2-CROHN Drug Discovery Project
Generates publication-quality figures for ISEF 2026

Author: Tanmay Vasudeva
"""

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
import numpy as np
from pathlib import Path

# Configuration
OUTPUT_DIR = Path(__file__).parent
DPI = 300

# Set publication style
plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
    'font.size': 12,
    'axes.labelsize': 14,
    'axes.titlesize': 16,
    'axes.titleweight': 'bold',
    'xtick.labelsize': 11,
    'ytick.labelsize': 11,
    'legend.fontsize': 10,
    'figure.facecolor': 'white',
    'axes.facecolor': 'white',
    'savefig.facecolor': 'white',
    'savefig.edgecolor': 'white',
    'axes.spines.top': False,
    'axes.spines.right': False,
})


def create_figure1_pipeline():
    """Create computational pipeline flowchart."""
    fig, ax = plt.subplots(figsize=(14, 8))
    ax.set_xlim(0, 14)
    ax.set_ylim(0, 8)
    ax.axis('off')

    # Color scheme
    colors = {
        'input': '#E3F2FD',      # Light blue
        'docking': '#BBDEFB',    # Blue
        'ml': '#90CAF9',         # Medium blue
        'md': '#64B5F6',         # Darker blue
        'fep': '#42A5F5',        # Dark blue
        'output': '#1E88E5',     # Darkest blue
        'arrow': '#424242',
        'text': '#212121',
    }

    # Box positions and sizes
    box_height = 1.2
    box_y = 4.0

    boxes = [
        {'x': 0.5, 'w': 2.2, 'label': 'Compound\nLibrary', 'sub': '9,566 compounds', 'color': colors['input']},
        {'x': 3.2, 'w': 2.2, 'label': 'GNINA\nDocking', 'sub': 'CNN scoring', 'color': colors['docking']},
        {'x': 5.9, 'w': 2.2, 'label': 'NOD2-Scout\nML Model', 'sub': 'AUC = 0.90', 'color': colors['ml']},
        {'x': 8.6, 'w': 2.2, 'label': 'MD\nSimulation', 'sub': '520 ns total', 'color': colors['md']},
        {'x': 11.3, 'w': 2.2, 'label': 'FEP\nValidation', 'sub': 'Binding ΔΔG', 'color': colors['fep']},
    ]

    # Draw boxes
    for box in boxes:
        rect = FancyBboxPatch(
            (box['x'], box_y - box_height/2), box['w'], box_height,
            boxstyle="round,pad=0.05,rounding_size=0.2",
            facecolor=box['color'], edgecolor='#1565C0', linewidth=2
        )
        ax.add_patch(rect)
        ax.text(box['x'] + box['w']/2, box_y + 0.15, box['label'],
                ha='center', va='center', fontsize=11, fontweight='bold', color=colors['text'])
        ax.text(box['x'] + box['w']/2, box_y - 0.35, box['sub'],
                ha='center', va='center', fontsize=9, color='#616161', style='italic')

    # Draw arrows
    for i in range(len(boxes) - 1):
        start_x = boxes[i]['x'] + boxes[i]['w']
        end_x = boxes[i+1]['x']
        ax.annotate('', xy=(end_x - 0.1, box_y), xytext=(start_x + 0.1, box_y),
                   arrowprops=dict(arrowstyle='->', color=colors['arrow'], lw=2))

    # Add filtering annotations (numbers below arrows)
    filter_texts = [
        {'x': 3.9, 'text': '→ 2,129'},
        {'x': 6.6, 'text': '→ 144'},
        {'x': 9.3, 'text': '→ 8 Tier-1'},
        {'x': 12.0, 'text': '→ Lead'},
    ]
    for ft in filter_texts:
        ax.text(ft['x'], box_y - 1.0, ft['text'], ha='center', va='center',
                fontsize=10, color='#C62828', fontweight='bold')

    # Add output box at bottom
    output_rect = FancyBboxPatch(
        (5, 1.0), 4, 1.0,
        boxstyle="round,pad=0.05,rounding_size=0.2",
        facecolor='#C8E6C9', edgecolor='#2E7D32', linewidth=2
    )
    ax.add_patch(output_rect)
    ax.text(7, 1.5, 'FEBUXOSTAT', ha='center', va='center',
            fontsize=14, fontweight='bold', color='#1B5E20')
    ax.text(7, 1.15, 'Lead Drug Candidate', ha='center', va='center',
            fontsize=10, color='#388E3C')

    # Arrow from FEP to output
    ax.annotate('', xy=(7, 2.0), xytext=(12.4, box_y - box_height/2 - 0.1),
               arrowprops=dict(arrowstyle='->', color='#2E7D32', lw=2,
                              connectionstyle="arc3,rad=-0.3"))

    # Title
    ax.text(7, 7.2, 'Computational Drug Discovery Pipeline for NOD2-Crohn\'s Disease',
            ha='center', va='center', fontsize=16, fontweight='bold', color='#0D47A1')

    # Legend/key
    ax.text(0.5, 0.5, 'Screening Funnel: 9,566 → 2,129 → 144 → 8 → 1 compound',
            fontsize=10, color='#616161', style='italic')

    plt.tight_layout()
    plt.savefig(OUTPUT_DIR / 'figure1_pipeline.png', dpi=DPI, bbox_inches='tight',
                facecolor='white', edgecolor='none')
    plt.close()
    print("Saved: figure1_pipeline.png")


def create_figure2_rmsd_comparison():
    """Create RMSD comparison bar chart."""
    fig, ax = plt.subplots(figsize=(8, 6))

    # Data
    groups = ['WT apo', 'R702W apo']
    means = [4.05, 3.27]
    stds = [0.69, 0.21]
    colors = ['#1976D2', '#D32F2F']

    x = np.arange(len(groups))
    bars = ax.bar(x, means, yerr=stds, capsize=8, color=colors,
                  edgecolor='black', linewidth=1.5, alpha=0.85, width=0.6)

    # Add value labels on bars
    for bar, mean, std in zip(bars, means, stds):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + std + 0.15,
                f'{mean:.2f} Å', ha='center', va='bottom', fontsize=12, fontweight='bold')

    # Add significance bracket
    ax.plot([0, 0, 1, 1], [5.2, 5.4, 5.4, 5.2], 'k-', lw=1.5)
    ax.text(0.5, 5.5, 'p = 0.17', ha='center', va='bottom', fontsize=10)

    # Annotations
    ax.annotate('Higher flexibility\n(12.6% >5Å)', xy=(0, 4.05), xytext=(-0.7, 5.5),
                fontsize=9, ha='center',
                arrowprops=dict(arrowstyle='->', color='#1976D2', lw=1.5))
    ax.annotate('Restricted\n(0% >5Å)', xy=(1, 3.27), xytext=(1.7, 4.5),
                fontsize=9, ha='center',
                arrowprops=dict(arrowstyle='->', color='#D32F2F', lw=1.5))

    ax.set_xticks(x)
    ax.set_xticklabels(groups, fontsize=13)
    ax.set_ylabel('Mean RMSD (Å)', fontsize=14)
    ax.set_title('R702W Mutation Restricts NOD2 Conformational Dynamics',
                 fontsize=14, fontweight='bold', pad=15)
    ax.set_ylim(0, 6.5)
    ax.axhline(5.0, color='red', linestyle='--', alpha=0.5, label='Extended state threshold')
    ax.legend(loc='upper right', fontsize=9)

    # Add sample sizes
    ax.text(0, -0.4, 'n=4', ha='center', fontsize=10, color='gray')
    ax.text(1, -0.4, 'n=3', ha='center', fontsize=10, color='gray')

    plt.tight_layout()
    plt.savefig(OUTPUT_DIR / 'figure2_rmsd_comparison.png', dpi=DPI, bbox_inches='tight',
                facecolor='white', edgecolor='none')
    plt.close()
    print("Saved: figure2_rmsd_comparison.png")


def create_figure3_heterogeneity():
    """Create conformational heterogeneity bar chart."""
    fig, ax = plt.subplots(figsize=(10, 6))

    # Data for all groups
    groups = ['WT\napo', 'R702W\napo', 'WT\n+Febuxostat', 'R702W\n+Febuxostat',
              'WT\n+Budesonide', 'WT\n+Bufadienolide', 'R702W\n+Bufadienolide']
    fracs = [12.6, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    colors = ['#1976D2', '#D32F2F', '#388E3C', '#FF8F00', '#7B1FA2', '#795548', '#E91E63']

    x = np.arange(len(groups))
    bars = ax.bar(x, fracs, color=colors, edgecolor='black', linewidth=1.5, alpha=0.85, width=0.7)

    # Highlight WT apo
    bars[0].set_edgecolor('#0D47A1')
    bars[0].set_linewidth(3)

    # Add value labels
    for bar, frac in zip(bars, fracs):
        height = max(frac, 0.3)  # Minimum height for visibility
        ax.text(bar.get_x() + bar.get_width()/2, height + 0.5,
                f'{frac:.1f}%', ha='center', va='bottom', fontsize=11, fontweight='bold')

    ax.set_xticks(x)
    ax.set_xticklabels(groups, fontsize=10)
    ax.set_ylabel('% Frames with RMSD > 5Å', fontsize=14)
    ax.set_title('R702W Mutation Eliminates Extended Conformational States',
                 fontsize=14, fontweight='bold', pad=15)
    ax.set_ylim(0, 18)

    # Add annotation
    ax.annotate('Only WT apo samples\nextended HD2 states',
                xy=(0, 12.6), xytext=(2.5, 15),
                fontsize=11, ha='center',
                bbox=dict(boxstyle='round,pad=0.3', facecolor='#FFECB3', edgecolor='#FF8F00'),
                arrowprops=dict(arrowstyle='->', color='#FF8F00', lw=2))

    # Add grouping brackets
    ax.axhline(y=-1.5, xmin=0.02, xmax=0.13, color='blue', linewidth=2, clip_on=False)
    ax.axhline(y=-1.5, xmin=0.16, xmax=0.84, color='green', linewidth=2, clip_on=False)

    plt.tight_layout()
    plt.savefig(OUTPUT_DIR / 'figure3_heterogeneity.png', dpi=DPI, bbox_inches='tight',
                facecolor='white', edgecolor='none')
    plt.close()
    print("Saved: figure3_heterogeneity.png")


def create_figure4_fep_results():
    """Create FEP binding energy comparison."""
    fig, ax = plt.subplots(figsize=(9, 6))

    # Data from MASTER_PROJECT_DATA
    compounds = ['Febuxostat', 'Bufadienolide\n(CID_10120)']
    wt_binding = [-22.5, -25.6]  # WT binding energy
    mut_binding = [-25.2, -29.5]  # R702W binding energy

    x = np.arange(len(compounds))
    width = 0.35

    bars1 = ax.bar(x - width/2, wt_binding, width, label='WT NOD2',
                   color='#1976D2', edgecolor='black', linewidth=1.5)
    bars2 = ax.bar(x + width/2, mut_binding, width, label='R702W NOD2',
                   color='#D32F2F', edgecolor='black', linewidth=1.5)

    # Add value labels
    for bars in [bars1, bars2]:
        for bar in bars:
            height = bar.get_height()
            ax.text(bar.get_x() + bar.get_width()/2, height - 1.5,
                    f'{height:.1f}', ha='center', va='top', fontsize=10,
                    fontweight='bold', color='white')

    # Add DDE annotations
    ax.annotate('ΔΔG = -2.7 kcal/mol\n(R702W binds STRONGER)',
                xy=(0, -25.2), xytext=(0.8, -15),
                fontsize=10, ha='center',
                bbox=dict(boxstyle='round,pad=0.3', facecolor='#FFCDD2', edgecolor='#D32F2F'),
                arrowprops=dict(arrowstyle='->', color='#D32F2F', lw=1.5))

    ax.set_xticks(x)
    ax.set_xticklabels(compounds, fontsize=12)
    ax.set_ylabel('Binding Energy (kcal/mol)', fontsize=14)
    ax.set_title('MM-GBSA Binding Energy Analysis:\nMutation Enhances Drug Binding',
                 fontsize=14, fontweight='bold', pad=15)
    ax.legend(loc='lower right', fontsize=11)
    ax.set_ylim(-35, 0)

    # Add "more negative = stronger" note
    ax.text(0.02, 0.02, '↓ More negative = Stronger binding',
            transform=ax.transAxes, fontsize=9, color='gray', style='italic')

    plt.tight_layout()
    plt.savefig(OUTPUT_DIR / 'figure4_fep_binding.png', dpi=DPI, bbox_inches='tight',
                facecolor='white', edgecolor='none')
    plt.close()
    print("Saved: figure4_fep_binding.png")


def create_figure5_drug_effects():
    """Create drug stabilization effects comparison."""
    fig, ax = plt.subplots(figsize=(10, 6))

    # Data: change in RMSD from apo state
    drugs = ['Febuxostat\n(on WT)', 'Febuxostat\n(on R702W)',
             'Budesonide\n(on WT)', 'Bufadienolide\n(on WT)', 'Bufadienolide\n(on R702W)']
    effects = [-0.66, -0.23, -0.34, -0.38, +0.09]
    colors = ['#388E3C' if e < 0 else '#D32F2F' for e in effects]

    x = np.arange(len(drugs))
    bars = ax.bar(x, effects, color=colors, edgecolor='black', linewidth=1.5, alpha=0.85, width=0.65)

    # Add value labels
    for bar, effect in zip(bars, effects):
        va = 'bottom' if effect >= 0 else 'top'
        offset = 0.03 if effect >= 0 else -0.03
        ax.text(bar.get_x() + bar.get_width()/2, effect + offset,
                f'{effect:+.2f} Å', ha='center', va=va, fontsize=11, fontweight='bold')

    ax.axhline(0, color='black', linewidth=1)
    ax.set_xticks(x)
    ax.set_xticklabels(drugs, fontsize=10)
    ax.set_ylabel('Change in Mean RMSD (Å)', fontsize=14)
    ax.set_title('Drug Stabilization Effects on NOD2 Conformational Dynamics',
                 fontsize=14, fontweight='bold', pad=15)
    ax.set_ylim(-0.9, 0.3)

    # Add legend
    green_patch = mpatches.Patch(color='#388E3C', label='Stabilizing (reduces RMSD)')
    red_patch = mpatches.Patch(color='#D32F2F', label='Destabilizing (increases RMSD)')
    ax.legend(handles=[green_patch, red_patch], loc='lower right', fontsize=10)

    # Highlight best result
    ax.annotate('Best\nstabilization', xy=(0, -0.66), xytext=(0, -0.85),
                fontsize=9, ha='center', fontweight='bold', color='#1B5E20')

    plt.tight_layout()
    plt.savefig(OUTPUT_DIR / 'figure5_drug_effects.png', dpi=DPI, bbox_inches='tight',
                facecolor='white', edgecolor='none')
    plt.close()
    print("Saved: figure5_drug_effects.png")


def create_figure6_precision_medicine():
    """Create precision medicine stratification figure."""
    fig, ax = plt.subplots(figsize=(11, 7))
    ax.axis('off')

    # Title
    ax.text(0.5, 0.95, 'Precision Medicine Stratification for NOD2-Crohn\'s Patients',
            ha='center', va='top', fontsize=16, fontweight='bold', transform=ax.transAxes)

    # Data for patient groups
    groups = [
        {'name': 'Wild-Type NOD2', 'prevalence': '~50%', 'recommendation': 'Febuxostat Candidate',
         'color': '#1976D2', 'evidence': 'Baseline binding', 'x': 0.15},
        {'name': 'R702W Mutation', 'prevalence': '~10%', 'recommendation': 'ENHANCED Response Expected',
         'color': '#388E3C', 'evidence': 'ΔΔG = -2.7 kcal/mol', 'x': 0.40},
        {'name': 'G908R Mutation', 'prevalence': '~5%', 'recommendation': 'Febuxostat Candidate',
         'color': '#FF8F00', 'evidence': 'Similar to WT', 'x': 0.65},
        {'name': 'L1007fs Mutation', 'prevalence': '~35%', 'recommendation': 'EXCLUDED',
         'color': '#D32F2F', 'evidence': 'Binding pocket deleted', 'x': 0.90},
    ]

    for g in groups:
        # Box
        rect = FancyBboxPatch(
            (g['x'] - 0.10, 0.35), 0.20, 0.45,
            boxstyle="round,pad=0.02,rounding_size=0.02",
            facecolor=g['color'], edgecolor='black', linewidth=2,
            alpha=0.2, transform=ax.transAxes
        )
        ax.add_patch(rect)

        # Header
        ax.text(g['x'], 0.75, g['name'], ha='center', va='center',
                fontsize=11, fontweight='bold', transform=ax.transAxes)
        ax.text(g['x'], 0.68, g['prevalence'], ha='center', va='center',
                fontsize=10, color='gray', transform=ax.transAxes)

        # Recommendation
        rec_color = '#1B5E20' if 'Candidate' in g['recommendation'] or 'ENHANCED' in g['recommendation'] else '#B71C1C'
        ax.text(g['x'], 0.55, g['recommendation'], ha='center', va='center',
                fontsize=10, fontweight='bold', color=rec_color, transform=ax.transAxes)

        # Evidence
        ax.text(g['x'], 0.42, g['evidence'], ha='center', va='center',
                fontsize=9, style='italic', color='#616161', transform=ax.transAxes)

    # Bottom summary
    ax.text(0.5, 0.15, '~65% of NOD2-Crohn\'s patients eligible for Febuxostat therapy',
            ha='center', va='center', fontsize=13, fontweight='bold',
            color='#0D47A1', transform=ax.transAxes,
            bbox=dict(boxstyle='round,pad=0.5', facecolor='#E3F2FD', edgecolor='#1976D2'))

    ax.text(0.5, 0.05, 'L1007fs patients (~35%) require alternative therapy - binding pocket is deleted',
            ha='center', va='center', fontsize=10, color='#B71C1C',
            style='italic', transform=ax.transAxes)

    plt.tight_layout()
    plt.savefig(OUTPUT_DIR / 'figure6_precision_medicine.png', dpi=DPI, bbox_inches='tight',
                facecolor='white', edgecolor='none')
    plt.close()
    print("Saved: figure6_precision_medicine.png")


def create_figure7_mechanism_summary():
    """Create mechanism summary schematic."""
    fig, axes = plt.subplots(1, 3, figsize=(14, 5))

    titles = ['A. Wild-Type NOD2', 'B. R702W Mutant', 'C. R702W + Febuxostat']
    colors = ['#1976D2', '#D32F2F', '#388E3C']

    for ax, title, color in zip(axes, titles, colors):
        ax.set_xlim(0, 10)
        ax.set_ylim(0, 10)
        ax.axis('off')
        ax.set_title(title, fontsize=13, fontweight='bold', color=color, pad=10)

        # Draw protein outline
        protein = plt.Circle((5, 5), 3.5, fill=False, color=color, linewidth=3)
        ax.add_patch(protein)

        if 'Wild-Type' in title:
            # Show extended state sampling (multiple conformations)
            for angle in [0, 45, 90, 135]:
                rad = np.radians(angle)
                x_end = 5 + 4.5 * np.cos(rad)
                y_end = 5 + 4.5 * np.sin(rad)
                ax.annotate('', xy=(x_end, y_end), xytext=(5, 5),
                           arrowprops=dict(arrowstyle='->', color=color, alpha=0.4, lw=1.5))
            ax.text(5, 1, 'Samples extended\nHD2 states (12.6%)',
                    ha='center', fontsize=10, color=color)

        elif 'R702W Mutant' in title:
            # Show restricted conformation
            ax.add_patch(plt.Circle((5, 5), 2.5, fill=True, color=color, alpha=0.3))
            ax.text(5, 5, 'RESTRICTED', ha='center', va='center',
                    fontsize=11, fontweight='bold', color='white')
            ax.text(5, 1, 'Confined to\nclosed state (0%)',
                    ha='center', fontsize=10, color=color)

        else:  # Drug bound
            # Show drug stabilization
            ax.add_patch(plt.Circle((5, 5), 2.5, fill=True, color=color, alpha=0.3))
            # Drug molecule
            drug = plt.Rectangle((4.2, 4.2), 1.6, 1.6, fill=True,
                                  color='#FFD700', edgecolor='black', linewidth=2)
            ax.add_patch(drug)
            ax.text(5, 5, 'FEB', ha='center', va='center',
                    fontsize=9, fontweight='bold')
            ax.text(5, 1, 'Drug stabilizes\nprotein (-0.23 Å)',
                    ha='center', fontsize=10, color=color)

    plt.suptitle('Mechanism: R702W Restricts Dynamics, Febuxostat Provides Additional Stabilization',
                 fontsize=14, fontweight='bold', y=1.02)

    plt.tight_layout()
    plt.savefig(OUTPUT_DIR / 'figure7_mechanism_summary.png', dpi=DPI, bbox_inches='tight',
                facecolor='white', edgecolor='none')
    plt.close()
    print("Saved: figure7_mechanism_summary.png")


def create_pymol_script():
    """Generate PyMOL visualization script."""
    script = '''# PyMOL Visualization Script for NOD2-Febuxostat Complex
# ISEF 2026 - NOD2-CROHN Drug Discovery Project
# Author: Tanmay Vasudeva

# Load structures
load nod2_alphafold.pdb, NOD2
load febuxostat_docked.pdb, febuxostat

# Basic setup
bg_color white
set ray_opaque_background, 1
set antialias, 2
set ray_trace_mode, 1

# Color scheme
color marine, NOD2
color gray80, NOD2 and name C*

# Highlight domains
# HD2 domain (contains R702 mutation site)
select HD2, resi 539-744
color lightblue, HD2

# LRR domain (drug binding region)
select LRR, resi 745-1040
color palegreen, LRR

# Mutation site R702
select mut_site, resi 702
show spheres, mut_site
color red, mut_site

# Binding pocket residues
select pocket, resi 1007+1008+1009+1010+1011+1012+1013+1014+1035+1037
show sticks, pocket
color yellow, pocket and name C*

# Febuxostat ligand
show sticks, febuxostat
color orange, febuxostat and name C*
set stick_radius, 0.2, febuxostat

# Surface for binding pocket
select pocket_surface, resi 1004-1040
show surface, pocket_surface
set transparency, 0.7, pocket_surface
color white, pocket_surface

# Labels
# Uncomment to add labels
# label resi 702 and name CA, "R702W"
# label resi 1008 and name CA, "E1008"

# Camera settings for different views
# View 1: Overall structure
orient
zoom NOD2, 10
ray 2400, 2400
png figure_nod2_overall.png, dpi=300

# View 2: Binding pocket close-up
zoom pocket, 15
turn y, 30
ray 2400, 2400
png figure_binding_pocket.png, dpi=300

# View 3: Mutation site
zoom mut_site, 20
turn x, 15
ray 2400, 2400
png figure_mutation_site.png, dpi=300

# Session save
save nod2_febuxostat_session.pse

print("PyMOL visualization complete!")
print("Generated: figure_nod2_overall.png")
print("Generated: figure_binding_pocket.png")
print("Generated: figure_mutation_site.png")
print("Session saved: nod2_febuxostat_session.pse")
'''

    with open(OUTPUT_DIR / 'pymol_visualization.pml', 'w') as f:
        f.write(script)
    print("Saved: pymol_visualization.pml")


def create_figure8_clinical_trial():
    """Create clinical trial design summary."""
    fig, ax = plt.subplots(figsize=(10, 7))
    ax.axis('off')

    # Title
    ax.text(0.5, 0.95, 'Phase 2a Clinical Trial Design', ha='center', va='top',
            fontsize=16, fontweight='bold', transform=ax.transAxes)
    ax.text(0.5, 0.89, 'Febuxostat for NOD2-Associated Crohn\'s Disease',
            ha='center', va='top', fontsize=12, color='gray', transform=ax.transAxes)

    # Trial parameters box
    params = [
        ('Design', 'Randomized, Double-Blind, Placebo-Controlled'),
        ('Arms', '3 (Placebo, 40mg, 80mg Febuxostat)'),
        ('Duration', '12 weeks treatment + 4 weeks follow-up'),
        ('Sample Size', '210 patients (70 per arm)'),
        ('Power', '88.1% (Monte Carlo validated)'),
        ('Primary Endpoint', 'Clinical response rate at Week 12'),
        ('Inclusion', 'NOD2 mutation (R702W, G908R, or WT)'),
        ('Exclusion', 'L1007fs mutation (binding site deleted)'),
    ]

    y_start = 0.78
    for i, (param, value) in enumerate(params):
        y = y_start - i * 0.08
        ax.text(0.15, y, param + ':', ha='left', va='center',
                fontsize=11, fontweight='bold', transform=ax.transAxes)
        ax.text(0.38, y, value, ha='left', va='center',
                fontsize=11, transform=ax.transAxes)

    # Power curve inset
    ax_inset = fig.add_axes([0.55, 0.25, 0.38, 0.35])
    effect_sizes = [10, 15, 20, 25, 30, 35]
    powers = [30.0, 45.2, 74.2, 89.0, 96.4, 98.6]
    ax_inset.plot(effect_sizes, powers, 'o-', color='#1976D2', linewidth=2, markersize=8)
    ax_inset.axhline(80, color='red', linestyle='--', label='80% power threshold')
    ax_inset.axvline(25, color='green', linestyle=':', alpha=0.7, label='Expected effect (25%)')
    ax_inset.fill_between(effect_sizes, powers, alpha=0.2, color='#1976D2')
    ax_inset.set_xlabel('Effect Size (%)', fontsize=10)
    ax_inset.set_ylabel('Statistical Power (%)', fontsize=10)
    ax_inset.set_title('Power Analysis', fontsize=11, fontweight='bold')
    ax_inset.legend(fontsize=8, loc='lower right')
    ax_inset.set_ylim(0, 105)
    ax_inset.grid(True, alpha=0.3)

    # Key result box
    rect = FancyBboxPatch(
        (0.1, 0.08), 0.35, 0.12,
        boxstyle="round,pad=0.02,rounding_size=0.02",
        facecolor='#E8F5E9', edgecolor='#2E7D32', linewidth=2,
        transform=ax.transAxes
    )
    ax.add_patch(rect)
    ax.text(0.275, 0.15, 'Expected Odds Ratio: 3.42', ha='center', va='center',
            fontsize=12, fontweight='bold', color='#1B5E20', transform=ax.transAxes)
    ax.text(0.275, 0.10, '(80mg vs Placebo)', ha='center', va='center',
            fontsize=10, color='#388E3C', transform=ax.transAxes)

    plt.savefig(OUTPUT_DIR / 'figure8_clinical_trial.png', dpi=DPI, bbox_inches='tight',
                facecolor='white', edgecolor='none')
    plt.close()
    print("Saved: figure8_clinical_trial.png")


def main():
    """Generate all ISEF figures."""
    print("=" * 60)
    print("ISEF Figure Generator - NOD2-CROHN Drug Discovery")
    print("=" * 60)
    print(f"\nOutput directory: {OUTPUT_DIR}")
    print(f"DPI: {DPI}\n")

    print("Generating figures...")
    create_figure1_pipeline()
    create_figure2_rmsd_comparison()
    create_figure3_heterogeneity()
    create_figure4_fep_results()
    create_figure5_drug_effects()
    create_figure6_precision_medicine()
    create_figure7_mechanism_summary()
    create_figure8_clinical_trial()
    create_pymol_script()

    print("\n" + "=" * 60)
    print("ALL FIGURES GENERATED SUCCESSFULLY")
    print("=" * 60)
    print(f"\nFiles saved to: {OUTPUT_DIR}")
    print("\nFigures:")
    print("  1. figure1_pipeline.png - Computational pipeline flowchart")
    print("  2. figure2_rmsd_comparison.png - RMSD WT vs R702W")
    print("  3. figure3_heterogeneity.png - Conformational heterogeneity")
    print("  4. figure4_fep_binding.png - MM-GBSA binding energies")
    print("  5. figure5_drug_effects.png - Drug stabilization effects")
    print("  6. figure6_precision_medicine.png - Patient stratification")
    print("  7. figure7_mechanism_summary.png - Mechanism schematic")
    print("  8. figure8_clinical_trial.png - Clinical trial design")
    print("  9. pymol_visualization.pml - PyMOL script for 3D figures")


if __name__ == "__main__":
    main()
