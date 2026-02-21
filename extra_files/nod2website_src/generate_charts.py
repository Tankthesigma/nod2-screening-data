"""
Generate scientific charts for NOD2 Drug Discovery website
"""
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import os

# Create output directory
os.makedirs('public/charts', exist_ok=True)

# Set style for publication-quality figures
plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
    'font.size': 12,
    'axes.labelsize': 14,
    'axes.titlesize': 16,
    'xtick.labelsize': 11,
    'ytick.labelsize': 11,
    'legend.fontsize': 11,
    'figure.dpi': 150,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
    'axes.spines.top': False,
    'axes.spines.right': False,
})

# ============================================================
# FIGURE 1: Binding Free Energy Bar Chart
# ============================================================
fig1, ax1 = plt.subplots(figsize=(10, 6))

compounds = ['Febuxostat\n(WT)', 'Febuxostat\n(R702W)', 'Bufadienolide\n(WT)', 'Bufadienolide\n(R702W)']
values = [-10.36, -8.02, -15.22, -15.66]
errors = [0.18, 0.19, 0.26, 0.26]
colors = ['#3b82f6', '#93c5fd', '#14b8a6', '#5eead4']

bars = ax1.bar(compounds, values, yerr=errors, capsize=8, color=colors,
               edgecolor='black', linewidth=1.5, error_kw={'linewidth': 2})

ax1.axhline(y=0, color='black', linewidth=0.5)
ax1.set_ylabel('Binding Free Energy (kcal/mol)', fontweight='bold')
ax1.set_title('Figure 1: FEP-Calculated Binding Affinities', fontweight='bold', pad=15)

# Add value labels on bars
for bar, val, err in zip(bars, values, errors):
    height = bar.get_height()
    ax1.annotate(f'{val:.2f} ± {err:.2f}',
                xy=(bar.get_x() + bar.get_width() / 2, height - 0.8),
                ha='center', va='top', fontsize=11, fontweight='bold', color='white')

ax1.set_ylim(-18, 0)
ax1.set_yticks([0, -5, -10, -15])

# Add significance annotation
ax1.annotate('', xy=(0.15, -8.02), xytext=(0.15, -10.36),
            arrowprops=dict(arrowstyle='<->', color='red', lw=2))
ax1.text(0.4, -9.2, 'ΔΔG = +2.34\n(50× weaker)', fontsize=10, color='red', fontweight='bold')

plt.tight_layout()
plt.savefig('public/charts/binding_energy.png', transparent=False, facecolor='white')
plt.close()
print("Created: binding_energy.png")

# ============================================================
# FIGURE 2: Delta-Delta-G Comparison
# ============================================================
fig2, ax2 = plt.subplots(figsize=(8, 5))

compounds_ddg = ['Febuxostat', 'Bufadienolide']
ddg_values = [2.34, -0.44]
ddg_errors = [0.26, 0.37]
colors_ddg = ['#ef4444', '#14b8a6']

bars2 = ax2.bar(compounds_ddg, ddg_values, yerr=ddg_errors, capsize=10,
                color=colors_ddg, edgecolor='black', linewidth=2,
                error_kw={'linewidth': 2.5}, width=0.5)

ax2.axhline(y=0, color='black', linewidth=2, linestyle='--', label='No effect from mutation')
ax2.set_ylabel('ΔΔG (kcal/mol)', fontweight='bold')
ax2.set_title('Figure 2: Mutation Effect on Binding', fontweight='bold', pad=15)

# Add value labels
for bar, val, err in zip(bars2, ddg_values, ddg_errors):
    height = bar.get_height()
    y_pos = height + 0.3 if height > 0 else height - 0.5
    ax2.annotate(f'{val:+.2f} ± {err:.2f}',
                xy=(bar.get_x() + bar.get_width() / 2, y_pos),
                ha='center', va='bottom' if height > 0 else 'top',
                fontsize=14, fontweight='bold')

# Add interpretation labels
ax2.text(0, 1.7, 'MUTATION\nSENSITIVE', ha='center', fontsize=10, color='#991b1b', fontweight='bold')
ax2.text(1, -1.3, 'MUTATION\nRESISTANT', ha='center', fontsize=10, color='#047857', fontweight='bold')

ax2.set_ylim(-2, 3.5)
ax2.legend(loc='upper right')

plt.tight_layout()
plt.savefig('public/charts/ddg_comparison.png', transparent=False, facecolor='white')
plt.close()
print("Created: ddg_comparison.png")

# ============================================================
# FIGURE 3: Screening Funnel
# ============================================================
fig3, ax3 = plt.subplots(figsize=(10, 7))

stages = ['Initial Library', 'GNINA Docking', 'ML Ranking', 'ADMET Filter', 'MD Simulation', 'FEP Validated']
counts = [9566, 500, 100, 8, 4, 2]
colors_funnel = ['#3b82f6', '#6366f1', '#8b5cf6', '#a855f7', '#d946ef', '#14b8a6']

# Create horizontal bar chart (funnel-like)
y_pos = np.arange(len(stages))
bars3 = ax3.barh(y_pos, counts, color=colors_funnel, edgecolor='black', linewidth=1.5, height=0.7)

ax3.set_yticks(y_pos)
ax3.set_yticklabels(stages)
ax3.invert_yaxis()  # Top to bottom
ax3.set_xlabel('Number of Compounds', fontweight='bold')
ax3.set_title('Figure 3: Computational Screening Funnel', fontweight='bold', pad=15)
ax3.set_xscale('log')

# Add count labels
for bar, count in zip(bars3, counts):
    width = bar.get_width()
    ax3.annotate(f'{count:,}',
                xy=(width, bar.get_y() + bar.get_height()/2),
                xytext=(5, 0), textcoords='offset points',
                ha='left', va='center', fontsize=12, fontweight='bold')

# Add enrichment annotation
ax3.text(5000, 5.5, f'Hit Rate: 2/9,566 = 0.021%\nEnrichment: 4,783× over random',
         fontsize=11, style='italic', bbox=dict(boxstyle='round', facecolor='#f0fdf4', edgecolor='#14b8a6'))

plt.tight_layout()
plt.savefig('public/charts/screening_funnel.png', transparent=False, facecolor='white')
plt.close()
print("Created: screening_funnel.png")

# ============================================================
# FIGURE 4: Binding Energy with WT vs MUT grouped
# ============================================================
fig4, ax4 = plt.subplots(figsize=(9, 6))

x = np.arange(2)
width = 0.35

wt_values = [-10.36, -15.22]
mut_values = [-8.02, -15.66]
wt_errors = [0.18, 0.26]
mut_errors = [0.19, 0.26]

bars_wt = ax4.bar(x - width/2, wt_values, width, yerr=wt_errors, label='Wild-Type (WT)',
                  color=['#3b82f6', '#14b8a6'], edgecolor='black', linewidth=1.5, capsize=6)
bars_mut = ax4.bar(x + width/2, mut_values, width, yerr=mut_errors, label='R702W Mutant',
                   color=['#93c5fd', '#5eead4'], edgecolor='black', linewidth=1.5, capsize=6,
                   hatch='///')

ax4.set_ylabel('ΔG (kcal/mol)', fontweight='bold')
ax4.set_title('Figure 4: Wild-Type vs Mutant Binding', fontweight='bold', pad=15)
ax4.set_xticks(x)
ax4.set_xticklabels(['Febuxostat', 'Bufadienolide'], fontweight='bold')
ax4.legend()
ax4.axhline(y=0, color='black', linewidth=0.5)
ax4.set_ylim(-18, 0)

# Add ΔΔG annotations
ax4.annotate('ΔΔG = +2.34\n(weaker)', xy=(0, -9), fontsize=10, ha='center',
            color='#dc2626', fontweight='bold',
            bbox=dict(boxstyle='round', facecolor='#fef2f2', edgecolor='#dc2626'))
ax4.annotate('ΔΔG = -0.44\n(unchanged)', xy=(1, -13), fontsize=10, ha='center',
            color='#047857', fontweight='bold',
            bbox=dict(boxstyle='round', facecolor='#f0fdf4', edgecolor='#047857'))

plt.tight_layout()
plt.savefig('public/charts/wt_vs_mut.png', transparent=False, facecolor='white')
plt.close()
print("Created: wt_vs_mut.png")

print("\n✓ All charts generated in public/charts/")
