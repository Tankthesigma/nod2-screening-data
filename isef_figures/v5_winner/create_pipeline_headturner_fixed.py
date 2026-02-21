"""
ISEF 2026 - NOD2 Screening Pipeline Image
Matching the headturner style but with CORRECT tools
High resolution for 48" wide poster
"""

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch, Rectangle, Circle
import numpy as np

# High DPI for print quality - matching headturner dimensions
fig, ax = plt.subplots(figsize=(32, 12), dpi=300)
ax.set_xlim(0, 32)
ax.set_ylim(0, 12)
ax.axis('off')
fig.patch.set_facecolor('white')

# Colors matching headturner
DARK_NAVY = '#1a202c'
TEAL = '#319795'
GREEN = '#276749'
CORAL = '#e53e3e'
GRAY_DARK = '#4a5568'
GRAY_LIGHT = '#a0aec0'
WHITE = '#ffffff'
BG_LIGHT = '#f7fafc'

# ============================================
# HEADER BAR
# ============================================
header = Rectangle((0, 10.5), 32, 1.5, facecolor=DARK_NAVY, edgecolor='none')
ax.add_patch(header)
ax.text(16, 11.25, 'VIRTUAL SCREENING PIPELINE', fontsize=28, fontweight='bold',
        color=WHITE, ha='center', va='center')

# ============================================
# COMPOUND LIBRARIES (Left side)
# ============================================
ax.text(2.5, 9.8, 'COMPOUND', fontsize=12, color=GRAY_DARK, ha='center', fontweight='bold')
ax.text(2.5, 9.3, 'LIBRARIES', fontsize=12, color=GRAY_DARK, ha='center', fontweight='bold')

# e-Drug3D box
edrug_box = FancyBboxPatch((0.5, 6.5), 4, 2.5,
                            boxstyle="round,pad=0.02,rounding_size=0.2",
                            facecolor=TEAL, edgecolor='none')
ax.add_patch(edrug_box)
ax.text(2.5, 8.3, 'e-Drug3D', fontsize=14, fontweight='bold',
        color=WHITE, ha='center', va='center')
ax.text(2.5, 7.5, '2,152', fontsize=28, fontweight='bold',
        color=WHITE, ha='center', va='center')
ax.text(2.5, 6.85, 'FDA drugs', fontsize=11, color='#b2f5ea',
        ha='center', va='center')

# COCONUT box
coconut_box = FancyBboxPatch((0.5, 3.3), 4, 2.5,
                              boxstyle="round,pad=0.02,rounding_size=0.2",
                              facecolor=GREEN, edgecolor='none')
ax.add_patch(coconut_box)
ax.text(2.5, 5.1, 'COCONUT', fontsize=14, fontweight='bold',
        color=WHITE, ha='center', va='center')
ax.text(2.5, 4.3, '7,414', fontsize=28, fontweight='bold',
        color=WHITE, ha='center', va='center')
ax.text(2.5, 3.65, 'Natural products', fontsize=11, color='#9ae6b4',
        ha='center', va='center')

# Arrows from libraries to first step
ax.annotate('', xy=(5.8, 6.2), xytext=(4.6, 7.5),
            arrowprops=dict(arrowstyle='-', color=GRAY_DARK, lw=2))
ax.annotate('', xy=(5.8, 6.2), xytext=(4.6, 4.8),
            arrowprops=dict(arrowstyle='-', color=GRAY_DARK, lw=2))
ax.annotate('', xy=(6.5, 6.2), xytext=(5.8, 6.2),
            arrowprops=dict(arrowstyle='->', color=GRAY_DARK, lw=2))

# ============================================
# PIPELINE STEPS - CORRECTED TOOLS
# ============================================
steps = [
    {'num': '1', 'tool': 'GNINA', 'name': 'Docking', 'n': '9,566', 'detail': 'CNN scoring', 'x': 6.2},
    {'num': '2', 'tool': 'NOD2-Scout', 'name': 'ML Filter', 'n': '2,129', 'detail': 'AUC = 0.89', 'x': 10.2},
    {'num': '3', 'tool': 'RDKit', 'name': 'ADMET', 'n': '148', 'detail': 'Drug-likeness', 'x': 14.2},
    {'num': '4', 'tool': 'Selection', 'name': 'Tier-1', 'n': '8', 'detail': 'Top picks', 'x': 18.2},
    {'num': '5', 'tool': 'OpenMM', 'name': 'MD Sim', 'n': '5', 'detail': '520 ns', 'x': 22.2},
    {'num': '6', 'tool': 'OpenMMTools', 'name': 'FEP', 'n': '2', 'detail': 'ΔΔG', 'x': 26.2},
]

for i, step in enumerate(steps):
    x = step['x']

    # Number circle (dark)
    circle = Circle((x + 1.5, 9.3), 0.5, facecolor=DARK_NAVY, edgecolor='none', zorder=10)
    ax.add_patch(circle)
    ax.text(x + 1.5, 9.3, step['num'], fontsize=16, fontweight='bold',
            color=WHITE, ha='center', va='center', zorder=11)

    # Tool name above circle
    ax.text(x + 1.5, 10.1, step['tool'], fontsize=11, color=GRAY_DARK,
            ha='center', va='center', style='italic')

    # Main white box
    box = FancyBboxPatch((x, 4.5), 3, 4.2,
                         boxstyle="round,pad=0.02,rounding_size=0.2",
                         facecolor=WHITE, edgecolor='#e2e8f0', linewidth=2)
    ax.add_patch(box)

    # Step name
    ax.text(x + 1.5, 7.8, step['name'], fontsize=16, fontweight='bold',
            color=DARK_NAVY, ha='center', va='center')

    # N value (big)
    ax.text(x + 1.5, 6.3, step['n'], fontsize=32, fontweight='bold',
            color=DARK_NAVY, ha='center', va='center')

    # Detail text
    ax.text(x + 1.5, 5.0, step['detail'], fontsize=11, color=GRAY_LIGHT,
            ha='center', va='center')

    # Arrow to next step
    if i < len(steps) - 1:
        ax.annotate('', xy=(x + 3.8, 6.6), xytext=(x + 3.1, 6.6),
                    arrowprops=dict(arrowstyle='->', color='#cbd5e0', lw=2.5))

# ============================================
# 99.98% REDUCTION BOX (Bottom)
# ============================================
reduction_box = FancyBboxPatch((8, 0.8), 16, 2.5,
                               boxstyle="round,pad=0.02,rounding_size=0.2",
                               facecolor=CORAL, edgecolor='none')
ax.add_patch(reduction_box)

ax.text(16, 2.4, '99.98% REDUCTION', fontsize=24, fontweight='bold',
        color=WHITE, ha='center', va='center')
ax.text(16, 1.4, '9,566 → 2 LEADS  •  Febuxostat & Bufadienolide',
        fontsize=14, color=WHITE, ha='center', va='center')

# Arrows pointing down to reduction box from ML and FEP
ax.plot([11.7, 11.7], [4.5, 3.5], color='#cbd5e0', lw=2, linestyle='--')
ax.plot([11.7, 12], [3.5, 3.3], color='#cbd5e0', lw=2, linestyle='--')

ax.plot([27.7, 27.7], [4.5, 3.5], color='#cbd5e0', lw=2, linestyle='--')
ax.plot([27.7, 24], [3.5, 3.3], color='#cbd5e0', lw=2, linestyle='--')

plt.tight_layout()

# Save to Downloads and project folder
plt.savefig(r'C:\Users\vasud\Downloads\pipeline_FIXED_v2.png',
            dpi=300, bbox_inches='tight', facecolor='white',
            edgecolor='none', pad_inches=0.2)
plt.savefig(r'C:\Users\vasud\nod2-screening-data\isef_figures\v5_winner\ppt_photos\pipeline_FIXED_v2.png',
            dpi=300, bbox_inches='tight', facecolor='white',
            edgecolor='none', pad_inches=0.2)

print("=" * 50)
print("SAVED: pipeline_FIXED_v2.png")
print("Location: Downloads folder")
print("Resolution: ~9600 x 3600 pixels (32in x 12in @ 300 DPI)")
print("=" * 50)
print("\nCORRECTED TOOLS:")
print("  Step 3: RDKit (was 'Filter')")
print("  Step 5: OpenMM (was 'GROMACS')")
print("  Step 6: OpenMMTools (was 'pmx')")
