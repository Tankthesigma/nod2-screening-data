"""
ISEF 2026 - NOD2 Screening Pipeline Image
High resolution for 48" wide poster
"""

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
import numpy as np

# High DPI for print quality
fig, ax = plt.subplots(figsize=(24, 6), dpi=300)
ax.set_xlim(0, 24)
ax.set_ylim(0, 6)
ax.axis('off')

# Colors
DARK_BLUE = '#1a365d'
TEAL = '#319795'
GREEN = '#276749'
CORAL = '#c53030'
GRAY = '#4a5568'
LIGHT_GRAY = '#f7fafc'
WHITE = '#ffffff'

# ============================================
# COMPOUND LIBRARIES (Left side)
# ============================================
# e-Drug3D box
edrug_box = FancyBboxPatch((0.3, 3.5), 2.5, 1.8,
                            boxstyle="round,pad=0.05,rounding_size=0.1",
                            facecolor=TEAL, edgecolor='none')
ax.add_patch(edrug_box)
ax.text(1.55, 4.7, 'e-Drug3D', fontsize=11, fontweight='bold',
        color=WHITE, ha='center', va='center')
ax.text(1.55, 4.2, 'n = 2,152', fontsize=14, fontweight='bold',
        color=WHITE, ha='center', va='center')
ax.text(1.55, 3.75, 'FDA-approved', fontsize=8, color=WHITE,
        ha='center', va='center', style='italic')

# COCONUT box
coconut_box = FancyBboxPatch((0.3, 1.2), 2.5, 1.8,
                              boxstyle="round,pad=0.05,rounding_size=0.1",
                              facecolor=GREEN, edgecolor='none')
ax.add_patch(coconut_box)
ax.text(1.55, 2.4, 'COCONUT', fontsize=11, fontweight='bold',
        color=WHITE, ha='center', va='center')
ax.text(1.55, 1.9, 'n = 7,414', fontsize=14, fontweight='bold',
        color=WHITE, ha='center', va='center')
ax.text(1.55, 1.45, 'Natural products', fontsize=8, color=WHITE,
        ha='center', va='center', style='italic')

# "COMPOUND LIBRARIES" label
ax.text(1.55, 5.6, 'COMPOUND', fontsize=8, color=GRAY, ha='center', fontweight='bold')
ax.text(1.55, 5.3, 'LIBRARIES', fontsize=8, color=GRAY, ha='center', fontweight='bold')

# Arrows from libraries
ax.annotate('', xy=(3.5, 3.3), xytext=(2.8, 4.2),
            arrowprops=dict(arrowstyle='->', color=GRAY, lw=1.5))
ax.annotate('', xy=(3.5, 3.3), xytext=(2.8, 2.2),
            arrowprops=dict(arrowstyle='->', color=GRAY, lw=1.5))

# ============================================
# PIPELINE STEPS
# ============================================
steps = [
    {'num': '1', 'tool': 'GNINA', 'name': 'Docking', 'n': '9,566', 'detail': 'CNN scoring', 'x': 4.2},
    {'num': '2', 'tool': 'NOD2-Scout', 'name': 'ML Filter', 'n': '2,129', 'detail': 'AUC = 0.89', 'x': 7.2},
    {'num': '3', 'tool': 'RDKit', 'name': 'ADMET', 'n': '148', 'detail': 'Drug-likeness', 'x': 10.2},
    {'num': '4', 'tool': 'Manual', 'name': 'Tier-1', 'n': '8', 'detail': 'Top picks', 'x': 13.2},
    {'num': '5', 'tool': 'OpenMM', 'name': 'MD Sim', 'n': '5', 'detail': '520 ns', 'x': 16.2},
    {'num': '6', 'tool': 'OpenMMTools', 'name': 'FEP', 'n': '2', 'detail': 'ΔΔG', 'x': 19.2},
]

for i, step in enumerate(steps):
    x = step['x']

    # Number circle
    circle = plt.Circle((x + 1, 5.2), 0.35, color=GRAY, zorder=10)
    ax.add_patch(circle)
    ax.text(x + 1, 5.2, step['num'], fontsize=12, fontweight='bold',
            color=WHITE, ha='center', va='center', zorder=11)

    # Tool name above
    ax.text(x + 1, 5.7, step['tool'], fontsize=9, color=GRAY,
            ha='center', va='center', style='italic')

    # Main box
    box = FancyBboxPatch((x, 2.3), 2, 2.5,
                         boxstyle="round,pad=0.05,rounding_size=0.15",
                         facecolor=WHITE, edgecolor='#e2e8f0', linewidth=2)
    ax.add_patch(box)

    # Step name
    ax.text(x + 1, 4.3, step['name'], fontsize=12, fontweight='bold',
            color=DARK_BLUE, ha='center', va='center')

    # N value
    ax.text(x + 1, 3.4, f"n = {step['n']}", fontsize=16, fontweight='bold',
            color=DARK_BLUE, ha='center', va='center')

    # Detail
    ax.text(x + 1, 2.6, step['detail'], fontsize=9, color=GRAY,
            ha='center', va='center')

    # Arrow to next step
    if i < len(steps) - 1:
        ax.annotate('', xy=(x + 2.7, 3.55), xytext=(x + 2.1, 3.55),
                    arrowprops=dict(arrowstyle='->', color='#cbd5e0', lw=2))

# ============================================
# 2 LEADS BOX (Bottom)
# ============================================
leads_box = FancyBboxPatch((6, 0.2), 12, 1.5,
                           boxstyle="round,pad=0.05,rounding_size=0.15",
                           facecolor=CORAL, edgecolor='none')
ax.add_patch(leads_box)

ax.text(12, 1.2, '2 LEADS IDENTIFIED', fontsize=16, fontweight='bold',
        color=WHITE, ha='center', va='center')
ax.text(12, 0.6, 'Febuxostat  •  Bufadienolide   (99.98% reduction)',
        fontsize=11, color=WHITE, ha='center', va='center')

# Arrows pointing to leads box
ax.annotate('', xy=(9, 1.7), xytext=(8, 2.3),
            arrowprops=dict(arrowstyle='->', color='#cbd5e0', lw=1.5,
                           connectionstyle='arc3,rad=-0.1'))
ax.annotate('', xy=(15, 1.7), xytext=(20, 2.3),
            arrowprops=dict(arrowstyle='->', color='#cbd5e0', lw=1.5,
                           connectionstyle='arc3,rad=0.1'))

plt.tight_layout()
plt.savefig(r'C:\Users\vasud\Downloads\pipeline_FIXED.png',
            dpi=300, bbox_inches='tight', facecolor='white',
            edgecolor='none', pad_inches=0.1)
plt.savefig(r'C:\Users\vasud\nod2-screening-data\isef_figures\v5_winner\ppt_photos\pipeline_FIXED.png',
            dpi=300, bbox_inches='tight', facecolor='white',
            edgecolor='none', pad_inches=0.1)
print("Saved: pipeline_FIXED.png")
print("Resolution: 7200 x 1800 pixels (24in x 6in @ 300 DPI)")
