#!/usr/bin/env python3
"""
NOD2 Domains Key/Legend SVG
"""

import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Rectangle

plt.rcParams["font.family"] = "DejaVu Sans"

fig, ax = plt.subplots(figsize=(6, 5))
ax.set_xlim(0, 6)
ax.set_ylim(0, 5)
ax.axis('off')
fig.patch.set_facecolor('none')
fig.patch.set_alpha(0)

# Title
ax.text(3, 4.7, 'NOD2 Domains', fontsize=18, fontweight='bold',
        color='#1e3a5f', ha='center', va='center')

# Domain colors
domains = [
    {'name': 'CARD', 'desc': 'Signaling', 'color': '#6366f1', 'y': 4.0},
    {'name': 'NBD', 'desc': 'Nucleotide binding', 'color': '#3b82f6', 'y': 3.3},
    {'name': 'HD2', 'desc': 'R702W mutation site', 'color': '#f97316', 'y': 2.6},
    {'name': 'LRR', 'desc': 'Ligand recognition', 'color': '#22c55e', 'y': 1.9},
]

for d in domains:
    # Color box
    box = FancyBboxPatch((0.5, d['y'] - 0.2), 0.6, 0.4,
                          boxstyle="round,pad=0.02,rounding_size=0.05",
                          facecolor=d['color'], edgecolor='none')
    ax.add_patch(box)

    # Name
    ax.text(1.4, d['y'], d['name'], fontsize=14, fontweight='bold',
            color='#1e3a5f', ha='left', va='center')

    # Description
    ax.text(2.5, d['y'], d['desc'], fontsize=12,
            color='#64748b', ha='left', va='center')

# Distance line
ax.plot([0.8, 5.2], [1.2, 1.2], color='#dc2626', linewidth=2)
ax.plot([0.8, 0.8], [1.1, 1.3], color='#dc2626', linewidth=2)
ax.plot([5.2, 5.2], [1.1, 1.3], color='#dc2626', linewidth=2)

# Distance text
ax.text(3, 1.0, '79.4 Å', fontsize=16, fontweight='bold',
        color='#dc2626', ha='center', va='center')

# Bottom text
ax.text(3, 0.6, 'R702W to binding pocket', fontsize=11,
        color='#1e3a5f', ha='center', va='center')
ax.text(3, 0.25, '(allosteric mechanism)', fontsize=10,
        color='#64748b', ha='center', va='center', style='italic')

plt.tight_layout()

plt.savefig("C:/Users/vasud/Desktop/poster photos/nod2_domains_key.svg", format="svg",
            bbox_inches="tight", pad_inches=0.1, transparent=True)
plt.savefig("C:/Users/vasud/Desktop/poster photos/nod2_domains_key.png", format="png",
            dpi=1200, bbox_inches="tight", pad_inches=0.1, transparent=True)

print("DONE - NOD2 Domains Key")
print("Files: nod2_domains_key.svg + nod2_domains_key.png")
