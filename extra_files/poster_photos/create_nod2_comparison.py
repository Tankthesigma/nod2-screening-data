#!/usr/bin/env python3
"""
NOD2 Structure: WT vs R702W Comparison
Using PyMOL renders from ppt_photos
"""

import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Rectangle, Circle
from matplotlib.offsetbox import OffsetImage, AnnotationBbox
import matplotlib.image as mpimg
from pathlib import Path

# Font
plt.rcParams["font.family"] = "DejaVu Sans"

# Paths
IMG_DIR = Path(r"C:\Users\vasud\nod2-screening-data\isef_figures\v5_winner\ppt_photos")

# Load images
wt_img = mpimg.imread(IMG_DIR / "figE1_wt_electrostatic.png")
r702w_img = mpimg.imread(IMG_DIR / "figE2_r702w_mutation_site.png")

# Figure - wide format
fig, ax = plt.subplots(figsize=(18, 8))
ax.set_xlim(0, 18)
ax.set_ylim(0, 8)
ax.axis('off')
fig.patch.set_facecolor('white')

# ===== HEADER =====
header = FancyBboxPatch((0.2, 7.1), 17.6, 0.75,
                         boxstyle="round,pad=0.02,rounding_size=0.12",
                         facecolor='#1e3a5f', edgecolor='none')
ax.add_patch(header)
ax.text(9, 7.48, 'NOD2 STRUCTURE: WILD-TYPE vs R702W', fontsize=20, fontweight='bold',
        color='white', ha='center', va='center')

# ===== WT (Left) =====
wt_card = FancyBboxPatch((0.3, 0.3), 8.2, 6.5,
                          boxstyle="round,pad=0.02,rounding_size=0.2",
                          facecolor='white', edgecolor='#3b82f6', linewidth=3)
ax.add_patch(wt_card)

# WT header
wt_header = FancyBboxPatch((0.3, 6.1), 8.2, 0.7,
                            boxstyle="round,pad=0.02,rounding_size=0.1",
                            facecolor='#3b82f6', edgecolor='none')
ax.add_patch(wt_header)
ax.text(4.4, 6.45, 'WILD-TYPE NOD2', fontsize=16, fontweight='bold',
        color='white', ha='center', va='center')

# WT structure image
wt_imgbox = OffsetImage(wt_img, zoom=0.28)
wt_ab = AnnotationBbox(wt_imgbox, (4.4, 3.3), frameon=False)
ax.add_artist(wt_ab)

# WT info
ax.text(4.4, 0.8, 'R702 in HD2 domain', fontsize=11, color='#1e40af',
        ha='center', va='center', fontweight='bold')
ax.text(4.4, 0.45, 'Normal conformational dynamics', fontsize=10, color='#64748b',
        ha='center', va='center')

# ===== R702W (Right) =====
r702w_card = FancyBboxPatch((9.5, 0.3), 8.2, 6.5,
                             boxstyle="round,pad=0.02,rounding_size=0.2",
                             facecolor='white', edgecolor='#dc2626', linewidth=3)
ax.add_patch(r702w_card)

# R702W header
r702w_header = FancyBboxPatch((9.5, 6.1), 8.2, 0.7,
                               boxstyle="round,pad=0.02,rounding_size=0.1",
                               facecolor='#dc2626', edgecolor='none')
ax.add_patch(r702w_header)
ax.text(13.6, 6.45, 'R702W MUTANT', fontsize=16, fontweight='bold',
        color='white', ha='center', va='center')

# R702W structure image
r702w_imgbox = OffsetImage(r702w_img, zoom=0.28)
r702w_ab = AnnotationBbox(r702w_imgbox, (13.6, 3.3), frameon=False)
ax.add_artist(r702w_ab)

# R702W info
ax.text(13.6, 0.8, 'W702 - Crohn\'s variant (3x risk)', fontsize=11, color='#991b1b',
        ha='center', va='center', fontweight='bold')
ax.text(13.6, 0.45, 'Reduced conformational sampling', fontsize=10, color='#64748b',
        ha='center', va='center')

# ===== VS Arrow =====
# Arrow from WT to R702W
ax.annotate('', xy=(9.3, 3.5), xytext=(8.7, 3.5),
            arrowprops=dict(arrowstyle='->', color='#1e3a5f', lw=3))
ax.text(9, 4.0, 'Mutation', fontsize=10, fontweight='bold', color='#1e3a5f',
        ha='center', va='center')

plt.tight_layout()

# Save
plt.savefig("C:/Users/vasud/Desktop/poster photos/nod2_comparison.svg", format="svg",
            bbox_inches="tight", pad_inches=0.05, facecolor="white")
plt.savefig("C:/Users/vasud/Desktop/poster photos/nod2_comparison.png", format="png",
            dpi=1200, bbox_inches="tight", pad_inches=0.05, facecolor="white")

print("DONE - NOD2 WT vs R702W Comparison")
print("Files: nod2_comparison.svg + nod2_comparison.png (1200 DPI)")
