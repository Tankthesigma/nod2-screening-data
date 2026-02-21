# ============================================================================
# PUBLICATION QUALITY - Nature/Cell Style
# Surface pocket + bright ligand, NO clutter
# Resolution: 4000x3000 (4K)
# ============================================================================

reinitialize

# High quality settings
set ray_trace_mode, 1
set antialias, 2
set ray_shadows, 1
set ambient, 0.4
set direct, 0.5
set specular, 0.3
set ray_opaque_background, 0
set orthoscopic, 1
set depth_cue, 0
set fog, 0
set surface_quality, 1
set solvent_radius, 1.6
set transparency_mode, 2

bg_color white
cd C:/Users/vasud/nod2-screening-data/isef_figures/pathway_diagram

# ============================================================================
# LOAD ALL STRUCTURES
# ============================================================================
print "Loading structures..."

# Load LRR domain for pocket
load C:/Users/vasud/nod2-screening-data/PHASE_6/structures/NOD2_LRR_clean.pdb, pocket

# Load ligands
load C:/Users/vasud/nod2-screening-data/PHASE_6/structures/natural_top_docked.sdf, bufa
load C:/Users/vasud/nod2-screening-data/PHASE_6/structures/febuxostat_docked.sdf, febu

# ============================================================================
# POCKET SURFACE SETUP (shared for all images)
# ============================================================================
hide all

# Show pocket as semi-transparent surface - GRAY/WHITE
show surface, pocket
set surface_color, gray80, pocket
set transparency, 0.75, pocket

# ============================================================================
# IMAGE 1: WT + BUFADIENOLIDE
# ============================================================================
print "=== IMAGE 1: WT + Bufadienolide ==="

# Hide febuxostat
disable febu

# Show bufadienolide as BRIGHT CYAN SPHERES
enable bufa
show spheres, bufa
set sphere_scale, 0.35, bufa
color cyan, bufa

# Camera: Center on ligand, pocket wraps around it
center bufa
zoom bufa, 3
orient bufa
turn y, 15
turn x, -10

# Store this EXACT view for all pocket images
scene pocket_master, store

ray 4000, 3000
png pocket1_wt_bufa.png
print "Saved: pocket1_wt_bufa.png"

# ============================================================================
# IMAGE 2: R702W + BUFADIENOLIDE (same view - LRR identical)
# ============================================================================
print "=== IMAGE 2: R702W + Bufadienolide ==="

# Same exact view - mutation is in HD2, not LRR
scene pocket_master, recall
ray 4000, 3000
png pocket2_r702w_bufa.png
print "Saved: pocket2_r702w_bufa.png"

# ============================================================================
# IMAGE 3: WT + FEBUXOSTAT
# ============================================================================
print "=== IMAGE 3: WT + Febuxostat ==="

# Hide bufadienolide, show febuxostat
disable bufa
enable febu

# Show febuxostat as BRIGHT ORANGE SPHERES
show spheres, febu
set sphere_scale, 0.35, febu
color orange, febu

# Camera: Same pocket view but centered on febuxostat
center febu
zoom febu, 3
orient febu
turn y, 15
turn x, -10

# Store febuxostat view
scene pocket_febu, store

ray 4000, 3000
png pocket3_wt_febu.png
print "Saved: pocket3_wt_febu.png"

# ============================================================================
# IMAGE 4: R702W + FEBUXOSTAT (same view)
# ============================================================================
print "=== IMAGE 4: R702W + Febuxostat ==="

scene pocket_febu, recall
ray 4000, 3000
png pocket4_r702w_febu.png
print "Saved: pocket4_r702w_febu.png"

# ============================================================================
# BONUS: STICKS VERSION (Nature style - element coloring)
# ============================================================================
print "=== BONUS: Sticks version ==="

# Reset
hide spheres, bufa
hide spheres, febu

# Bufadienolide as CYAN sticks
enable bufa
disable febu
show sticks, bufa
set stick_radius, 0.2, bufa
util.cbac bufa
scene pocket_master, recall
ray 4000, 3000
png pocket1_wt_bufa_sticks.png
print "Saved: pocket1_wt_bufa_sticks.png"

# R702W version
ray 4000, 3000
png pocket2_r702w_bufa_sticks.png
print "Saved: pocket2_r702w_bufa_sticks.png"

# Febuxostat as ORANGE sticks
disable bufa
enable febu
show sticks, febu
set stick_radius, 0.2, febu
util.cbao febu
scene pocket_febu, recall
ray 4000, 3000
png pocket3_wt_febu_sticks.png
print "Saved: pocket3_wt_febu_sticks.png"

# R702W version
ray 4000, 3000
png pocket4_r702w_febu_sticks.png
print "Saved: pocket4_r702w_febu_sticks.png"

# ============================================================================
# DONE
# ============================================================================
print ""
print "=============================================="
print "PUBLICATION QUALITY IMAGES COMPLETE"
print "=============================================="
print ""
print "SPHERES VERSION (bright, obvious):"
print "  pocket1_wt_bufa.png"
print "  pocket2_r702w_bufa.png"
print "  pocket3_wt_febu.png"
print "  pocket4_r702w_febu.png"
print ""
print "STICKS VERSION (Nature style):"
print "  pocket1_wt_bufa_sticks.png"
print "  pocket2_r702w_bufa_sticks.png"
print "  pocket3_wt_febu_sticks.png"
print "  pocket4_r702w_febu_sticks.png"
print ""
print "Features:"
print "  - Semi-transparent gray surface (pocket)"
print "  - Bright ligand (cyan=Bufa, orange=Febu)"
print "  - Transparent background"
print "  - NO clutter, NO residue sticks"
print "  - Ligand centered, fills frame"
print "=============================================="
