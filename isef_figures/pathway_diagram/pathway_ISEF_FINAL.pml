# ============================================================================
# ISEF 2026 - NATURE PAPER STYLE
# Resolution: 4000x3000 (4K)
# Based on actual Nature Communications NOD2 figures
# ============================================================================

reinitialize

set ray_trace_mode, 1
set antialias, 2
set cartoon_sampling, 14
set ambient, 0.45
set ray_opaque_background, 0
set orthoscopic, 1
set depth_cue, 0
set fog, 0
set spec_reflect, 0.25
set spec_power, 150
set cartoon_fancy_helices, 1
set cartoon_smooth_loops, 1
set stick_radius, 0.15
set sphere_scale, 0.3

bg_color white
cd C:/Users/vasud/nod2-screening-data/isef_figures/pathway_diagram

# ============================================================================
# DOMAIN COLORS (professional palette)
# ============================================================================
set_color card_gray, [0.6, 0.6, 0.6]
set_color nbd_blue, [0.2, 0.5, 0.85]
set_color hd1_teal, [0.1, 0.6, 0.6]
set_color hd2_orange, [1.0, 0.55, 0.1]
set_color lrr_green, [0.3, 0.7, 0.3]
set_color mutation_red, [0.95, 0.2, 0.2]
set_color pocket_magenta, [0.85, 0.2, 0.85]

# ============================================================================
# IMAGE 1: WT NOD2 FULL STRUCTURE
# ============================================================================
print "=== IMAGE 1: WT NOD2 Full ==="

load C:/Users/vasud/nod2-screening-data/PHASE_6/structures/NOD2_alphafold_full.pdb, wt

hide all
show cartoon, wt

# Domain colors
color card_gray, wt and resi 1-191
color nbd_blue, wt and resi 273-577
color hd1_teal, wt and resi 578-628
color hd2_orange, wt and resi 629-743
color lrr_green, wt and resi 744-1040

# R702 = GIGANTIC CYAN - IN YOUR FACE!!!
select mut_wt, wt and resi 702
show spheres, mut_wt
set sphere_scale, 1.8, mut_wt
color cyan, mut_wt

# Binding pocket = MAGENTA SPHERES
select pocket_wt, wt and resi 1008+1011+1037
show spheres, pocket_wt
set sphere_scale, 1.2, pocket_wt
color pocket_magenta, pocket_wt

# NO DISTANCE LINE - let spatial gap speak

# Camera
orient wt
turn y, -50
turn x, 10
zoom wt, 0

ray 4000, 3000
png img1_wt_full_ISEF.png
print "Saved: img1_wt_full_ISEF.png"

# ============================================================================
# IMAGE 2: R702W NOD2 FULL - SAME ANGLE
# ============================================================================
print "=== IMAGE 2: R702W NOD2 Full ==="

delete all
load C:/Users/vasud/nod2-screening-data/PHASE_A1_mutant_MD/structures/NOD2_R702W.pdb, r702w

hide all
show cartoon, r702w

# Same domain colors
color card_gray, r702w and resi 1-191
color nbd_blue, r702w and resi 273-577
color hd1_teal, r702w and resi 578-628
color hd2_orange, r702w and resi 629-743
color lrr_green, r702w and resi 744-1040

# W702 = GIGANTIC CYAN - IN YOUR FACE!!!
select mut_r702w, r702w and resi 702
show spheres, mut_r702w
set sphere_scale, 4.0, mut_r702w
color cyan, mut_r702w

# Binding pocket = MAGENTA SPHERES
select pocket_r702w, r702w and resi 1008+1011+1037
show spheres, pocket_r702w
set sphere_scale, 1.2, pocket_r702w
color pocket_magenta, pocket_r702w

# SAME camera as WT
orient r702w
turn y, -50
turn x, 10
zoom r702w, 0

ray 4000, 3000
png img2_r702w_full_ISEF.png
print "Saved: img2_r702w_full_ISEF.png"

# ============================================================================
# POCKET SETUP - shared for all 4 pocket images
# ============================================================================
print "=== Setting up pocket ==="

delete all
load C:/Users/vasud/nod2-screening-data/PHASE_6/structures/NOD2_LRR_clean.pdb, lrr

hide all
show cartoon, lrr
color lrr_green, lrr

# Key residues as sticks
select key_res, lrr and resi 1008+1011+1037
show sticks, key_res
color pocket_magenta, key_res
set stick_radius, 0.18, key_res

# ============================================================================
# IMAGE 3: WT + BUFADIENOLIDE (CYAN sticks)
# ============================================================================
print "=== IMAGE 3: WT + Bufadienolide ==="

load C:/Users/vasud/nod2-screening-data/PHASE_6/structures/natural_top_docked.sdf, bufa

# Ligand as STICKS - CYAN carbons (Nature style)
show sticks, bufa
set stick_radius, 0.22, bufa
util.cbac bufa

# ZOOM on ligand - fills frame
center bufa
orient bufa
turn y, 25
turn x, -15
zoom bufa, 4

# Store this EXACT view
scene pocket_view, store

ray 4000, 3000
png img3_wt_bufa_ISEF.png
print "Saved: img3_wt_bufa_ISEF.png"

# ============================================================================
# IMAGE 4: R702W + BUFADIENOLIDE (same view)
# ============================================================================
print "=== IMAGE 4: R702W + Bufadienolide ==="

scene pocket_view, recall
ray 4000, 3000
png img4_r702w_bufa_ISEF.png
print "Saved: img4_r702w_bufa_ISEF.png"

# ============================================================================
# IMAGE 5: WT + FEBUXOSTAT (ORANGE sticks)
# ============================================================================
print "=== IMAGE 5: WT + Febuxostat ==="

delete bufa
load C:/Users/vasud/nod2-screening-data/PHASE_6/structures/febuxostat_docked.sdf, febu

# Ligand as STICKS - ORANGE carbons
show sticks, febu
set stick_radius, 0.22, febu
util.cbao febu

# ZOOM on ligand
center febu
orient febu
turn y, 25
turn x, -15
zoom febu, 4

# Store this view
scene pocket_febu, store

ray 4000, 3000
png img5_wt_febu_ISEF.png
print "Saved: img5_wt_febu_ISEF.png"

# ============================================================================
# IMAGE 6: R702W + FEBUXOSTAT (same view)
# ============================================================================
print "=== IMAGE 6: R702W + Febuxostat ==="

scene pocket_febu, recall
ray 4000, 3000
png img6_r702w_febu_ISEF.png
print "Saved: img6_r702w_febu_ISEF.png"

# ============================================================================
# DONE
# ============================================================================
print ""
print "=============================================="
print "NATURE PAPER STYLE IMAGES COMPLETE"
print "=============================================="
print ""
print "FULL STRUCTURES:"
print "  - Cartoon by domain"
print "  - R702/W702 = RED SPHERES (scale 1.8)"
print "  - Pocket = MAGENTA SPHERES"
print "  - NO distance line"
print ""
print "POCKET VIEWS:"
print "  - Ligand = STICKS (util.cbac/cbao)"
print "  - Bufadienolide = CYAN carbons"
print "  - Febuxostat = ORANGE carbons"
print "  - ZOOMED on ligand"
print "  - Key residues as magenta sticks"
print "=============================================="
