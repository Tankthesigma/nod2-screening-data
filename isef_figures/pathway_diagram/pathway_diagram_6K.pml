# PyMOL Pathway Diagram Script - 6K Resolution (IMPROVED VIEWS)
# ISEF 2026 - NOD2 Project
# Author: Tanmay Vasudeva
# Resolution: 6000x4500 (6K)

# ============================================================================
# INITIALIZATION AND SETTINGS
# ============================================================================

reinitialize

# High-quality rendering settings
set ray_trace_mode, 1
set antialias, 2
set cartoon_sampling, 20
set ambient, 0.4
set ray_opaque_background, 1
set orthoscopic, 1
set depth_cue, 0
set fog, 0
set spec_reflect, 0.3
set spec_power, 200
set cartoon_fancy_helices, 1
set cartoon_highlight_color, gray75
set stick_radius, 0.25
set sphere_scale, 0.4
set surface_quality, 1
set transparency, 0.5

bg_color white
cd C:/Users/vasud/nod2-screening-data/isef_figures/pathway_diagram

# ============================================================================
# COLOR DEFINITIONS
# ============================================================================
set_color domain_gray, [0.7, 0.7, 0.7]
set_color domain_marine, [0.0, 0.4, 0.8]
set_color domain_cyan, [0.0, 0.8, 0.8]
set_color domain_orange, [1.0, 0.55, 0.0]
set_color domain_lime, [0.2, 0.8, 0.2]
set_color mutation_red, [1.0, 0.0, 0.0]
set_color pocket_magenta, [0.8, 0.2, 0.8]
set_color ligand_cyan, [0.0, 0.9, 0.9]
set_color ligand_orange, [1.0, 0.5, 0.0]
set_color pocket_surface, [0.6, 0.9, 0.6]

# ============================================================================
# IMAGE 1: WT NOD2 FULL STRUCTURE - CLEAN SIDE VIEW
# ============================================================================
print "=== Generating Image 1: WT NOD2 Full Structure (Clean Side View) ==="

load C:/Users/vasud/nod2-screening-data/PHASE_6/structures/NOD2_alphafold_full.pdb, wt_full

hide all
show cartoon, wt_full

# Color by domain
select card1_wt, wt_full and resi 1-93
select card2_wt, wt_full and resi 104-191
select nbd_wt, wt_full and resi 273-577
select hd1_wt, wt_full and resi 578-628
select hd2_wt, wt_full and resi 629-743
select lrr_wt, wt_full and resi 744-1040

color domain_gray, card1_wt
color domain_gray, card2_wt
color domain_marine, nbd_wt
color domain_cyan, hd1_wt
color domain_orange, hd2_wt
color domain_lime, lrr_wt

# Highlight R702 position with larger spheres
select r702_wt, wt_full and resi 702
show spheres, r702_wt
set sphere_scale, 0.6, r702_wt
color mutation_red, r702_wt

# Mark binding pocket location
select pocket_wt, wt_full and resi 1007+1008+1011+1034+1037
show spheres, pocket_wt
set sphere_scale, 0.5, pocket_wt
color pocket_magenta, pocket_wt

# Add distance line showing allosteric distance (~80 Angstroms)
distance allosteric_dist, r702_wt and name CA, pocket_wt and resi 1008 and name CA
hide labels, allosteric_dist
set dash_color, yellow, allosteric_dist
set dash_width, 3, allosteric_dist
set dash_gap, 0.3, allosteric_dist

# CLEAN SIDE VIEW - domains stacked vertically, CARD top, LRR bottom
orient wt_full
turn x, 90
turn y, -30
zoom wt_full, 0

ray 6000, 4500
png img1_wt_nod2_full.png

print "Saved: img1_wt_nod2_full.png"

# ============================================================================
# IMAGE 2: R702W NOD2 FULL STRUCTURE (SAME ANGLE)
# ============================================================================
print "=== Generating Image 2: R702W NOD2 Full Structure ==="

delete all
load C:/Users/vasud/nod2-screening-data/PHASE_A1_mutant_MD/structures/NOD2_R702W.pdb, r702w_full

hide all
show cartoon, r702w_full

# Color by domain (same as WT)
select card1_mut, r702w_full and resi 1-93
select card2_mut, r702w_full and resi 104-191
select nbd_mut, r702w_full and resi 273-577
select hd1_mut, r702w_full and resi 578-628
select hd2_mut, r702w_full and resi 629-743
select lrr_mut, r702w_full and resi 744-1040

color domain_gray, card1_mut
color domain_gray, card2_mut
color domain_marine, nbd_mut
color domain_cyan, hd1_mut
color domain_orange, hd2_mut
color domain_lime, lrr_mut

# Highlight W702 position (mutant Tryptophan) - larger spheres
select w702_mut, r702w_full and resi 702
show spheres, w702_mut
set sphere_scale, 0.6, w702_mut
color mutation_red, w702_mut

# Mark binding pocket
select pocket_mut, r702w_full and resi 1007+1008+1011+1034+1037
show spheres, pocket_mut
set sphere_scale, 0.5, pocket_mut
color pocket_magenta, pocket_mut

# Add distance line
distance allosteric_dist2, w702_mut and name CA, pocket_mut and resi 1008 and name CA
hide labels, allosteric_dist2
set dash_color, yellow, allosteric_dist2
set dash_width, 3, allosteric_dist2
set dash_gap, 0.3, allosteric_dist2

# SAME clean side view as WT
orient r702w_full
turn x, 90
turn y, -30
zoom r702w_full, 0

ray 6000, 4500
png img2_r702w_nod2_full.png

print "Saved: img2_r702w_nod2_full.png"

# ============================================================================
# IMAGE 3: WT LRR + BUFADIENOLIDE - TOP-DOWN POCKET VIEW
# ============================================================================
print "=== Generating Image 3: WT LRR + Bufadienolide (Top-Down Pocket) ==="

delete all

load C:/Users/vasud/nod2-screening-data/PHASE_6/structures/NOD2_LRR_clean.pdb, wt_lrr
load C:/Users/vasud/nod2-screening-data/PHASE_6/structures/natural_top_docked.sdf, bufadienolide

hide all

# Show LRR as surface (the "bowl")
show surface, wt_lrr
color pocket_surface, wt_lrr
set transparency, 0.7, wt_lrr

# Show ligand prominently
show sticks, bufadienolide
show spheres, bufadienolide
set sphere_scale, 0.25, bufadienolide
set stick_radius, 0.2, bufadienolide
color ligand_cyan, bufadienolide

# Key pocket residues as sticks (visible through surface)
select pocket_res, wt_lrr and resi 1007+1008+1009+1010+1011+1034+1035+1036+1037
show sticks, pocket_res
color pocket_magenta, pocket_res
set stick_radius, 0.15, pocket_res

# TOP-DOWN VIEW - looking into the pocket like a bowl
center bufadienolide
orient bufadienolide
turn x, -90
zoom bufadienolide, 8

ray 6000, 4500
png img3_wt_lrr_bufadienolide.png

print "Saved: img3_wt_lrr_bufadienolide.png"

# Store view for image 4
view_bufa = cmd.get_view()

# ============================================================================
# IMAGE 4: R702W LRR + BUFADIENOLIDE (SAME TOP-DOWN VIEW)
# ============================================================================
print "=== Generating Image 4: R702W LRR + Bufadienolide ==="

# Same structure, same view - mutation is in HD2 not LRR
ray 6000, 4500
png img4_r702w_lrr_bufadienolide.png

print "Saved: img4_r702w_lrr_bufadienolide.png"

# ============================================================================
# IMAGE 5: WT LRR + FEBUXOSTAT - TOP-DOWN POCKET VIEW
# ============================================================================
print "=== Generating Image 5: WT LRR + Febuxostat (Top-Down Pocket) ==="

delete bufadienolide

load C:/Users/vasud/nod2-screening-data/PHASE_6/structures/febuxostat_docked.sdf, febuxostat

# Show ligand prominently
show sticks, febuxostat
show spheres, febuxostat
set sphere_scale, 0.25, febuxostat
set stick_radius, 0.2, febuxostat
color ligand_orange, febuxostat

# TOP-DOWN VIEW
center febuxostat
orient febuxostat
turn x, -90
zoom febuxostat, 8

ray 6000, 4500
png img5_wt_lrr_febuxostat.png

print "Saved: img5_wt_lrr_febuxostat.png"

# ============================================================================
# IMAGE 6: R702W LRR + FEBUXOSTAT (SAME TOP-DOWN VIEW)
# ============================================================================
print "=== Generating Image 6: R702W LRR + Febuxostat ==="

ray 6000, 4500
png img6_r702w_lrr_febuxostat.png

print "Saved: img6_r702w_lrr_febuxostat.png"

# ============================================================================
# BONUS: FULL STRUCTURE OVERVIEW WITH ANNOTATIONS
# ============================================================================
print "=== Generating Bonus: Full Structure Overview ==="

delete all
load C:/Users/vasud/nod2-screening-data/PHASE_6/structures/NOD2_alphafold_full.pdb, nod2_overview

hide all
show cartoon, nod2_overview

# Color domains
select card1, nod2_overview and resi 1-93
select card2, nod2_overview and resi 104-191
select nbd, nod2_overview and resi 273-577
select hd1, nod2_overview and resi 578-628
select hd2, nod2_overview and resi 629-743
select lrr, nod2_overview and resi 744-1040

color domain_gray, card1
color domain_gray, card2
color domain_marine, nbd
color domain_cyan, hd1
color domain_orange, hd2
color domain_lime, lrr

# Mutation site
select r702, nod2_overview and resi 702
show spheres, r702
set sphere_scale, 0.6, r702
color mutation_red, r702

# Binding pocket
select pocket, nod2_overview and resi 1007+1008+1011+1034+1037
show spheres, pocket
set sphere_scale, 0.5, pocket
color pocket_magenta, pocket

# Distance line
distance allosteric_dist3, r702 and name CA, pocket and resi 1008 and name CA
hide labels, allosteric_dist3
set dash_color, yellow, allosteric_dist3
set dash_width, 3, allosteric_dist3

# Clean side view
orient nod2_overview
turn x, 90
turn y, -30
zoom nod2_overview, 0

ray 6000, 4500
png img_bonus_overview.png

print "Saved: img_bonus_overview.png"

# ============================================================================
# SUMMARY
# ============================================================================
print ""
print "=============================================="
print "ALL PATHWAY DIAGRAM IMAGES GENERATED"
print "=============================================="
print "Resolution: 6000 x 4500 (6K)"
print ""
print "IMPROVED VIEWS:"
print "  Full structures: Clean side view, domains stacked vertically"
print "  Pocket views: Top-down 'bowl' view with transparent surface"
print "  Distance line: Yellow dashed line shows ~80A allosteric distance"
print ""
print "Files created:"
print "  1. img1_wt_nod2_full.png       - WT full (side view + distance)"
print "  2. img2_r702w_nod2_full.png    - R702W full (same angle)"
print "  3. img3_wt_lrr_bufadienolide.png   - WT pocket (top-down bowl)"
print "  4. img4_r702w_lrr_bufadienolide.png - R702W pocket (same angle)"
print "  5. img5_wt_lrr_febuxostat.png  - WT pocket (top-down bowl)"
print "  6. img6_r702w_lrr_febuxostat.png - R702W pocket (same angle)"
print "  7. img_bonus_overview.png      - Full overview"
print "=============================================="
