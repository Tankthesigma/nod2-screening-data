# ============================================================
# ISEF 2026 - IMPROVED POSTER FIGURES
# Better quality, larger markers, 4K resolution
# ============================================================

reinitialize

# ============================================================
# LOAD STRUCTURES
# ============================================================

load C:/Users/vasud/nod2-screening-data/PHASE_6/structures/NOD2_alphafold_full.pdb, nod2_full

# Domain selections
select CARD1, nod2_full and resi 1-93
select CARD2, nod2_full and resi 104-191
select NBD, nod2_full and resi 273-577
select HD1, nod2_full and resi 578-628
select HD2, nod2_full and resi 629-743
select LRR, nod2_full and resi 744-1040

# R702 mutation site
select R702_site, nod2_full and resi 702

# Binding pocket - key residues
select binding_pocket, nod2_full and resi 1007+1008+1011+1034+1036+1037

# ============================================================
# HIGH QUALITY SETTINGS
# ============================================================
set antialias, 2
set ray_trace_mode, 1
set ray_shadows, 1
set ray_trace_gain, 0.1
set ray_trace_fog, 0
set depth_cue, 0
set specular, 0.3
set cartoon_fancy_helices, 1
set cartoon_smooth_loops, 1
set cartoon_sampling, 14
set ribbon_sampling, 10

# ============================================================
# FIGURE 1: FULL NOD2 WITH VISIBLE POCKET (WHITE BG)
# ============================================================
bg_color white
hide everything

show cartoon, nod2_full

# Domain colors (matching key)
color gray70, CARD1
color gray70, CARD2
color 0x0066cc, NBD
color 0x00cccc, HD1
color 0xff8c00, HD2
color 0x32cd32, LRR

# R702 - LARGE RED SPHERES
show spheres, R702_site
color red, R702_site
set sphere_scale, 1.5, R702_site

# Binding pocket - LARGE MAGENTA SPHERES (MORE VISIBLE!)
show spheres, binding_pocket
color magenta, binding_pocket
set sphere_scale, 1.2, binding_pocket

# Add distance line between R702 and pocket
distance dist_allosteric, R702_site and name CA, binding_pocket and resi 1008 and name CA
color yellow, dist_allosteric
set dash_width, 4
set dash_gap, 0.3
hide labels, dist_allosteric

orient nod2_full
zoom nod2_full, 5

# Ray trace at 4K
ray 3840, 2880
png C:/Users/vasud/nod2-screening-data/isef_figures/v5_winner/ppt_photos/figA1_full_domains_IMPROVED.png, dpi=300
print "Saved: figA1_full_domains_IMPROVED.png"

# ============================================================
# FIGURE 2: FULL NOD2 WITH DISTANCE LABEL
# ============================================================
# Show distance label
set label_size, 24
set label_color, black
label dist_allosteric, "79.4 A"

ray 3840, 2880
png C:/Users/vasud/nod2-screening-data/isef_figures/v5_winner/ppt_photos/figA1_full_domains_LABELED.png, dpi=300
print "Saved: figA1_full_domains_LABELED.png"

hide labels

# ============================================================
# FIGURE 3: R702W MUTATION CLOSEUP (IMPROVED)
# ============================================================
hide everything
bg_color white

# Show region around R702
select R702_region, nod2_full and byres (resi 702 around 15)
show cartoon, R702_region
color 0xff8c00, R702_region

# R702 as sticks - BRIGHT RED
show sticks, R702_site
color red, R702_site and elem C
color 0xff6666, R702_site and elem O
color 0x6666ff, R702_site and elem N
set stick_radius, 0.35, R702_site

# Also show as transparent sphere for emphasis
show spheres, R702_site
set sphere_transparency, 0.6, R702_site
set sphere_scale, 1.0, R702_site

# Nearby residues as thin sticks
select nearby_702, byres (R702_site around 6) and not R702_site
show sticks, nearby_702
color gray50, nearby_702 and elem C
util.cnc nearby_702
set stick_radius, 0.15, nearby_702

center R702_site
zoom R702_site, 10

ray 3840, 2880
png C:/Users/vasud/nod2-screening-data/isef_figures/v5_winner/ppt_photos/figE2_r702w_IMPROVED.png, dpi=300
print "Saved: figE2_r702w_IMPROVED.png"

# ============================================================
# FIGURE 4: BINDING POCKET CLOSEUP
# ============================================================
hide everything
bg_color white

# Show LRR around pocket
select pocket_region, nod2_full and byres (binding_pocket around 12)
show cartoon, pocket_region
color 0x32cd32, pocket_region
set cartoon_transparency, 0.3

# Pocket residues as sticks
show sticks, binding_pocket
color magenta, binding_pocket and elem C
util.cnc binding_pocket
set stick_radius, 0.25, binding_pocket

# Labels for key residues
set label_size, 18
set label_color, 0x1e3a5f
label binding_pocket and name CA and resi 1008, "GLU1008"
label binding_pocket and name CA and resi 1011, "ASP1011"
label binding_pocket and name CA and resi 1037, "ARG1037"

center binding_pocket
zoom binding_pocket, 8

ray 3840, 2880
png C:/Users/vasud/nod2-screening-data/isef_figures/v5_winner/ppt_photos/figE4_pocket_IMPROVED.png, dpi=300
print "Saved: figE4_pocket_IMPROVED.png"

hide labels

# ============================================================
# FIGURE 5: SIDE-BY-SIDE VIEW (R702 and Pocket visible)
# ============================================================
hide everything
bg_color white

show cartoon, nod2_full

# Domain colors
color gray70, CARD1
color gray70, CARD2
color 0x0066cc, NBD
color 0x00cccc, HD1
color 0xff8c00, HD2
color 0x32cd32, LRR

# Both markers visible
show spheres, R702_site
color red, R702_site
set sphere_scale, 2.0, R702_site

show spheres, binding_pocket
color magenta, binding_pocket
set sphere_scale, 1.5, binding_pocket

# Distance line
show dashes, dist_allosteric
color yellow, dist_allosteric

# Rotate for best view of both sites
turn y, 45
turn x, -15

zoom nod2_full, 5

ray 3840, 2880
png C:/Users/vasud/nod2-screening-data/isef_figures/v5_winner/ppt_photos/figA1_full_ROTATED.png, dpi=300
print "Saved: figA1_full_ROTATED.png"

# ============================================================
print ""
print "============================================"
print "IMPROVED POSTER FIGURES DONE!"
print ""
print "NEW FILES:"
print "  figA1_full_domains_IMPROVED.png  - 4K, larger spheres"
print "  figA1_full_domains_LABELED.png   - With 79.4A label"
print "  figE2_r702w_IMPROVED.png         - Better mutation closeup"
print "  figE4_pocket_IMPROVED.png        - Pocket with labels"
print "  figA1_full_ROTATED.png           - Rotated view"
print ""
print "SPHERE SIZES:"
print "  R702W = 1.5-2.0 (LARGE RED)"
print "  Pocket = 1.2-1.5 (LARGE MAGENTA)"
print ""
print "Run this in PyMOL to generate!"
print "============================================"
