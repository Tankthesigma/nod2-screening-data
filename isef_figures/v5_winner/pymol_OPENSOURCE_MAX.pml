# ============================================================
# ISEF 2026 - MAXIMUM QUALITY FOR OPEN-SOURCE PYMOL
# Codex-verified - All settings compatible!
# ============================================================

reinitialize
set auto_zoom, off
set retain_order, 1

# Load structure
load C:/Users/vasud/nod2-screening-data/PHASE_6/structures/NOD2_alphafold_full.pdb, nod2_full

# Domain selections
select CARD1, nod2_full and resi 1-93
select CARD2, nod2_full and resi 104-191
select NBD,  nod2_full and resi 273-577
select HD1,  nod2_full and resi 578-628
select HD2,  nod2_full and resi 629-743
select LRR,  nod2_full and resi 744-1040
select R702_site, nod2_full and resi 702
select binding_pocket, nod2_full and resi 1007+1008+1011+1034+1036+1037

# ============================================================
# MAXIMUM QUALITY SETTINGS (2025 BEST PRACTICES)
# Source: PyMOLWiki + compchems.com
# ============================================================
bg_color white
set orthoscopic, 1
set depth_cue, 0

# Ray tracing - mode 1 = color + black outline (publication standard)
set ray_trace_mode, 1
set ray_shadow, 1
set ray_opaque_background, 0
set ray_trace_fog, 0
set ray_trace_gain, 0.1

# Antialiasing - MAXIMUM (2 = best quality, 4x render time)
set antialias, 2

# Geometry quality - CRANKED UP
set stick_quality, 25
set sphere_quality, 4
set surface_quality, 2
set cartoon_sampling, 20
set cartoon_smooth_loops, 1
set cartoon_fancy_helices, 1
set cartoon_fancy_sheets, 1
set cartoon_flat_sheets, 0
set cartoon_loop_quality, 10

# Lighting - BRIGHT (so whites look white)
set ambient, 0.4
set direct, 0.6
set reflect, 0.3
set shininess, 10
set specular, 0.25
set spec_reflect, 0.3

# ============================================================
# FIGURE 1: FULL NOD2 STRUCTURE
# ============================================================
hide everything
show cartoon, nod2_full

# Domain colors (matching your key)
color gray70, CARD1
color gray70, CARD2
color 0x0066cc, NBD
color 0x00cccc, HD1
color 0xff8c00, HD2
color 0x32cd32, LRR

# R702W mutation - LARGE RED SPHERES
show spheres, R702_site
color red, R702_site
set sphere_scale, 2.0, R702_site

# Binding pocket - LARGE MAGENTA SPHERES
show spheres, binding_pocket
color magenta, binding_pocket
set sphere_scale, 1.8, binding_pocket

orient nod2_full
zoom nod2_full, 5

# 6K OUTPUT - MAXIMUM RESOLUTION (20" x 15" at 300dpi = poster size!)
ray 6000, 4500
png C:/Users/vasud/nod2-screening-data/isef_figures/v5_winner/ppt_photos/fig1_BEST_full.png, dpi=300

echo SAVED: fig1_BEST_full.png (6K - 20x15 inches at 300dpi)

# ============================================================
# FIGURE 2: R702W MUTATION CLOSEUP
# ============================================================
hide everything

# Show region around R702
select R702_region, nod2_full and byres (resi 702 around 12)
show cartoon, R702_region
color 0xff8c00, R702_region

# R702 as thick sticks - RED
show sticks, R702_site
color red, R702_site and elem C
color 0xff6666, R702_site and elem O
color 0x6666ff, R702_site and elem N
set stick_radius, 0.30, R702_site

# Nearby residues as thin sticks
select nearby_702, byres (R702_site around 5) and not R702_site
show sticks, nearby_702
color gray50, nearby_702 and elem C
util.cnc nearby_702
set stick_radius, 0.15, nearby_702

center R702_site
zoom R702_site, 10

# 6K OUTPUT
ray 6000, 4500
png C:/Users/vasud/nod2-screening-data/isef_figures/v5_winner/ppt_photos/fig2_BEST_r702w.png, dpi=300

echo SAVED: fig2_BEST_r702w.png (6K - 20x15 inches at 300dpi)

# ============================================================
echo
echo ============================================
echo MAXIMUM QUALITY RENDERS COMPLETE!
echo
echo Resolution: 6000 x 4500 (6K)
echo Print size: 20 x 15 inches at 300dpi
echo Ray trace mode: 1 (publication standard)
echo Antialias: 2 (maximum)
echo Cartoon sampling: 20
echo
echo Files saved to ppt_photos/
echo ============================================
