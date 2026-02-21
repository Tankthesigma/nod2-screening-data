# ============================================================
# ISEF 2026 - MAXIMUM QUALITY PYMOL SCRIPT
# 5K resolution, highest ray tracing, publication-grade
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
# MAXIMUM QUALITY SETTINGS
# ============================================================
bg_color white
set orthoscopic, 1
set depth_cue, 0

# Ray tracing - HIGHEST quality
set ray_trace_mode, 1
set ray_trace_color, black
set ray_shadow, 1
set ray_shadows, 1
set ray_opaque_background, 0
set ray_trace_fog, 0
set ray_trace_gain, 0.12
set ray_trace_disco_factor, 0
set ray_trace_frames, 1
set ray_trace_thickness, 0.2
set ray_trace_dither, 0
set ray_trace_antialias, 2

# Antialiasing - MAXIMUM
set antialias, 4
set smooth_lines, 1
set line_smooth, 1
set stick_quality, 24
set sphere_quality, 4
set cartoon_sampling, 16
set cartoon_smooth_loops, 1
set cartoon_fancy_helices, 1
set cartoon_fancy_sheets, 1
set cartoon_flat_sheets, 0
set cartoon_transparency, 0

# Lighting - OPTIMAL for print
set ambient, 0.25
set direct, 0.75
set reflect, 0.10
set shininess, 10
set specular, 0.25
set specular_intensity, 0.3
set light_count, 6
set light2, 1
set light3, 1
set light4, 1
set light5, 1
set light6, 1

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

# Binding pocket - LARGE MAGENTA SPHERES (VISIBLE!)
show spheres, binding_pocket
color magenta, binding_pocket
set sphere_scale, 1.8, binding_pocket

orient nod2_full
zoom nod2_full, 5

# 5K OUTPUT - HIGHEST RESOLUTION
ray 5120, 3840
png C:/Users/vasud/nod2-screening-data/isef_figures/v5_winner/ppt_photos/fig1_MAXQUALITY_full.png, dpi=300

print "SAVED: fig1_MAXQUALITY_full.png (5K)"

# ============================================================
# FIGURE 2: R702W MUTATION CLOSEUP
# ============================================================
hide everything

select R702_region, nod2_full and byres (resi 702 around 12)
show cartoon, R702_region
color 0xff8c00, R702_region

# R702 as thick sticks
show sticks, R702_site
color red, R702_site and elem C
color 0xff6666, R702_site and elem O
color 0x6666ff, R702_site and elem N
set stick_radius, 0.30, R702_site

# Nearby residues
select nearby_702, byres (R702_site around 5) and not R702_site
show sticks, nearby_702
color gray50, nearby_702 and elem C
util.cnc nearby_702
set stick_radius, 0.15, nearby_702

center R702_site
zoom R702_site, 10

# 5K OUTPUT
ray 5120, 3840
png C:/Users/vasud/nod2-screening-data/isef_figures/v5_winner/ppt_photos/fig2_MAXQUALITY_r702w.png, dpi=300

print "SAVED: fig2_MAXQUALITY_r702w.png (5K)"

# ============================================================
print ""
print "============================================"
print "MAXIMUM QUALITY RENDERS COMPLETE!"
print ""
print "Resolution: 5120 x 3840 (5K)"
print "Antialias: 4 (maximum)"
print "Ray trace antialias: 2"
print "Sphere quality: 4"
print "Cartoon sampling: 16"
print ""
print "Files saved to ppt_photos/"
print "============================================"
