# ============================================================
# ISEF 2026 - ONLY 2 FIGURES
# 1. Full NOD2 with domains + visible pocket
# 2. R702W mutation closeup
# ============================================================

reinitialize

# Load structure
load C:/Users/vasud/nod2-screening-data/PHASE_6/structures/NOD2_alphafold_full.pdb, nod2_full

# Domain selections
select CARD1, nod2_full and resi 1-93
select CARD2, nod2_full and resi 104-191
select NBD, nod2_full and resi 273-577
select HD1, nod2_full and resi 578-628
select HD2, nod2_full and resi 629-743
select LRR, nod2_full and resi 744-1040
select R702_site, nod2_full and resi 702
select binding_pocket, nod2_full and resi 1007+1008+1011+1034+1036+1037

# High quality settings
set antialias, 2
set ray_trace_mode, 1
set ray_shadows, 1
set cartoon_fancy_helices, 1
set cartoon_smooth_loops, 1

# ============================================================
# FIGURE 1: FULL NOD2 STRUCTURE
# ============================================================
bg_color white
hide everything
show cartoon, nod2_full

# Domain colors
color gray70, CARD1
color gray70, CARD2
color 0x0066cc, NBD
color 0x00cccc, HD1
color 0xff8c00, HD2
color 0x32cd32, LRR

# R702 - LARGE RED SPHERES
show spheres, R702_site
color red, R702_site
set sphere_scale, 1.8, R702_site

# Binding pocket - LARGE MAGENTA SPHERES (VISIBLE!)
show spheres, binding_pocket
color magenta, binding_pocket
set sphere_scale, 1.5, binding_pocket

orient nod2_full
zoom nod2_full, 5

ray 3840, 2880
png C:/Users/vasud/nod2-screening-data/isef_figures/v5_winner/ppt_photos/fig1_full_structure.png, dpi=300
print "Saved: fig1_full_structure.png"

# ============================================================
# FIGURE 2: R702W MUTATION CLOSEUP
# ============================================================
hide everything
bg_color white

select R702_region, nod2_full and byres (resi 702 around 12)
show cartoon, R702_region
color 0xff8c00, R702_region

# R702 as sticks
show sticks, R702_site
color red, R702_site and elem C
color red, R702_site and elem O
color blue, R702_site and elem N
set stick_radius, 0.35, R702_site

# Nearby residues
select nearby_702, byres (R702_site around 5) and not R702_site
show sticks, nearby_702
color gray50, nearby_702 and elem C
util.cnc nearby_702
set stick_radius, 0.15, nearby_702

center R702_site
zoom R702_site, 8

ray 3840, 2880
png C:/Users/vasud/nod2-screening-data/isef_figures/v5_winner/ppt_photos/fig2_r702w_mutation.png, dpi=300
print "Saved: fig2_r702w_mutation.png"

# ============================================================
print ""
print "DONE! 2 FIGURES CREATED:"
print "  1. fig1_full_structure.png - Full NOD2 with BIG magenta pocket"
print "  2. fig2_r702w_mutation.png - R702W closeup"
print ""
print "Run in PyMOL!"
# ============================================================
