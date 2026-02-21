# ============================================================
# ISEF 2026 - FINAL CLEAN EXPORTS (FIXED)
# SOLID backgrounds, high contrast, professional
# ============================================================

reinitialize

# Load complex
load C:/Users/vasud/nod2-screening-data/MAC_ALL_FILES/md_simulation/structures/complex_febuxostat.pdb, complex

select protein, polymer
select ligand, organic
select R702, resi 702
select binding_res, resi 1007+1008+1011+1034+1036+1037

# ============================================================
# PROFESSIONAL SETTINGS
# ============================================================
set cartoon_fancy_helices, 1
set cartoon_smooth_loops, 1
set cartoon_tube_radius, 0.8
set stick_radius, 0.25
set sphere_scale, 0.4
set antialias, 2
set spec_reflect, 0.3
set specular, 0.3
set ray_trace_mode, 1
set ray_shadows, 1

# CRITICAL: Solid background (NOT transparent)
set ray_opaque_background, 1

# ============================================================
# FIGURE 1: OVERVIEW - SOLID BLACK BG
# ============================================================
bg_color black

hide everything

# Protein cartoon - GOLDEN/TAN color (like ISEF winners)
show cartoon, protein
color wheat, protein

# Semi-transparent white surface
show surface, protein
color white, protein
set transparency, 0.85, protein

# R702W - BRIGHT MAGENTA spheres
show spheres, R702
color magenta, R702
set sphere_scale, 0.6, R702

# Ligand - BRIGHT ORANGE (better than yellow on black)
show sticks, ligand
color tv_orange, ligand and elem C
color red, ligand and elem O
color blue, ligand and elem N
color yellow, ligand and elem S
set stick_radius, 0.3, ligand

orient complex
zoom complex, 5

ray 2400, 1800
png C:/Users/vasud/nod2-screening-data/isef_figures/v5_winner/final_overview_black.png, dpi=300
print "Saved: final_overview_black.png"

# ============================================================
# FIGURE 2: BINDING SITE - SOLID BLACK BG, HIGH CONTRAST
# ============================================================
hide everything
bg_color black

# Protein cartoon - dim gray (context only)
show cartoon, protein
color gray40, protein
set cartoon_transparency, 0.6

# Binding residues - BRIGHT GREEN (visible on black)
show sticks, binding_res
color limegreen, binding_res and elem C
color red, binding_res and elem O
color tv_blue, binding_res and elem N
set stick_radius, 0.2, binding_res

# Ligand - HOT ORANGE (pops!)
show sticks, ligand
color tv_orange, ligand and elem C
color red, ligand and elem O
color tv_blue, ligand and elem N
set stick_radius, 0.3, ligand

# R702W - MAGENTA
show sticks, R702
color hotpink, R702 and elem C
set stick_radius, 0.25, R702

center ligand
zoom ligand, 5

ray 2400, 1800
png C:/Users/vasud/nod2-screening-data/isef_figures/v5_winner/final_binding_black.png, dpi=300
print "Saved: final_binding_black.png"

# ============================================================
# FIGURE 3: POCKET SURFACE - SOLID BLACK BG
# ============================================================
hide everything
bg_color black

# Dark surface
show surface, protein
color gray20, protein
set transparency, 0.5, protein

# Ligand as bright spheres
show spheres, ligand
color tv_orange, ligand and elem C
color red, ligand and elem O
set sphere_scale, 0.5, ligand

# Highlight binding pocket area
color deepteal, binding_res

center ligand
zoom ligand, 8

ray 2400, 1800
png C:/Users/vasud/nod2-screening-data/isef_figures/v5_winner/final_pocket_black.png, dpi=300
print "Saved: final_pocket_black.png"

# ============================================================
# FIGURE 4: WHITE BACKGROUND VERSION (for printing)
# ============================================================
hide everything
bg_color white

show cartoon, protein
color palecyan, protein

show sticks, ligand
color tv_orange, ligand and elem C
color red, ligand and elem O
color tv_blue, ligand and elem N
set stick_radius, 0.3, ligand

show sticks, binding_res
color gray50, binding_res and elem C
util.cnc binding_res
set stick_radius, 0.2, binding_res

show spheres, R702
color magenta, R702

center ligand
zoom ligand, 6

ray 2400, 1800
png C:/Users/vasud/nod2-screening-data/isef_figures/v5_winner/final_binding_white.png, dpi=300
print "Saved: final_binding_white.png"

# ============================================================
# FIGURE 5: RIBBON + SURFACE (like winning ISEF posters)
# ============================================================
hide everything
bg_color white

# Nice blue ribbon
show cartoon, protein
color lightblue, protein

# Transparent surface overlay
show surface, protein
color white, protein
set transparency, 0.8, protein

# Ligand visible inside
show sticks, ligand
color tv_orange, ligand and elem C
set stick_radius, 0.25, ligand

# R702W highlighted
show spheres, R702
color magenta, R702

orient complex
zoom complex, 5

ray 2400, 1800
png C:/Users/vasud/nod2-screening-data/isef_figures/v5_winner/final_ribbon_surface.png, dpi=300
print "Saved: final_ribbon_surface.png"

# ============================================================
# FIGURE 6: CLOSE-UP R702W MUTATION SITE
# ============================================================
hide everything
bg_color black

select R702_region, byres (R702 expand 10)
show cartoon, R702_region
color gray50, R702_region

show sticks, R702
color hotpink, R702 and elem C
util.cnc R702
set stick_radius, 0.3, R702

# Nearby residues
show sticks, byres (R702 around 5) and not R702
color gray70, byres (R702 around 5) and not R702

center R702
zoom R702, 6

ray 2400, 1800
png C:/Users/vasud/nod2-screening-data/isef_figures/v5_winner/final_r702w_closeup.png, dpi=300
print "Saved: final_r702w_closeup.png"

# ============================================================
print ""
print "============================================"
print "FINAL EXPORTS COMPLETE!"
print ""
print "Files saved to v5_winner folder:"
print "  final_overview_black.png    - Full structure (BLACK bg)"
print "  final_binding_black.png     - Binding site (BLACK bg)"
print "  final_pocket_black.png      - Pocket surface (BLACK bg)"
print "  final_binding_white.png     - Binding site (WHITE bg)"
print "  final_ribbon_surface.png    - Ribbon + surface"
print "  final_r702w_closeup.png     - R702W mutation detail"
print ""
print "COLOR KEY:"
print "  ORANGE = Febuxostat (ligand)"
print "  GREEN  = Binding residues"
print "  MAGENTA/PINK = R702W mutation"
print "============================================"
