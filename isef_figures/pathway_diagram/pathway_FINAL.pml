# ============================================================================
# ISEF 2026 - PUBLICATION POCKET IMAGES
# 4 images: Pocket surface + ligand spheres
# ============================================================================

reinitialize

# Quality settings
set ray_trace_mode, 1
set antialias, 3
set ray_shadows, 1
set line_smooth, 1
set sphere_mode, 5
set ambient, 0.4
set specular, 0.3
set ray_opaque_background, 1
set orthoscopic, 1
set depth_cue, 0
set fog, 0
set surface_quality, 2
set sphere_quality, 2
set solvent_radius, 1.4

bg_color white
cd C:/Users/vasud/nod2-screening-data/isef_figures/pathway_diagram

# ============================================================================
# LOAD STRUCTURES
# ============================================================================
load C:/Users/vasud/nod2-screening-data/PHASE_6/structures/NOD2_LRR_clean.pdb, protein
load C:/Users/vasud/nod2-screening-data/PHASE_6/structures/natural_top_docked.sdf, bufa
load C:/Users/vasud/nod2-screening-data/PHASE_6/structures/febuxostat_docked.sdf, febu

hide all

# ============================================================================
# IMAGE 1 & 2: BUFADIENOLIDE (CYAN)
# ============================================================================
print "=== Rendering Bufadienolide images ==="

# Disable febuxostat
disable febu

# Select pocket residues around bufadienolide (5A)
select pocket_bufa, byres (protein within 5 of bufa)

# Show pocket surface
show surface, pocket_bufa
set surface_color, gray90, pocket_bufa
set transparency, 0.3, pocket_bufa

# Show bufadienolide as CYAN SPHERES
show spheres, bufa
set sphere_scale, 1.0, bufa
color cyan, bufa

# Camera setup
center bufa
orient bufa
turn y, 20
zoom bufa, 4

# Save this view matrix
get_view

ray 1500, 1500
png pocket_bufa_wt.png
print "Saved: pocket_bufa_wt.png"

# Image 2 is identical (same LRR pocket)
ray 1500, 1500
png pocket_bufa_r702w.png
print "Saved: pocket_bufa_r702w.png"

# ============================================================================
# IMAGE 3 & 4: FEBUXOSTAT (ORANGE)
# ============================================================================
print "=== Rendering Febuxostat images ==="

# Hide bufadienolide surface and spheres
hide surface, pocket_bufa
hide spheres, bufa
disable bufa

# Enable febuxostat
enable febu

# Select pocket residues around febuxostat (5A)
select pocket_febu, byres (protein within 5 of febu)

# Show pocket surface
show surface, pocket_febu
set surface_color, gray90, pocket_febu
set transparency, 0.3, pocket_febu

# Show febuxostat as ORANGE SPHERES
show spheres, febu
set sphere_scale, 1.0, febu
color orange, febu

# Camera setup
center febu
orient febu
turn y, 20
zoom febu, 4

ray 1500, 1500
png pocket_febu_wt.png
print "Saved: pocket_febu_wt.png"

# Image 4 is identical (same LRR pocket)
ray 1500, 1500
png pocket_febu_r702w.png
print "Saved: pocket_febu_r702w.png"

# ============================================================================
# DONE
# ============================================================================
print ""
print "=============================================="
print "4 POCKET IMAGES COMPLETE"
print "=============================================="
print "  pocket_bufa_wt.png      - Bufadienolide (cyan)"
print "  pocket_bufa_r702w.png   - Bufadienolide (cyan)"
print "  pocket_febu_wt.png      - Febuxostat (orange)"
print "  pocket_febu_r702w.png   - Febuxostat (orange)"
print ""
print "Features:"
print "  - Pocket surface (5A around ligand)"
print "  - Opaque white background"
print "  - Ligand as bright spheres (scale 1.0)"
print "  - 1500x1500 resolution"
print "=============================================="
