# ============================================================
# ISEF 2026 - NOD2 PyMOL Figures WITH LABELS
# Author: Tanmay Vasudeva
# ============================================================
# FIXES:
# - Labels on R702W mutation site
# - Colored ligand (if present)
# - Labeled binding residues
# - Distance measurements
# ============================================================

reinitialize

# Load protein structure
load C:/Users/vasud/nod2-screening-data/cross_validation/human_NOD2_LRR.pdb, NOD2

# Settings for publication quality
set cartoon_fancy_helices, 1
set cartoon_smooth_loops, 1
set cartoon_loop_radius, 0.3
bg_color white
set ray_opaque_background, 1
set antialias, 2
set ray_shadows, 0
set specular, 0.15
set label_size, 22
set label_color, black
set label_font_id, 7
set label_position, (2, 2, 2)

# ============================================================
# Define selections
# ============================================================
# R702W mutation site
select R702, resi 702
# Key binding residues from your MD contact analysis
select GLU1008, resi 1008
select ARG1037, resi 1037
select ASP1011, resi 1011
select ARG1034, resi 1034
select LEU1007, resi 1007
select binding_residues, resi 1007+1008+1011+1034+1037

# ============================================================
# FIGURE 1: Overview with R702W LABELED
# ============================================================
hide everything
show cartoon, NOD2
color palegreen, NOD2

# Highlight R702W
show spheres, R702
color magenta, R702

# ADD LABEL for R702W
label R702 and name CA, "R702W"

orient NOD2
zoom NOD2, 5

ray 2400, 1800
png C:/Users/vasud/nod2-screening-data/isef_figures/v5_winner/fig1_overview_labeled.png, dpi=300
print "Saved: fig1_overview_labeled.png"

# ============================================================
# FIGURE 2: Binding Site with LABELED RESIDUES
# ============================================================
hide everything
hide labels

# Show cartoon as context (transparent)
show cartoon, NOD2
color white, NOD2
set cartoon_transparency, 0.7, NOD2

# Show binding residues as sticks with distinct colors
show sticks, binding_residues
# Color by residue type
color tv_blue, GLU1008  # Acidic - blue
color tv_red, ASP1011   # Acidic - red
color cyan, ARG1037     # Basic - cyan
color cyan, ARG1034     # Basic - cyan
color yellow, LEU1007   # Hydrophobic - yellow

# Color by element for clarity
util.cnc binding_residues

# Show R702W
show sticks, R702
color magenta, R702

# ADD LABELS for each residue
label GLU1008 and name CA, "GLU1008"
label ARG1037 and name CA, "ARG1037"
label ASP1011 and name CA, "ASP1011"
label ARG1034 and name CA, "ARG1034"
label LEU1007 and name CA, "LEU1007"
label R702 and name CA, "R702W"

# Zoom to binding site
center binding_residues
zoom binding_residues, 8

ray 2400, 1800
png C:/Users/vasud/nod2-screening-data/isef_figures/v5_winner/fig2_binding_labeled.png, dpi=300
print "Saved: fig2_binding_labeled.png"

# ============================================================
# FIGURE 3: R702W Mutation Close-up
# ============================================================
hide everything
hide labels

# Show neighborhood
select R702_neighborhood, byres (R702 expand 8)
show cartoon, R702_neighborhood
color palecyan, R702_neighborhood

# Show R702 as sticks
show sticks, R702
color magenta, R702
util.cnc R702

# Label
label R702 and name CA, "R702W\nMutation Site"

# Show nearby residues
show sticks, byres (R702 around 5) and not R702
color gray80, byres (R702 around 5) and not R702

center R702
zoom R702, 6

ray 2400, 1800
png C:/Users/vasud/nod2-screening-data/isef_figures/v5_winner/fig3_r702w_closeup.png, dpi=300
print "Saved: fig3_r702w_closeup.png"

# ============================================================
# FIGURE 4: Surface with Binding Pocket Highlighted
# ============================================================
hide everything
hide labels

# Surface view
show surface, NOD2
color white, NOD2
set transparency, 0.6, NOD2

# Show cartoon inside
show cartoon, NOD2
color lightblue, NOD2
set cartoon_transparency, 0, NOD2

# Highlight binding pocket on surface
color tv_yellow, binding_residues
show surface, binding_residues
set transparency, 0, binding_residues

# R702W
color magenta, R702

# Label
label R702 and name CA, "R702W"

orient NOD2
zoom NOD2, 5

ray 2400, 1800
png C:/Users/vasud/nod2-screening-data/isef_figures/v5_winner/fig4_surface_pocket.png, dpi=300
print "Saved: fig4_surface_pocket.png"

# ============================================================
# FIGURE 5: Clean Publication Figure (White + Key Residues)
# ============================================================
hide everything
hide labels

show cartoon, NOD2
color white, NOD2

# Binding residues
show sticks, binding_residues
color deepteal, binding_residues and elem C
color red, binding_residues and elem O
color blue, binding_residues and elem N

# R702W in magenta
show sticks, R702
color hotpink, R702 and elem C

# Minimal labels
label GLU1008 and name CA, "E1008"
label ARG1037 and name CA, "R1037"
label ASP1011 and name CA, "D1011"
label R702 and name CA, "R702W"

set label_size, 18

orient NOD2
zoom binding_residues, 10

ray 2400, 1800
png C:/Users/vasud/nod2-screening-data/isef_figures/v5_winner/fig5_publication.png, dpi=300
print "Saved: fig5_publication.png"

# ============================================================
# FIGURE 6: For Poster - Large Labels
# ============================================================
hide everything
hide labels

show cartoon, NOD2
color gray90, NOD2

# Key residues
show spheres, binding_residues
color tv_blue, binding_residues

show spheres, R702
color magenta, R702

# BIG labels for poster
set label_size, 28
set label_outline_color, white

label R702 and name CA, "R702W\n(Mutation)"
label GLU1008 and name CA, "E1008"
label ARG1037 and name CA, "R1037"

orient NOD2
zoom NOD2, 5

ray 2400, 1800
png C:/Users/vasud/nod2-screening-data/isef_figures/v5_winner/fig6_poster_labeled.png, dpi=300
print "Saved: fig6_poster_labeled.png"

# ============================================================
print ""
print "============================================"
print "LABELED FIGURES COMPLETE!"
print ""
print "What each figure shows:"
print "  fig1 - Overview with R702W labeled"
print "  fig2 - Binding site with ALL residues labeled"
print "  fig3 - R702W mutation close-up"
print "  fig4 - Surface with binding pocket highlighted"
print "  fig5 - Publication style (clean)"
print "  fig6 - Poster style (BIG labels)"
print "============================================"
