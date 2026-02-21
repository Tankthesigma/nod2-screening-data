# ============================================================
# ISEF 2026 - CLEAN EXPORTS (Add labels in PowerPoint!)
# Dark background + No PyMOL labels
# ============================================================

reinitialize

# Load complex
load C:/Users/vasud/nod2-screening-data/MAC_ALL_FILES/md_simulation/structures/complex_febuxostat.pdb, complex

select protein, polymer
select ligand, organic
select R702, resi 702
select binding_res, resi 1007+1008+1011+1034+1037

# ============================================================
# SETTINGS FOR PROFESSIONAL LOOK
# ============================================================
set cartoon_fancy_helices, 1
set cartoon_smooth_loops, 1
set cartoon_tube_radius, 0.8
set stick_radius, 0.2
set sphere_scale, 0.3
set ray_shadows, 1
set ray_shadow_decay_factor, 0.1
set antialias, 2
set spec_reflect, 0.5
set specular, 0.5

# ============================================================
# FIGURE 1: Overview - DARK BACKGROUND (like ISEF winners)
# ============================================================
bg_color black

hide everything
show cartoon, protein
color slate, protein

# Surface (semi-transparent white)
show surface, protein
set transparency, 0.8, protein
color white, protein and surface

# R702W in magenta spheres
show spheres, R702
color magenta, R702

# Ligand in YELLOW (visible on dark bg)
show sticks, ligand
color tv_yellow, ligand and elem C
color red, ligand and elem O
color blue, ligand and elem N
set stick_radius, 0.3, ligand

orient complex
zoom complex, 5

ray 2400, 1800
png C:/Users/vasud/nod2-screening-data/isef_figures/v5_winner/clean_overview_dark.png, dpi=300
print "Saved: clean_overview_dark.png"

# ============================================================
# FIGURE 2: Binding Site Close-up - DARK BG
# ============================================================
hide everything
bg_color black

# Cartoon context
show cartoon, protein
color gray50, protein
set cartoon_transparency, 0.7

# Binding residues - CYAN sticks
show sticks, binding_res
color cyan, binding_res and elem C
color red, binding_res and elem O
color blue, binding_res and elem N

# Ligand - YELLOW (pops on dark)
show sticks, ligand
color tv_yellow, ligand and elem C
set stick_radius, 0.25, ligand

# R702W - MAGENTA
show sticks, R702
color magenta, R702 and elem C

center ligand
zoom ligand, 5

ray 2400, 1800
png C:/Users/vasud/nod2-screening-data/isef_figures/v5_winner/clean_binding_dark.png, dpi=300
print "Saved: clean_binding_dark.png"

# ============================================================
# FIGURE 3: Surface Pocket - DARK BG
# ============================================================
hide everything
bg_color black

show surface, protein
color gray30, protein
set transparency, 0.6, protein

# Ligand as spheres INSIDE pocket
show spheres, ligand
color tv_yellow, ligand and elem C
color red, ligand and elem O
set sphere_scale, 0.4, ligand

# Highlight pocket region
color tv_blue, binding_res

center ligand
zoom ligand, 8

ray 2400, 1800
png C:/Users/vasud/nod2-screening-data/isef_figures/v5_winner/clean_pocket_dark.png, dpi=300
print "Saved: clean_pocket_dark.png"

# ============================================================
# FIGURE 4: Same views with WHITE background (for variety)
# ============================================================
bg_color white

hide everything
show cartoon, protein
color palecyan, protein

show sticks, ligand
color orange, ligand and elem C
set stick_radius, 0.25

show sticks, binding_res
color gray70, binding_res and elem C
util.cnc binding_res

show spheres, R702
color magenta, R702
set sphere_scale, 0.5, R702

center ligand
zoom ligand, 6

ray 2400, 1800
png C:/Users/vasud/nod2-screening-data/isef_figures/v5_winner/clean_binding_white.png, dpi=300
print "Saved: clean_binding_white.png"

# ============================================================
# FIGURE 5: Ribbon with surface overlay (like ISEF example)
# ============================================================
bg_color white

hide everything

# Cartoon in nice blue
show cartoon, protein
color lightblue, protein

# Semi-transparent surface
show surface, protein
color white, protein
set transparency, 0.75, protein

# Ligand
show sticks, ligand
color tv_yellow, ligand and elem C

# R702W
show spheres, R702
color magenta, R702

orient complex
zoom complex, 5

ray 2400, 1800
png C:/Users/vasud/nod2-screening-data/isef_figures/v5_winner/clean_ribbon_surface.png, dpi=300
print "Saved: clean_ribbon_surface.png"

# ============================================================
print ""
print "============================================"
print "CLEAN EXPORTS DONE!"
print ""
print "NOW ADD LABELS IN POWERPOINT:"
print "  1. Insert image"
print "  2. Insert > Shapes > Line with arrow"
print "  3. Add text boxes: 'R702W', 'Febuxostat', 'E1008', etc."
print "  4. Add caption: 'Image created by Tanmay Vasudeva'"
print "============================================"
