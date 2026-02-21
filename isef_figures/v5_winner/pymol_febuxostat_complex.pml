# ============================================================
# ISEF 2026 - FEBUXOSTAT BINDING SITE (with ligand!)
# THE MONEY SHOT for your poster
# ============================================================

reinitialize

# Load the complex structure
load C:/Users/vasud/nod2-screening-data/MAC_ALL_FILES/md_simulation/structures/complex_febuxostat.pdb, complex

# Separate protein and ligand
select protein, polymer
select ligand, organic

# ============================================================
# FIGURE: Febuxostat in Binding Pocket - POSTER QUALITY
# ============================================================
hide everything

# Protein as white cartoon
show cartoon, protein
color white, protein

# LIGAND - Make it POP (cyan/yellow)
show sticks, ligand
color cyan, ligand and elem C
color red, ligand and elem O
color blue, ligand and elem N
color yellow, ligand and elem S
set stick_radius, 0.25, ligand

# Key binding residues (from your contact analysis)
select binding_res, resi 1007+1008+1011+1034+1037
show sticks, binding_res
color gray70, binding_res and elem C
util.cnc binding_res

# R702W mutation site
select R702, resi 702
show sticks, R702
color magenta, R702 and elem C
util.cnc R702

# Settings
bg_color white
set ray_shadows, 0
set specular, 0.1
set antialias, 2
set label_size, 20
set label_color, black
set label_font_id, 7

# LABELS
label ligand and name C1, "Febuxostat"
label resi 1008 and name CA, "E1008"
label resi 1037 and name CA, "R1037"
label resi 1011 and name CA, "D1011"
label resi 702 and name CA, "R702W"

# Zoom to binding site
center ligand
zoom ligand, 8

ray 2400, 1800
png C:/Users/vasud/nod2-screening-data/isef_figures/v5_winner/fig_febuxostat_binding.png, dpi=300
print "Saved: fig_febuxostat_binding.png"

# ============================================================
# FIGURE: Surface view with ligand visible
# ============================================================
hide everything
hide labels

# Protein surface (transparent)
show surface, protein
color white, protein
set transparency, 0.7, protein

# Ligand as spheres inside
show spheres, ligand
color cyan, ligand and elem C
color red, ligand and elem O
color blue, ligand and elem N
set sphere_scale, 0.3, ligand

# R702W on surface
color magenta, R702

orient complex
zoom ligand, 12

ray 2400, 1800
png C:/Users/vasud/nod2-screening-data/isef_figures/v5_winner/fig_febuxostat_surface.png, dpi=300
print "Saved: fig_febuxostat_surface.png"

# ============================================================
# FIGURE: Show potential H-bonds (dashed lines)
# ============================================================
hide everything
hide labels

show cartoon, protein
color gray90, protein
set cartoon_transparency, 0.5

show sticks, ligand
color cyan, ligand and elem C
set stick_radius, 0.2, ligand

show sticks, binding_res
color gray60, binding_res and elem C
util.cnc binding_res

# Find polar contacts (H-bonds)
distance hbonds, ligand, binding_res, mode=2
color yellow, hbonds
set dash_width, 3
set dash_gap, 0.3
hide labels, hbonds

# Labels
set label_size, 18
label ligand and name C1, "Febuxostat"
label resi 1008 and name CA, "E1008"
label resi 1037 and name CA, "R1037"

center ligand
zoom ligand, 6

ray 2400, 1800
png C:/Users/vasud/nod2-screening-data/isef_figures/v5_winner/fig_febuxostat_hbonds.png, dpi=300
print "Saved: fig_febuxostat_hbonds.png"

# ============================================================
# FIGURE: Final POSTER figure (big, clear, labeled)
# ============================================================
hide everything
delete hbonds

show cartoon, protein
color palegreen, protein

# Big ligand
show spheres, ligand
color orange, ligand and elem C
color red, ligand and elem O
color blue, ligand and elem N
set sphere_scale, 0.5, ligand

# R702W
show spheres, R702
color magenta, R702

# BIG labels
set label_size, 28
set label_outline_color, white
label ligand and name C1, "FEBUXOSTAT"
label R702 and name CA, "R702W"

orient complex
zoom ligand, 10

ray 2400, 1800
png C:/Users/vasud/nod2-screening-data/isef_figures/v5_winner/fig_poster_final.png, dpi=300
print "Saved: fig_poster_final.png"

# ============================================================
print ""
print "============================================"
print "FEBUXOSTAT BINDING FIGURES COMPLETE!"
print ""
print "  fig_febuxostat_binding.png - Sticks view with labels"
print "  fig_febuxostat_surface.png - Surface with ligand inside"
print "  fig_febuxostat_hbonds.png  - H-bond visualization"
print "  fig_poster_final.png       - Big poster figure"
print "============================================"
