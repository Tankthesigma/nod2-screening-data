# PyMOL Script for ISEF 2026 NOD2 Project
# Author: Tanmay Vasudeva
# Texas Virtual Academy at Hallsville
#
# Run: pymol isef_pymol_script.pml

# Load NOD2 structure
load ../PHASE_6/structures/NOD2_alphafold.pdb, NOD2

# Basic settings for publication quality
bg_color white
set ray_opaque_background, 1
set antialias, 2
set ray_trace_mode, 1
set ray_shadows, 0
set depth_cue, 0

# Color by domain
# CARD1: 1-93 (gray)
select CARD1, resi 1-93
color gray60, CARD1

# CARD2: 104-191 (gray)
select CARD2, resi 104-191
color gray60, CARD2

# NBD: 273-577 (marine blue)
select NBD, resi 273-577
color marine, NBD

# HD1: 578-628 (cyan)
select HD1, resi 578-628
color cyan, HD1

# HD2: 629-743 (orange - contains R702W mutation)
select HD2, resi 629-743
color orange, HD2

# LRR: 744-1040 (forest green - ligand binding)
select LRR, resi 744-1040
color forest, LRR

# Highlight R702 mutation site
select mut_site, resi 702
show spheres, mut_site
color red, mut_site
set sphere_scale, 1.5, mut_site

# Key binding pocket residues
select pocket_residues, resi 1007+1008+1011+1034+1037
show sticks, pocket_residues
color yellow, pocket_residues and name C*
set stick_radius, 0.25, pocket_residues

# Show cartoon for overall structure
hide everything
show cartoon, NOD2
show spheres, mut_site
show sticks, pocket_residues

# ======= VIEW 1: Overall Structure =======
orient
zoom NOD2, 10
set_view (\
     0.9,    0.0,    0.4,\
     0.0,    1.0,    0.0,\
    -0.4,    0.0,    0.9,\
     0.0,    0.0, -300.0,\
     0.0,    0.0,    0.0,\
   200.0,  400.0,  -20.0 )

# Add labels for domains
# pseudoatom label_hd2, pos=[X, Y, Z]  # Position near HD2
# label label_hd2, "HD2 (R702W)"

ray 2400, 2400
png figure6_nod2_structure.png, dpi=300

# ======= VIEW 2: Mutation Site Close-up =======
zoom mut_site, 25
turn y, 30
ray 2400, 2400
png figure_mutation_site.png, dpi=300

# ======= VIEW 3: Binding Pocket =======
zoom pocket_residues, 20
turn x, 15
turn y, -20
ray 2400, 2400
png figure7_binding_pocket.png, dpi=300

# ======= VIEW 4: Distance measurement =======
# Show distance from R702 to binding pocket
orient
zoom NOD2, 10
distance dist_mut_pocket, resi 702 and name CA, resi 1008 and name CA
hide labels, dist_mut_pocket
set dash_color, red, dist_mut_pocket
set dash_width, 3, dist_mut_pocket

ray 2400, 2400
png figure_distance.png, dpi=300

# Save session
save nod2_isef_session.pse

print("PyMOL figures generated:")
print("  - figure6_nod2_structure.png")
print("  - figure_mutation_site.png")
print("  - figure7_binding_pocket.png")
print("  - figure_distance.png")
print("  - nod2_isef_session.pse (session file)")
