# PyMOL Visualization Script for NOD2-Febuxostat Complex
# ISEF 2026 - NOD2-CROHN Drug Discovery Project
# Author: Tanmay Vasudeva

# Load structures
load nod2_alphafold.pdb, NOD2
load febuxostat_docked.pdb, febuxostat

# Basic setup
bg_color white
set ray_opaque_background, 1
set antialias, 2
set ray_trace_mode, 1

# Color scheme
color marine, NOD2
color gray80, NOD2 and name C*

# Highlight domains
# HD2 domain (contains R702 mutation site)
select HD2, resi 539-744
color lightblue, HD2

# LRR domain (drug binding region)
select LRR, resi 745-1040
color palegreen, LRR

# Mutation site R702
select mut_site, resi 702
show spheres, mut_site
color red, mut_site

# Binding pocket residues
select pocket, resi 1007+1008+1009+1010+1011+1012+1013+1014+1035+1037
show sticks, pocket
color yellow, pocket and name C*

# Febuxostat ligand
show sticks, febuxostat
color orange, febuxostat and name C*
set stick_radius, 0.2, febuxostat

# Surface for binding pocket
select pocket_surface, resi 1004-1040
show surface, pocket_surface
set transparency, 0.7, pocket_surface
color white, pocket_surface

# Labels
# Uncomment to add labels
# label resi 702 and name CA, "R702W"
# label resi 1008 and name CA, "E1008"

# Camera settings for different views
# View 1: Overall structure
orient
zoom NOD2, 10
ray 2400, 2400
png figure_nod2_overall.png, dpi=300

# View 2: Binding pocket close-up
zoom pocket, 15
turn y, 30
ray 2400, 2400
png figure_binding_pocket.png, dpi=300

# View 3: Mutation site
zoom mut_site, 20
turn x, 15
ray 2400, 2400
png figure_mutation_site.png, dpi=300

# Session save
save nod2_febuxostat_session.pse

print("PyMOL visualization complete!")
print("Generated: figure_nod2_overall.png")
print("Generated: figure_binding_pocket.png")
print("Generated: figure_mutation_site.png")
print("Session saved: nod2_febuxostat_session.pse")
