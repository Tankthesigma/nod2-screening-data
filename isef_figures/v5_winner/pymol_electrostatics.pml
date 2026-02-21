# PyMOL Script: Electrostatic Surface Comparison (WT vs R702W)
# ISEF 2026 - NOD2 Project
# Author: Tanmay Vasudeva
#
# NOTE: This shows vacuum electrostatics for visual comparison.
# R702W mutation is ~80A from binding site - differences are subtle/allosteric.
#
# Run: pymol pymol_electrostatics.pml

# Reset and configure
reinitialize
bg_color white
set ray_trace_mode, 1
set ray_shadows, 0
set antialias, 2
set hash_max, 300
set cartoon_fancy_helices, 1
set cartoon_smooth_loops, 1
set spec_reflect, 0.3
set ambient, 0.3

# ==============================================================================
# FIGURE E1: WT NOD2 Electrostatic Surface (LRR domain)
# ==============================================================================

# Load WT structure (LRR domain from complex)
load C:/Users/vasud/nod2-screening-data/PHASE_6/structures/complex_febuxostat.pdb, wt_complex
remove resn HOH or resn NA or resn CL
remove not polymer

# Color by electrostatics (vacuum)
util.protein_vacuum_esp all, mode=2

# Set surface
hide all
show surface, wt_complex

# Focus on LRR binding region
orient wt_complex
zoom wt_complex, 5

# Render
ray 2400, 2400
png C:/Users/vasud/nod2-screening-data/isef_figures/v5_winner/ppt_photos/figE1_wt_electrostatic.png, dpi=300

# ==============================================================================
# FIGURE E2: R702W Mutation Site Close-up
# ==============================================================================

reinitialize
bg_color white
set ray_trace_mode, 1
set ray_shadows, 0

# Load full NOD2 structure for R702W visualization
load C:/Users/vasud/nod2-screening-data/PHASE_6/structures/NOD2_alphafold_full.pdb, nod2_full

# Color by domain
color gray80, nod2_full
# HD2 domain (contains R702) - orange
select hd2, resi 629-743
color orange, hd2
# LRR domain (binding site) - forest green
select lrr, resi 744-1040
color forest, lrr

# Highlight R702 mutation site
select r702_site, resi 702
show sticks, r702_site
color red, r702_site
show spheres, r702_site and name CA
set sphere_scale, 0.8, r702_site

# Label
#label r702_site and name CA, "R702W"

# Show cartoon
show cartoon, nod2_full
hide surface

# Focus on mutation region
center r702_site
zoom r702_site, 15

# Render
ray 2400, 2400
png C:/Users/vasud/nod2-screening-data/isef_figures/v5_winner/ppt_photos/figE2_r702w_mutation_site.png, dpi=300

# ==============================================================================
# FIGURE E3: Distance Visualization (R702 to Binding Pocket)
# ==============================================================================

reinitialize
bg_color white
set ray_trace_mode, 1
set ray_shadows, 0
set dash_gap, 0.3
set dash_width, 3
set dash_color, red

# Load full structure
load C:/Users/vasud/nod2-screening-data/PHASE_6/structures/NOD2_alphafold_full.pdb, nod2

# Color scheme
color gray70, nod2
select hd2, resi 629-743
color orange, hd2
select lrr, resi 744-1040
color forest, lrr

# Show cartoon
show cartoon, nod2

# Highlight R702 (mutation site)
select r702, resi 702
show spheres, r702 and name CA
color red, r702

# Highlight binding pocket residues
select pocket, resi 1008+1037+1011 and name CA
show spheres, pocket
color cyan, pocket

# Draw distance line
distance dist_r702_pocket, r702 and name CA, resi 1008 and name CA
hide labels, dist_r702_pocket
color red, dist_r702_pocket

# Orient to show full protein
orient nod2
turn y, 30

# Render
ray 2400, 2400
png C:/Users/vasud/nod2-screening-data/isef_figures/v5_winner/ppt_photos/figE3_distance_r702_to_pocket.png, dpi=300

# ==============================================================================
# FIGURE E4: Binding Pocket Surface with Key Residues
# ==============================================================================

reinitialize
bg_color white
set ray_trace_mode, 1
set ray_shadows, 0
set transparency, 0.3

# Load complex with febuxostat
load C:/Users/vasud/nod2-screening-data/PHASE_6/structures/complex_febuxostat.pdb, complex

# Clean up
remove resn HOH or resn NA or resn CL

# Separate protein and ligand
select protein_sel, polymer
select ligand_sel, organic

# Color protein
color palecyan, protein_sel

# Show surface with transparency
show surface, protein_sel
set surface_color, palecyan, protein_sel

# Show ligand as sticks
show sticks, ligand_sel
color orange, ligand_sel
set stick_radius, 0.25, ligand_sel

# Highlight key binding residues
select key_res, resi 1008+1037+1011+1010+1007+1014 and protein_sel
show sticks, key_res
color yellow, key_res and elem C
color red, key_res and elem O
color blue, key_res and elem N

# Focus on binding site
center ligand_sel
zoom ligand_sel, 8

# Render
ray 2400, 2400
png C:/Users/vasud/nod2-screening-data/isef_figures/v5_winner/ppt_photos/figE4_binding_pocket_surface.png, dpi=300

# ==============================================================================
# DONE
# ==============================================================================

print ""
print "============================================================"
print "ELECTROSTATIC FIGURES GENERATED:"
print "  figE1_wt_electrostatic.png - WT vacuum electrostatics"
print "  figE2_r702w_mutation_site.png - R702W location"
print "  figE3_distance_r702_to_pocket.png - 80A distance"
print "  figE4_binding_pocket_surface.png - Binding site surface"
print "============================================================"
print ""
print "All saved to: ppt_photos/"
print ""

quit
