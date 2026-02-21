# PyMOL Script: Transparent Background Renders
# Output: transparent_bg/ folder

# Global settings for ALL renders
set ray_opaque_background, 0
set ray_shadows, 0
set antialias, 2
set ray_trace_mode, 1

# ==============================================================================
# FIGURE A1: Full NOD2 Structure with Domains
# ==============================================================================
reinitialize
set ray_opaque_background, 0

load C:/Users/vasud/nod2-screening-data/PHASE_6/structures/NOD2_alphafold_full.pdb, nod2_full

# Domain coloring
color gray80, nod2_full
select card1, resi 1-93
color gray60, card1
select card2, resi 104-191
color gray60, card2
select nbd, resi 273-577
color marine, nbd
select hd1, resi 578-628
color cyan, hd1
select hd2, resi 629-743
color orange, hd2
select lrr, resi 744-1040
color forest, lrr

show cartoon, nod2_full
orient nod2_full

ray 2400, 2400
png C:/Users/vasud/nod2-screening-data/isef_figures/v5_winner/ppt_photos/transparent_bg/figA1_transparent.png

# ==============================================================================
# FIGURE E3: 80A Distance (R702W to pocket)
# ==============================================================================
reinitialize
set ray_opaque_background, 0
set dash_gap, 0.3
set dash_width, 3
set dash_color, red

load C:/Users/vasud/nod2-screening-data/PHASE_6/structures/NOD2_alphafold_full.pdb, nod2

color gray70, nod2
select hd2, resi 629-743
color orange, hd2
select lrr, resi 744-1040
color forest, lrr

show cartoon, nod2

select r702, resi 702
show spheres, r702 and name CA
color red, r702

select pocket, resi 1008 and name CA
show spheres, pocket
color cyan, pocket

distance dist_line, r702 and name CA, resi 1008 and name CA
hide labels, dist_line
color red, dist_line

orient nod2
turn y, 30

ray 2400, 2400
png C:/Users/vasud/nod2-screening-data/isef_figures/v5_winner/ppt_photos/transparent_bg/figE3_transparent.png

# ==============================================================================
# FIGURE E2: R702W Mutation Site Closeup
# ==============================================================================
reinitialize
set ray_opaque_background, 0

load C:/Users/vasud/nod2-screening-data/PHASE_6/structures/NOD2_alphafold_full.pdb, nod2_full

color gray80, nod2_full
select hd2, resi 629-743
color orange, hd2
select lrr, resi 744-1040
color forest, lrr

select r702_site, resi 702
show sticks, r702_site
color red, r702_site
show spheres, r702_site and name CA
set sphere_scale, 0.8, r702_site

show cartoon, nod2_full

center r702_site
zoom r702_site, 15

ray 2400, 2400
png C:/Users/vasud/nod2-screening-data/isef_figures/v5_winner/ppt_photos/transparent_bg/figE2_transparent.png

# ==============================================================================
# FIGURE E4: Binding Pocket with Surface
# ==============================================================================
reinitialize
set ray_opaque_background, 0
set transparency, 0.3

load C:/Users/vasud/nod2-screening-data/PHASE_6/structures/complex_febuxostat.pdb, complex
remove resn HOH or resn NA or resn CL

select protein_sel, polymer
select ligand_sel, organic

color palecyan, protein_sel
show surface, protein_sel
set surface_color, palecyan, protein_sel

show sticks, ligand_sel
color orange, ligand_sel
set stick_radius, 0.25, ligand_sel

select key_res, resi 1008+1037+1011+1010+1007+1014 and protein_sel
show sticks, key_res
color yellow, key_res and elem C
color red, key_res and elem O
color blue, key_res and elem N

center ligand_sel
zoom ligand_sel, 8

ray 2400, 2400
png C:/Users/vasud/nod2-screening-data/isef_figures/v5_winner/ppt_photos/transparent_bg/figE4_transparent.png

# ==============================================================================
# FIGURE B1: Ligand in Pocket with H-bonds
# ==============================================================================
reinitialize
set ray_opaque_background, 0
set h_bond_cutoff_center, 3.6
set h_bond_cutoff_edge, 3.2

load C:/Users/vasud/nod2-screening-data/PHASE_6/structures/complex_febuxostat.pdb, complex
remove resn HOH or resn NA or resn CL

select protein_sel, polymer
select ligand_sel, organic

color palegreen, protein_sel
show cartoon, protein_sel

show sticks, ligand_sel
color orange, ligand_sel
set stick_radius, 0.2

select pocket_res, byres ligand_sel around 5 and polymer
show sticks, pocket_res
color lightteal, pocket_res and elem C

distance hbonds, ligand_sel, pocket_res, 3.5, mode=2
color yellow, hbonds
hide labels, hbonds

center ligand_sel
zoom ligand_sel, 6

ray 2400, 2400
png C:/Users/vasud/nod2-screening-data/isef_figures/v5_winner/ppt_photos/transparent_bg/figB1_transparent.png

# ==============================================================================
# DONE
# ==============================================================================
print ""
print "============================================"
print "TRANSPARENT PYMOL RENDERS COMPLETE:"
print "  figA1_transparent.png - Full structure"
print "  figE3_transparent.png - 80A distance"
print "  figE2_transparent.png - R702W closeup"
print "  figE4_transparent.png - Binding pocket"
print "  figB1_transparent.png - H-bonds"
print "============================================"
print "Saved to: transparent_bg/"
print ""

quit
