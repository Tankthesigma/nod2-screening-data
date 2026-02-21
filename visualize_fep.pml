# PyMOL visualization script for NOD2 FEP analysis
# Run with: pymol visualize_fep.pml

# Load the WT complex
load C:/Users/vasud/nod2-screening-data/fep_complete/fep_pmx/wt_complex/topology.pdb, wt_complex

# Hide everything first
hide all

# NOD2 Domain Boundaries (residue ranges)
# CARD1: 1-93
# CARD2: 104-191
# NBD:   273-577
# HD1:   578-628
# HD2:   629-743  <-- R702W is here
# LRR:   744-1040 <-- Ligand binds here

# Color protein by domain
select card1, resi 1-93 and chain A
select card2, resi 104-191 and chain A
select nbd, resi 273-577 and chain A
select hd1, resi 578-628 and chain A
select hd2, resi 629-743 and chain A
select lrr, resi 744-1040 and chain A

color gray50, card1
color gray60, card2
color marine, nbd
color cyan, hd1
color orange, hd2
color forest, lrr

# Show cartoon representation
show cartoon, chain A

# Highlight R702 (mutation site in HD2)
select res702, resi 702 and chain A
show sticks, res702
color red, res702
label res702 and name CA, "R702"

# Show LRR binding pocket residues
select pocket, resi 1007+1008+1009+1010+1011+1034+1035+1036+1037 and chain A
show sticks, pocket
color yellow, pocket

# Show ligand (UNK residue)
select ligand, resn UNK
show sticks, ligand
color magenta, ligand

# Add distance measurement between R702 and ligand
distance r702_to_ligand, res702 and name CA, ligand

# Center view on ligand
center ligand
zoom ligand, 30

# Set nice visualization parameters
set cartoon_transparency, 0.3
set stick_radius, 0.2
bg_color white
set ray_shadows, 0

# Create a legend (group objects)
group domains, card1 card2 nbd hd1 hd2 lrr

# Print summary
print "=== NOD2 FEP Visualization ==="
print "Domain colors:"
print "  CARD1 (1-93): gray"
print "  CARD2 (104-191): gray"
print "  NBD (273-577): marine blue"
print "  HD1 (578-628): cyan"
print "  HD2 (629-743): ORANGE - contains R702W mutation"
print "  LRR (744-1040): forest green - contains ligand binding site"
print ""
print "Key features:"
print "  R702 (mutation site): RED sticks"
print "  LRR pocket residues: YELLOW sticks"
print "  Ligand: MAGENTA sticks"
print ""
print "Distance from R702 to ligand: ~80 Angstroms"
print "R702W is in HD2 domain, NOT in LRR binding site!"
