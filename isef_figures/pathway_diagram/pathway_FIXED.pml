# ISEF 2026 - NOD2 PATHWAY - FIXED VERSION
# Fixes: Mutation spheres visible, Febuxostat orange (not cyan)
reinitialize

# Settings
set ray_trace_mode, 1
set antialias, 2
set cartoon_sampling, 14
set orthoscopic, 1
set ray_opaque_background, 0

# Load
load C:/Users/vasud/nod2-screening-data/PHASE_6/structures/NOD2_alphafold_full.pdb, wt_nod2
load C:/Users/vasud/nod2-screening-data/PHASE_A1_mutant_MD/structures/NOD2_R702W.pdb, r702w_nod2
load C:/Users/vasud/nod2-screening-data/PHASE_6/structures/natural_top_docked.sdf, bufadienolide
load C:/Users/vasud/nod2-screening-data/PHASE_6/structures/febuxostat_docked.sdf, febuxostat

# Selections
select wt_card, wt_nod2 and resi 1-191
select wt_nbd, wt_nod2 and resi 273-577
select wt_hd1, wt_nod2 and resi 578-628
select wt_hd2, wt_nod2 and resi 629-743
select wt_lrr, wt_nod2 and resi 744-1040
select mut_card, r702w_nod2 and resi 1-191
select mut_nbd, r702w_nod2 and resi 273-577
select mut_hd1, r702w_nod2 and resi 578-628
select mut_hd2, r702w_nod2 and resi 629-743
select mut_lrr, r702w_nod2 and resi 744-1040
select w702, r702w_nod2 and resi 702
select wt_pocket, wt_nod2 and resi 1007+1008+1011+1034+1037
select mut_pocket, r702w_nod2 and resi 1007+1008+1011+1034+1037

# === IMG1: WT FULL ===
print "IMG1: WT Full"
hide everything
bg_color white
show cartoon, wt_nod2
color gray70, wt_card
color marine, wt_nbd
color cyan, wt_hd1
color tv_orange, wt_hd2
color lime, wt_lrr
orient wt_nod2
zoom wt_nod2, 3
ray 4000, 3000
png C:/Users/vasud/nod2-screening-data/isef_figures/pathway_diagram/img1_wt_full.png, dpi=300
print "DONE img1"

# === IMG2: R702W FULL WITH MUTATION HIGHLIGHTED ===
print "IMG2: R702W Full + Mutation"
hide everything
show cartoon, r702w_nod2
color gray70, mut_card
color marine, mut_nbd
color cyan, mut_hd1
color tv_orange, mut_hd2
color lime, mut_lrr
# BIG RED SPHERES for mutation
show spheres, w702
color red, w702
set sphere_scale, 2.0, w702
set sphere_transparency, 0, w702
# Same view
orient r702w_nod2
zoom r702w_nod2, 3
ray 4000, 3000
png C:/Users/vasud/nod2-screening-data/isef_figures/pathway_diagram/img2_r702w_full.png, dpi=300
print "DONE img2"

# === IMG3: WT + BUFADIENOLIDE (CYAN) ===
print "IMG3: WT + Bufadienolide"
hide everything
bg_color white
# Hide febuxostat completely
disable febuxostat
# Show WT LRR
show cartoon, wt_lrr
color lime, wt_lrr
set cartoon_transparency, 0.3, wt_lrr
# Pocket residues
show sticks, wt_pocket
color lime, wt_pocket and elem C
util.cnc wt_pocket
# Bufadienolide - CYAN
enable bufadienolide
show sticks, bufadienolide
color cyan, bufadienolide and elem C
util.cnc bufadienolide
set stick_radius, 0.25, bufadienolide
center bufadienolide
zoom bufadienolide, 8
ray 4000, 3000
png C:/Users/vasud/nod2-screening-data/isef_figures/pathway_diagram/img3_wt_bufa.png, dpi=300
print "DONE img3"

# === IMG4: R702W + BUFADIENOLIDE (CYAN) ===
print "IMG4: R702W + Bufadienolide"
hide everything
bg_color white
disable febuxostat
show cartoon, mut_lrr
color lime, mut_lrr
set cartoon_transparency, 0.3, mut_lrr
show sticks, mut_pocket
color lime, mut_pocket and elem C
util.cnc mut_pocket
enable bufadienolide
show sticks, bufadienolide
color cyan, bufadienolide and elem C
util.cnc bufadienolide
set stick_radius, 0.25, bufadienolide
center bufadienolide
zoom bufadienolide, 8
ray 4000, 3000
png C:/Users/vasud/nod2-screening-data/isef_figures/pathway_diagram/img4_r702w_bufa.png, dpi=300
print "DONE img4"

# === IMG5: WT + FEBUXOSTAT (ORANGE) ===
print "IMG5: WT + Febuxostat"
hide everything
bg_color white
# Hide bufadienolide, show febuxostat
disable bufadienolide
enable febuxostat
show cartoon, wt_lrr
color lime, wt_lrr
set cartoon_transparency, 0.3, wt_lrr
show sticks, wt_pocket
color lime, wt_pocket and elem C
util.cnc wt_pocket
# Febuxostat - ORANGE
show sticks, febuxostat
color tv_orange, febuxostat and elem C
color red, febuxostat and elem O
color blue, febuxostat and elem N
color yellow, febuxostat and elem S
set stick_radius, 0.25, febuxostat
center febuxostat
zoom febuxostat, 8
ray 4000, 3000
png C:/Users/vasud/nod2-screening-data/isef_figures/pathway_diagram/img5_wt_febu.png, dpi=300
print "DONE img5"

# === IMG6: R702W + FEBUXOSTAT (ORANGE) ===
print "IMG6: R702W + Febuxostat"
hide everything
bg_color white
disable bufadienolide
enable febuxostat
show cartoon, mut_lrr
color lime, mut_lrr
set cartoon_transparency, 0.3, mut_lrr
show sticks, mut_pocket
color lime, mut_pocket and elem C
util.cnc mut_pocket
show sticks, febuxostat
color tv_orange, febuxostat and elem C
color red, febuxostat and elem O
color blue, febuxostat and elem N
color yellow, febuxostat and elem S
set stick_radius, 0.25, febuxostat
center febuxostat
zoom febuxostat, 8
ray 4000, 3000
png C:/Users/vasud/nod2-screening-data/isef_figures/pathway_diagram/img6_r702w_febu.png, dpi=300
print "DONE img6"

print "=== ALL 6 IMAGES COMPLETE ==="
quit
