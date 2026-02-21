# WT vs R702W - FIXED ANGLE + HIGH QUALITY
# Mutation clearly visible!

reinitialize

set ray_trace_mode, 1
set antialias, 3
set ray_shadows, 1
set ray_opaque_background, 0
set orthoscopic, 1
set depth_cue, 0
set cartoon_sampling, 14
set ambient, 0.4
set specular, 0.3

bg_color white
cd C:/Users/vasud/nod2-screening-data/isef_figures/pathway_diagram

# ============================================================================
# IMAGE 1: WILD-TYPE NOD2
# ============================================================================
print "=== WT NOD2 ==="

load C:/Users/vasud/nod2-screening-data/PHASE_6/structures/NOD2_alphafold_full.pdb, wt

hide all
show cartoon, wt

# Domain colors
color gray, wt and resi 1-191
color blue, wt and resi 273-577
color teal, wt and resi 578-628
color orange, wt and resi 629-743
color green, wt and resi 744-1040

# EXACT ANGLE from interactive session - mutation visible!
set_view (\
    -0.760378778,   -0.241969660,   -0.602679789,\
    -0.310938239,    0.950324833,    0.010744601,\
     0.570164680,    0.195578158,   -0.797883332,\
     0.000000000,    0.000000000, -287.975341797,\
     4.254291534,   -6.865200043,   -0.104904175,\
  -3751.645263672, 4327.597167969,  -20.000000000 )

ray 3000, 3000
png wt_full_FINAL.png
print "Saved: wt_full_FINAL.png"

# ============================================================================
# IMAGE 2: R702W MUTANT
# ============================================================================
print "=== R702W NOD2 ==="

delete all
load C:/Users/vasud/nod2-screening-data/PHASE_A1_mutant_MD/structures/NOD2_R702W.pdb, r702w

hide all
show cartoon, r702w

# Same domain colors
color gray, r702w and resi 1-191
color blue, r702w and resi 273-577
color teal, r702w and resi 578-628
color orange, r702w and resi 629-743
color green, r702w and resi 744-1040

# SOLID MAGENTA MUTATION - NO GLOW (glow causes black render bug)
select mutation, r702w and resi 702
show spheres, mutation
set sphere_scale, 2.0, mutation
color magenta, mutation

# SAME EXACT ANGLE as WT - mutation visible!
set_view (\
    -0.760378778,   -0.241969660,   -0.602679789,\
    -0.310938239,    0.950324833,    0.010744601,\
     0.570164680,    0.195578158,   -0.797883332,\
     0.000000000,    0.000000000, -287.975341797,\
     4.254291534,   -6.865200043,   -0.104904175,\
  -3751.645263672, 4327.597167969,  -20.000000000 )

ray 3000, 3000
png r702w_full_FINAL.png
print "Saved: r702w_full_FINAL.png"

print ""
print "DONE!"
