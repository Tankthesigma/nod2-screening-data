# R702W ONLY - FIXED MAGENTA GLOW
reinitialize

set ray_trace_mode, 1
set antialias, 3
set ray_opaque_background, 0
set orthoscopic, 1
set depth_cue, 0
set cartoon_sampling, 14

# FIX FOR MAGENTA GLOW (prevents black absorption)
set ray_transparency_mode, 3
set ray_shadow, off
set ray_trace_fog, off
set two_sided_lighting, on
set ambient, 0.6
set specular, 0.2

bg_color white
cd C:/Users/vasud/nod2-screening-data/isef_figures/pathway_diagram

print "=== R702W NOD2 with MAGENTA GLOW ==="

load C:/Users/vasud/nod2-screening-data/PHASE_A1_mutant_MD/structures/NOD2_R702W.pdb, r702w

hide all
show cartoon, r702w

# Domain colors
color gray, r702w and resi 1-191
color blue, r702w and resi 273-577
color teal, r702w and resi 578-628
color orange, r702w and resi 629-743
color green, r702w and resi 744-1040

# MAGENTA MUTATION with GLOW
select mutation, r702w and resi 702
show spheres, mutation
set sphere_scale, 2.0, mutation
color magenta, mutation

# Glow halo - now renders MAGENTA not black!
create glow, mutation
show spheres, glow
set sphere_scale, 2.8, glow
set sphere_transparency, 0.65, glow
color magenta, glow

# EXACT ANGLE from interactive session
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
print "DONE!"
