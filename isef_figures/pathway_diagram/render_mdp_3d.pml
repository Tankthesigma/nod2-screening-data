# ============================================================================
# MDP 3D Render - Muramyl Dipeptide (bacterial trigger for NOD2)
# ============================================================================

reinitialize

# Quality settings
set antialias, 2
set ray_trace_mode, 1
set ray_shadows, 1
set ray_trace_fog, 0
set depth_cue, 0
set ambient, 0.4
set spec_reflect, 0.5
set spec_power, 200
set ray_opaque_background, 0
bg_color white

# Load MDP structure
load C:/Users/vasud/nod2-screening-data/isef_figures/pathway_diagram/MDP.sdf, mdp

# Show as ball-and-stick (spheres + sticks)
hide everything
show sticks, mdp
show spheres, mdp

# Set sphere size for ball-and-stick look
set sphere_scale, 0.25, mdp
set stick_radius, 0.15, mdp

# Color by element - orange carbons to match pathway theme
color orange, mdp and elem C
color red, mdp and elem O
color blue, mdp and elem N
color white, mdp and elem H

# Orient and zoom
orient mdp
zoom mdp, 2

# Render
ray 1200, 1200
png C:/Users/vasud/nod2-screening-data/isef_figures/pathway_diagram/MDP_3D_ISEF.png, dpi=300
