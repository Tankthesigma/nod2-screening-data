# ============================================================================
# ISEF 2026 - NOD2 Pathway Diagram PyMOL Script (10K Resolution)
# Tanmay Vasudeva - Texas Virtual Academy at Hallsville
# ============================================================================
#
# This script generates 6 high-resolution images for the pathway diagram:
# 1. img1_wt_nod2_full.png - WT NOD2 full structure with domain coloring
# 2. img2_r702w_nod2_full.png - R702W NOD2 full structure with mutation marker
# 3. img3_wt_lrr_bufadienolide.png - WT LRR pocket with Bufadienolide
# 4. img4_r702w_lrr_bufadienolide.png - R702W LRR pocket with Bufadienolide
# 5. img5_wt_lrr_febuxostat.png - WT LRR pocket with Febuxostat
# 6. img6_r702w_lrr_febuxostat.png - R702W LRR pocket with Febuxostat
#
# Resolution: 10000x7500 (10K) for maximum quality
# Run in PyMOL: @pathway_diagram_10K.pml
# ============================================================================

# ============================================================================
# IMAGE 1: WT NOD2 Full Structure
# ============================================================================
reinitialize

# QUALITY SETTINGS (must be AFTER reinitialize)
set antialias, 2
set ray_trace_mode, 1
set ray_shadows, 1
set ray_trace_fog, 0
set depth_cue, 0
set cartoon_fancy_helices, 1
set cartoon_smooth_loops, 1
set cartoon_sampling, 20
set ambient, 0.4
set spec_reflect, 0.5
set spec_power, 200
set ray_opaque_background, 0
bg_color white

# Load WT structure
load C:/Users/vasud/nod2-screening-data/PHASE_6/structures/NOD2_alphafold_full.pdb, wt_nod2

# Show as cartoon
hide everything
show cartoon, wt_nod2

# Color domains
color gray70, wt_nod2 and resi 1-93
color gray70, wt_nod2 and resi 104-191
color marine, wt_nod2 and resi 273-577
color cyan, wt_nod2 and resi 578-628
color orange, wt_nod2 and resi 629-743
color forest, wt_nod2 and resi 744-1040

# Mark R702 location (WT - show as reference)
select r702_site, wt_nod2 and resi 702
show spheres, r702_site and name CA
color yellow, r702_site and name CA
set sphere_scale, 1.5, r702_site

# Mark binding pocket residues
select pocket_residues, wt_nod2 and resi 1007+1008+1011+1034+1037
show spheres, pocket_residues and name CA
color magenta, pocket_residues and name CA
set sphere_scale, 0.8, pocket_residues

# Set view - full structure, slightly rotated
orient wt_nod2
turn y, 30
turn x, 10
zoom wt_nod2, 5

# Render at 10K resolution
ray 10000, 7500
png C:/Users/vasud/nod2-screening-data/isef_figures/pathway_diagram/img1_wt_nod2_full.png, dpi=300

# ============================================================================
# IMAGE 2: R702W NOD2 Full Structure
# ============================================================================
reinitialize
set antialias, 2
set ray_trace_mode, 1
set ray_shadows, 1
set ray_trace_fog, 0
set depth_cue, 0
set cartoon_fancy_helices, 1
set cartoon_smooth_loops, 1
set cartoon_sampling, 20
set ambient, 0.4
set ray_opaque_background, 0
bg_color white

# Load R702W structure
load C:/Users/vasud/nod2-screening-data/PHASE_A1_mutant_MD/structures/NOD2_R702W.pdb, r702w_nod2

# Show as cartoon
hide everything
show cartoon, r702w_nod2

# Color domains (same as WT)
color gray70, r702w_nod2 and resi 1-93
color gray70, r702w_nod2 and resi 104-191
color marine, r702w_nod2 and resi 273-577
color cyan, r702w_nod2 and resi 578-628
color orange, r702w_nod2 and resi 629-743
color forest, r702w_nod2 and resi 744-1040

# Mark W702 mutation site with RED sphere (mutant)
select w702_site, r702w_nod2 and resi 702
show spheres, w702_site and name CA
color red, w702_site and name CA
set sphere_scale, 2.0, w702_site

# Mark binding pocket residues
select pocket_residues, r702w_nod2 and resi 1007+1008+1011+1034+1037
show spheres, pocket_residues and name CA
color magenta, pocket_residues and name CA
set sphere_scale, 0.8, pocket_residues

# Set SAME view as WT for direct comparison
orient r702w_nod2
turn y, 30
turn x, 10
zoom r702w_nod2, 5

# Render at 10K resolution
ray 10000, 7500
png C:/Users/vasud/nod2-screening-data/isef_figures/pathway_diagram/img2_r702w_nod2_full.png, dpi=300

# ============================================================================
# IMAGE 3: WT LRR + Bufadienolide
# ============================================================================
reinitialize
set antialias, 2
set ray_trace_mode, 1
set ray_shadows, 1
set ray_trace_fog, 0
set depth_cue, 0
set cartoon_fancy_helices, 1
set cartoon_smooth_loops, 1
set cartoon_sampling, 20
set ambient, 0.4
set ray_opaque_background, 0
bg_color white

# Load WT structure and ligand
load C:/Users/vasud/nod2-screening-data/PHASE_6/structures/NOD2_alphafold_full.pdb, wt_nod2
load C:/Users/vasud/nod2-screening-data/PHASE_6/structures/natural_top_docked.sdf, bufadienolide

# Show LRR domain only
hide everything
show cartoon, wt_nod2 and resi 744-1040
color forest, wt_nod2 and resi 744-1040

# Show key binding residues as sticks
select key_residues, wt_nod2 and resi 1007+1008+1011+1034+1037
show sticks, key_residues
color yellow, key_residues and elem C
util.cnc key_residues

# Show ligand
show sticks, bufadienolide
color cyan, bufadienolide and elem C
util.cnc bufadienolide

# Zoom to pocket
select pocket_region, wt_nod2 and resi 1000-1040
orient pocket_region
zoom pocket_region, 10

# Render at 10K resolution
ray 10000, 7500
png C:/Users/vasud/nod2-screening-data/isef_figures/pathway_diagram/img3_wt_lrr_bufadienolide.png, dpi=300

# ============================================================================
# IMAGE 4: R702W LRR + Bufadienolide
# ============================================================================
reinitialize
set antialias, 2
set ray_trace_mode, 1
set ray_shadows, 1
set ray_trace_fog, 0
set depth_cue, 0
set cartoon_fancy_helices, 1
set cartoon_smooth_loops, 1
set cartoon_sampling, 20
set ambient, 0.4
set ray_opaque_background, 0
bg_color white

# Load R702W structure and ligand
load C:/Users/vasud/nod2-screening-data/PHASE_A1_mutant_MD/structures/NOD2_R702W.pdb, r702w_nod2
load C:/Users/vasud/nod2-screening-data/PHASE_6/structures/natural_top_docked.sdf, bufadienolide

# Show LRR domain only
hide everything
show cartoon, r702w_nod2 and resi 744-1040
color forest, r702w_nod2 and resi 744-1040

# Show key binding residues as sticks
select key_residues, r702w_nod2 and resi 1007+1008+1011+1034+1037
show sticks, key_residues
color yellow, key_residues and elem C
util.cnc key_residues

# Show ligand
show sticks, bufadienolide
color cyan, bufadienolide and elem C
util.cnc bufadienolide

# SAME zoom as WT for direct comparison
select pocket_region, r702w_nod2 and resi 1000-1040
orient pocket_region
zoom pocket_region, 10

# Render at 10K resolution
ray 10000, 7500
png C:/Users/vasud/nod2-screening-data/isef_figures/pathway_diagram/img4_r702w_lrr_bufadienolide.png, dpi=300

# ============================================================================
# IMAGE 5: WT LRR + Febuxostat
# ============================================================================
reinitialize
set antialias, 2
set ray_trace_mode, 1
set ray_shadows, 1
set ray_trace_fog, 0
set depth_cue, 0
set cartoon_fancy_helices, 1
set cartoon_smooth_loops, 1
set cartoon_sampling, 20
set ambient, 0.4
set ray_opaque_background, 0
bg_color white

# Load WT structure and ligand
load C:/Users/vasud/nod2-screening-data/PHASE_6/structures/NOD2_alphafold_full.pdb, wt_nod2
load C:/Users/vasud/nod2-screening-data/PHASE_6/structures/febuxostat_docked.sdf, febuxostat

# Show LRR domain only
hide everything
show cartoon, wt_nod2 and resi 744-1040
color forest, wt_nod2 and resi 744-1040

# Show key binding residues as sticks
select key_residues, wt_nod2 and resi 1007+1008+1011+1034+1037
show sticks, key_residues
color yellow, key_residues and elem C
util.cnc key_residues

# Show ligand
show sticks, febuxostat
color orange, febuxostat and elem C
util.cnc febuxostat

# Zoom to pocket
select pocket_region, wt_nod2 and resi 1000-1040
orient pocket_region
zoom pocket_region, 10

# Render at 10K resolution
ray 10000, 7500
png C:/Users/vasud/nod2-screening-data/isef_figures/pathway_diagram/img5_wt_lrr_febuxostat.png, dpi=300

# ============================================================================
# IMAGE 6: R702W LRR + Febuxostat
# ============================================================================
reinitialize
set antialias, 2
set ray_trace_mode, 1
set ray_shadows, 1
set ray_trace_fog, 0
set depth_cue, 0
set cartoon_fancy_helices, 1
set cartoon_smooth_loops, 1
set cartoon_sampling, 20
set ambient, 0.4
set ray_opaque_background, 0
bg_color white

# Load R702W structure and ligand
load C:/Users/vasud/nod2-screening-data/PHASE_A1_mutant_MD/structures/NOD2_R702W.pdb, r702w_nod2
load C:/Users/vasud/nod2-screening-data/PHASE_6/structures/febuxostat_docked.sdf, febuxostat

# Show LRR domain only
hide everything
show cartoon, r702w_nod2 and resi 744-1040
color forest, r702w_nod2 and resi 744-1040

# Show key binding residues as sticks
select key_residues, r702w_nod2 and resi 1007+1008+1011+1034+1037
show sticks, key_residues
color yellow, key_residues and elem C
util.cnc key_residues

# Show ligand
show sticks, febuxostat
color orange, febuxostat and elem C
util.cnc febuxostat

# SAME zoom as WT for direct comparison
select pocket_region, r702w_nod2 and resi 1000-1040
orient pocket_region
zoom pocket_region, 10

# Render at 10K resolution
ray 10000, 7500
png C:/Users/vasud/nod2-screening-data/isef_figures/pathway_diagram/img6_r702w_lrr_febuxostat.png, dpi=300

# ============================================================================
# DONE - Run this script in PyMOL: @pathway_diagram_10K.pml
# ============================================================================
