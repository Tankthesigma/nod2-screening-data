# ============================================================================
# ISEF 2026 - NOD2 Signaling Pathway Protein Renders
# Tanmay Vasudeva - Texas Virtual Academy at Hallsville
# ============================================================================
#
# This script generates renders of downstream signaling proteins:
# - RIPK2 (PDB: 4C8B) - full and dim versions
# - NF-kB (PDB: 1VKX) - p50/p65 with DNA, full and dim versions
#
# Note: MDP molecule should be created separately with RDKit or ChemDraw
# - TAK1 (PDB: 5JH6) - TAK1-TAB1 complex, full and dim versions
#
# Run in PyMOL: @pathway_proteins.pml
# ============================================================================

# ============================================================================
# RIPK2 - Full (Active pathway - bright)
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
set ambient, 0.4
set spec_reflect, 0.5
set spec_power, 200
set ray_opaque_background, 0
bg_color white

# Fetch RIPK2 structure from PDB
fetch 4C8B, ripk2, async=0

# Clean up - remove waters and ligands
remove solvent
remove organic

# Show as cartoon
hide everything
show cartoon, ripk2

# Color cyan (cleaner look)
color cyan, ripk2

# Orient and zoom
orient ripk2
zoom ripk2, 5

# Render
ray 1200, 1200
png C:/Users/vasud/nod2-screening-data/isef_figures/pathway_diagram/RIPK2_full_ISEF.png, dpi=300

# ============================================================================
# RIPK2 - Dim (Inactive/R702W pathway - faded but VISIBLE)
# ============================================================================
# Lighter color, less transparency so it's still visible

set cartoon_transparency, 0.3, ripk2
color palecyan, ripk2

# Render dim version
ray 1200, 1200
png C:/Users/vasud/nod2-screening-data/isef_figures/pathway_diagram/RIPK2_dim_ISEF.png, dpi=300

# ============================================================================
# NF-kB - Full (Active pathway - bright)
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
set ambient, 0.4
set spec_reflect, 0.5
set spec_power, 200
set ray_opaque_background, 0
bg_color white

# Fetch NF-kB structure (p50-p65 heterodimer with DNA - ACTIVE state)
fetch 1VKX, nfkb, async=0

# Remove waters
remove solvent

# Show DNA as sticks (optional - keep for context)
show sticks, nfkb and resn DA+DT+DG+DC
color white, nfkb and resn DA+DT+DG+DC

# Show as cartoon
hide everything
show cartoon, nfkb

# Color purple
color purple, nfkb

# Remove DNA if desired (optional - keep for context)
# remove nfkb and chain C

# Orient and zoom
orient nfkb
zoom nfkb, 5

# Render
ray 1200, 1200
png C:/Users/vasud/nod2-screening-data/isef_figures/pathway_diagram/NFkB_full_ISEF.png, dpi=300

# ============================================================================
# NF-kB - Dim (Inactive/R702W pathway - faded but VISIBLE)
# ============================================================================

set cartoon_transparency, 0.3, nfkb
color lightpink, nfkb

# Render dim version
ray 1200, 1200
png C:/Users/vasud/nod2-screening-data/isef_figures/pathway_diagram/NFkB_dim_ISEF.png, dpi=300

# ============================================================================
# TAK1 - Full (Active pathway - bright)
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
set ambient, 0.4
set spec_reflect, 0.5
set spec_power, 200
set ray_opaque_background, 0
bg_color white

# Fetch TAK1 structure from PDB
fetch 5JH6, tak1, async=0

# Remove waters and ligands
remove solvent
remove organic

# Show as cartoon
hide everything
show cartoon, tak1

# Color gold/orange
color orange, tak1

# Orient and zoom
orient tak1
zoom tak1, 5

# Render
ray 1200, 1200
png C:/Users/vasud/nod2-screening-data/isef_figures/pathway_diagram/TAK1_full_ISEF.png, dpi=300

# ============================================================================
# TAK1 - Dim (Inactive/R702W pathway - faded but VISIBLE)
# ============================================================================

set cartoon_transparency, 0.3, tak1
color lightorange, tak1

# Render dim version
ray 1200, 1200
png C:/Users/vasud/nod2-screening-data/isef_figures/pathway_diagram/TAK1_dim_ISEF.png, dpi=300

# ============================================================================
# DONE
# ============================================================================
# Done - check output directory for images
