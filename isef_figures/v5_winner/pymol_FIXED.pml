# ============================================================
# ISEF 2026 - FIXED PYMOL SCRIPT
# Structure: complex_febuxostat.pdb (LRR domain 744-1040 + Febuxostat)
# NOTE: R702W mutation is NOT in this structure (HD2 domain not included)
# ============================================================

reinitialize

# Load the complex structure
load C:/Users/vasud/nod2-screening-data/MAC_ALL_FILES/md_simulation/structures/complex_febuxostat.pdb, complex

# ============================================================
# STEP 1: CREATE PROPER SELECTIONS
# ============================================================

# Protein = chain A (LRR domain, residues 744-1040)
select protein, chain A

# Ligand = Febuxostat specifically (chain L, resn FEB)
select ligand, resn FEB

# Key binding residues from PLIP analysis
select binding_res, resi 1007+1008+1011+1034+1036+1037

# Active site = residues within 5A of ligand, EXCLUDING ligand itself
select active_site, (byres chain A within 5 of ligand) and not ligand

# ============================================================
# STEP 2: GLOBAL SETTINGS FOR PUBLICATION QUALITY
# ============================================================

# High quality rendering
set antialias, 2
set ray_trace_mode, 1
set ray_shadows, 1
set ray_opaque_background, 1
set orthoscopic, 1

# Stick and sphere sizes
set stick_radius, 0.2
set sphere_scale, 0.3

# Cartoon settings
set cartoon_fancy_helices, 1
set cartoon_smooth_loops, 1
set cartoon_tube_radius, 0.6

# Label settings
set label_size, 22
set label_color, black
set label_font_id, 7

# Dash settings for polar contacts
set dash_width, 2.5
set dash_gap, 0.4
set dash_radius, 0.15

# ============================================================
# FIGURE 1: FULL STRUCTURE OVERVIEW (WHITE BG)
# ============================================================
bg_color white

hide everything

# Protein as light cyan cartoon
show cartoon, protein
color cyan, protein

# Ligand as magenta sticks (high contrast)
show sticks, ligand
color magenta, ligand and elem C
color red, ligand and elem O
color blue, ligand and elem N
color yellow, ligand and elem S
set stick_radius, 0.25, ligand

orient complex
zoom complex, 5

scene overview, store

ray 2400, 1800
png C:/Users/vasud/nod2-screening-data/isef_figures/v5_winner/fig1_overview_white.png, dpi=300
print "Saved: fig1_overview_white.png"

# ============================================================
# FIGURE 2: BINDING SITE WITH POLAR CONTACTS (KEY FIGURE!)
# ============================================================
hide everything
bg_color white

# Protein cartoon - light gray (context)
show cartoon, protein
color gray80, protein
set cartoon_transparency, 0.5

# Active site residues as sticks
show sticks, active_site
color gray50, active_site and elem C
util.cnc active_site
set stick_radius, 0.18, active_site

# Binding residues highlighted in green
color green, binding_res and elem C

# Ligand as orange sticks (pops!)
show sticks, ligand
color orange, ligand and elem C
color red, ligand and elem O
color blue, ligand and elem N
color yellow, ligand and elem S
set stick_radius, 0.22, ligand

# POLAR CONTACTS (H-bonds as yellow dashes)
distance hbonds, ligand, active_site, mode=2
color yellow, hbonds
hide labels, hbonds

center ligand
zoom ligand, 6

scene binding_hbonds, store

ray 2400, 1800
png C:/Users/vasud/nod2-screening-data/isef_figures/v5_winner/fig2_binding_hbonds.png, dpi=300
print "Saved: fig2_binding_hbonds.png"

# ============================================================
# FIGURE 3: BINDING SITE - NO LABELS (for PowerPoint labeling)
# ============================================================
# Same view, no PyMOL labels - add labels in PowerPoint!
hide labels
delete hbonds

ray 2400, 1800
png C:/Users/vasud/nod2-screening-data/isef_figures/v5_winner/fig3_binding_clean.png, dpi=300
print "Saved: fig3_binding_clean.png"

# ============================================================
# FIGURE 4: SURFACE + LIGAND IN POCKET
# ============================================================
hide everything
bg_color white

# Protein surface - white, semi-transparent
show surface, protein
color white, protein
set transparency, 0.7, protein

# Ligand as spheres INSIDE pocket
show spheres, ligand
color orange, ligand and elem C
color red, ligand and elem O
color blue, ligand and elem N
set sphere_scale, 0.4, ligand

center ligand
zoom ligand, 10

scene surface, store

ray 2400, 1800
png C:/Users/vasud/nod2-screening-data/isef_figures/v5_winner/fig4_pocket_surface.png, dpi=300
print "Saved: fig4_pocket_surface.png"

# ============================================================
# FIGURE 5: DARK BACKGROUND VERSION (POSTER CONTRAST)
# ============================================================
hide everything
bg_color black

# Protein cartoon - blue tint
show cartoon, protein
color marine, protein
set cartoon_transparency, 0.3

# Active site residues - BRIGHT GREEN (visible on black!)
show sticks, active_site
color tv_green, active_site and elem C
color red, active_site and elem O
color tv_blue, active_site and elem N
set stick_radius, 0.18, active_site

# Ligand - BRIGHT ORANGE (pops on black)
show sticks, ligand
color tv_orange, ligand and elem C
color red, ligand and elem O
color tv_blue, ligand and elem N
set stick_radius, 0.25, ligand

# Polar contacts - YELLOW dashes
distance hbonds_dark, ligand, active_site, mode=2
color yellow, hbonds_dark
hide labels, hbonds_dark

center ligand
zoom ligand, 6

scene dark_binding, store

ray 2400, 1800
png C:/Users/vasud/nod2-screening-data/isef_figures/v5_winner/fig5_binding_dark.png, dpi=300
print "Saved: fig5_binding_dark.png"

# ============================================================
# FIGURE 6: KEY RESIDUES CLOSEUP (PLIP RESIDUES)
# ============================================================
hide everything
delete hbonds_dark
bg_color white

# Just show binding residues and ligand
show sticks, binding_res
color tv_green, binding_res and elem C
util.cnc binding_res
set stick_radius, 0.2, binding_res

show sticks, ligand
color orange, ligand and elem C
color red, ligand and elem O
color blue, ligand and elem N
set stick_radius, 0.25, ligand

# Polar contacts between ligand and key residues
distance key_contacts, ligand, binding_res, mode=2
color yellow, key_contacts
hide labels, key_contacts

center ligand
zoom ligand, 5

scene key_residues, store

ray 2400, 1800
png C:/Users/vasud/nod2-screening-data/isef_figures/v5_winner/fig6_key_residues.png, dpi=300
print "Saved: fig6_key_residues.png"

# ============================================================
# FIGURE 7: OVERVIEW WITH BINDING SITE HIGHLIGHTED
# ============================================================
hide everything
delete key_contacts
bg_color white

# Full protein cartoon
show cartoon, protein
color cyan, protein

# Highlight binding pocket region in green
color tv_green, binding_res

# Show ligand as spheres
show spheres, ligand
color orange, ligand
set sphere_scale, 0.5, ligand

orient protein
zoom protein, 3

scene full_highlighted, store

ray 2400, 1800
png C:/Users/vasud/nod2-screening-data/isef_figures/v5_winner/fig7_full_highlighted.png, dpi=300
print "Saved: fig7_full_highlighted.png"

# ============================================================
print ""
print "============================================"
print "FIXED PYMOL EXPORTS COMPLETE!"
print ""
print "Files saved to v5_winner folder:"
print "  fig1_overview_white.png   - Full LRR structure"
print "  fig2_binding_hbonds.png   - Binding site + H-bonds"
print "  fig3_binding_clean.png    - Binding site (no labels)"
print "  fig4_pocket_surface.png   - Surface with ligand"
print "  fig5_binding_dark.png     - Binding site (dark bg)"
print "  fig6_key_residues.png     - Key PLIP residues"
print "  fig7_full_highlighted.png - Full structure + pocket"
print ""
print "NOTE: R702W mutation is NOT in this structure!"
print "      This file only contains LRR domain (744-1040)"
print "      For R702W visualization, need full NOD2 structure"
print "============================================"
