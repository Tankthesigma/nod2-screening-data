# ============================================================
# ISEF 2026 - PUBLICATION QUALITY PYMOL FIGURES
# Based on proper PyMOL techniques (polar contacts, selections, contrast)
# ============================================================

reinitialize

# Load the complex structure
load C:/Users/vasud/nod2-screening-data/MAC_ALL_FILES/md_simulation/structures/complex_febuxostat.pdb, complex

# ============================================================
# STEP 1: CREATE PROPER SELECTIONS
# ============================================================

# Separate protein and ligand
select protein, polymer
select ligand, organic

# R702W mutation site (HD2 domain)
select R702, resi 702

# Key binding residues from PLIP analysis
select binding_res, resi 1007+1008+1011+1034+1036+1037

# Active site = all residues within 5A of ligand (like tutorial)
select active_site, byres all within 5 of ligand

# Water molecules near ligand (if any exist in structure)
select ligand_water, resn HOH within 3.2 of ligand

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

# Label settings (for measurements)
set label_size, 22
set label_color, black
set label_font_id, 7
set label_outline_color, white

# Dash settings for polar contacts
set dash_width, 2.5
set dash_gap, 0.4
set dash_radius, 0.15

# ============================================================
# FIGURE 1: OVERVIEW - WHITE BACKGROUND (PUBLICATION STYLE)
# ============================================================
bg_color white

hide everything

# Protein as cyan cartoon (good contrast)
show cartoon, protein
color palecyan, protein

# Ligand as MAGENTA sticks (high contrast on white)
show sticks, ligand
color magenta, ligand and elem C
color red, ligand and elem O
color blue, ligand and elem N
color yellow, ligand and elem S
set stick_radius, 0.25, ligand

# R702W as orange spheres (visible, distinct)
show spheres, R702
color orange, R702
set sphere_scale, 0.5, R702

orient complex
zoom complex, 5

# Store scene
scene 001, store

ray 2400, 1800
png C:/Users/vasud/nod2-screening-data/isef_figures/v5_winner/pub_overview_white.png, dpi=300
print "Saved: pub_overview_white.png"

# ============================================================
# FIGURE 2: BINDING SITE WITH POLAR CONTACTS (KEY FIGURE!)
# ============================================================
hide everything
bg_color white

# Protein cartoon - light gray background context
show cartoon, protein
color gray90, protein
set cartoon_transparency, 0.5

# Active site residues as sticks (darker gray carbons)
show sticks, active_site
color gray50, active_site and elem C
util.cnc active_site

# Ligand as MAGENTA (pops against gray)
show sticks, ligand
color deeppurple, ligand and elem C
color red, ligand and elem O
color blue, ligand and elem N
color yellow, ligand and elem S
set stick_radius, 0.22, ligand

# POLAR CONTACTS - This is the key visualization!
# Find all polar contacts from ligand to nearby atoms
distance hbonds, ligand, active_site, mode=2
color yellow, hbonds
hide labels, hbonds

# Key binding residues highlighted
color tv_green, binding_res and elem C

center ligand
zoom ligand, 6

scene 002, store

ray 2400, 1800
png C:/Users/vasud/nod2-screening-data/isef_figures/v5_winner/pub_binding_hbonds.png, dpi=300
print "Saved: pub_binding_hbonds.png"

# ============================================================
# FIGURE 3: BINDING SITE - LABELED VERSION
# ============================================================
# Keep same view, add labels

# Label key residues (one atom per residue)
label binding_res and name CA, resn+resi

# Label ligand
label ligand and name C1, "Febuxostat"

ray 2400, 1800
png C:/Users/vasud/nod2-screening-data/isef_figures/v5_winner/pub_binding_labeled.png, dpi=300
print "Saved: pub_binding_labeled.png"

# Hide labels for next figure
hide labels

# ============================================================
# FIGURE 4: R702W MUTATION CLOSEUP (FIXED - HIGH CONTRAST!)
# ============================================================
hide everything
bg_color white

# Select region around R702 (10A radius)
select R702_region, byres all within 10 of R702

# Show cartoon for context - LIGHT BLUE (visible on white)
show cartoon, R702_region
color lightblue, R702_region

# R702W as BRIGHT ORANGE sticks (HIGH CONTRAST)
show sticks, R702
color tv_orange, R702 and elem C
color red, R702 and elem O
color blue, R702 and elem N
set stick_radius, 0.3, R702

# Nearby residues as thin gray sticks
select nearby_R702, byres all within 5 of R702 and not R702
show sticks, nearby_R702
color gray60, nearby_R702 and elem C
util.cnc nearby_R702
set stick_radius, 0.15, nearby_R702

# Show any polar contacts from R702
distance R702_contacts, R702, nearby_R702, mode=2
color deepteal, R702_contacts
hide labels, R702_contacts

# Label R702
label R702 and name CA, "R702W"

center R702
zoom R702, 8

scene 003, store

ray 2400, 1800
png C:/Users/vasud/nod2-screening-data/isef_figures/v5_winner/pub_r702w_closeup.png, dpi=300
print "Saved: pub_r702w_closeup.png"

# ============================================================
# FIGURE 5: SURFACE + LIGAND IN POCKET
# ============================================================
hide everything
hide labels
bg_color white

# Protein surface - white/gray, semi-transparent
show surface, protein
color white, protein
set transparency, 0.7, protein

# Ligand as spheres INSIDE pocket (visible through surface)
show spheres, ligand
color magenta, ligand and elem C
color red, ligand and elem O
color blue, ligand and elem N
set sphere_scale, 0.4, ligand

# R702W location on surface - orange highlight
color tv_orange, R702

center ligand
zoom ligand, 10

scene 004, store

ray 2400, 1800
png C:/Users/vasud/nod2-screening-data/isef_figures/v5_winner/pub_pocket_surface.png, dpi=300
print "Saved: pub_pocket_surface.png"

# ============================================================
# FIGURE 6: DARK BACKGROUND VERSION (FOR POSTER VARIETY)
# ============================================================
hide everything
bg_color black

# Protein cartoon - slate blue
show cartoon, protein
color slate, protein
set cartoon_transparency, 0.3

# Active site residues - BRIGHT GREEN (visible on black!)
show sticks, active_site
color limegreen, active_site and elem C
color red, active_site and elem O
color tv_blue, active_site and elem N

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

scene 005, store

ray 2400, 1800
png C:/Users/vasud/nod2-screening-data/isef_figures/v5_winner/pub_binding_dark.png, dpi=300
print "Saved: pub_binding_dark.png"

# ============================================================
# FIGURE 7: FULL STRUCTURE WITH DOMAIN COLORS
# ============================================================
hide everything
delete hbonds*
delete R702_contacts
bg_color white

show cartoon, protein

# Color by domain (NOD2 domains)
# CARD1: 1-93 (gray)
# CARD2: 104-191 (gray)
# NBD: 273-577 (marine blue)
# HD1: 578-628 (cyan)
# HD2: 629-743 (orange) - contains R702W!
# LRR: 744-1040 (forest green) - binding site

color gray70, resi 1-191
color marine, resi 273-577
color cyan, resi 578-628
color tv_orange, resi 629-743
color forest, resi 744-1040

# Highlight R702W
show spheres, R702
color red, R702
set sphere_scale, 0.6, R702

# Show ligand location
show spheres, ligand
color magenta, ligand
set sphere_scale, 0.4, ligand

orient protein
zoom protein, 3

scene 006, store

ray 2400, 1800
png C:/Users/vasud/nod2-screening-data/isef_figures/v5_winner/pub_domains_colored.png, dpi=300
print "Saved: pub_domains_colored.png"

# ============================================================
print ""
print "============================================"
print "PUBLICATION QUALITY EXPORTS COMPLETE!"
print ""
print "Files saved:"
print "  pub_overview_white.png     - Full structure (white bg)"
print "  pub_binding_hbonds.png     - Binding site + H-bonds"
print "  pub_binding_labeled.png    - Binding site with labels"
print "  pub_r702w_closeup.png      - R702W mutation detail (FIXED!)"
print "  pub_pocket_surface.png     - Surface with ligand"
print "  pub_binding_dark.png       - Binding site (dark bg)"
print "  pub_domains_colored.png    - Full structure with domains"
print ""
print "KEY TECHNIQUES USED:"
print "  - Polar contacts (H-bonds as yellow dashes)"
print "  - byres within selections"
print "  - High contrast colors"
print "  - Scene storage"
print "============================================"
