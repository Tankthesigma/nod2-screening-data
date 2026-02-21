# ============================================================
# ISEF 2026 - COMPLETE POSTER FIGURES
# Uses BOTH structures:
#   1. NOD2_alphafold_full.pdb (full protein with R702)
#   2. complex_febuxostat.pdb (LRR + Febuxostat ligand)
# ============================================================

reinitialize

# ============================================================
# PART A: FULL NOD2 STRUCTURE WITH DOMAIN COLORING
# Shows where R702W mutation is relative to binding site
# ============================================================

load C:/Users/vasud/nod2-screening-data/PHASE_6/structures/NOD2_alphafold_full.pdb, nod2_full

# Create selections for domains
select CARD1, nod2_full and resi 1-93
select CARD2, nod2_full and resi 104-191
select NBD, nod2_full and resi 273-577
select HD1, nod2_full and resi 578-628
select HD2, nod2_full and resi 629-743
select LRR, nod2_full and resi 744-1040

# R702 mutation site (in HD2 domain!)
select R702_site, nod2_full and resi 702

# Binding pocket location (in LRR domain)
select binding_pocket, nod2_full and resi 1007+1008+1011+1034+1036+1037

# ============================================================
# GLOBAL SETTINGS
# ============================================================
set antialias, 2
set ray_trace_mode, 1
set ray_shadows, 1
set ray_opaque_background, 1
set cartoon_fancy_helices, 1
set cartoon_smooth_loops, 1

# ============================================================
# FIGURE A1: FULL NOD2 WITH DOMAIN COLORS (WHITE BG)
# ============================================================
bg_color white
hide everything

show cartoon, nod2_full

# Color by domain
color gray70, CARD1
color gray70, CARD2
color marine, NBD
color cyan, HD1
color tv_orange, HD2
color tv_green, LRR

# Highlight R702 as RED spheres
show spheres, R702_site
color red, R702_site
set sphere_scale, 1.0, R702_site

# Highlight binding pocket as MAGENTA spheres
show spheres, binding_pocket
color magenta, binding_pocket
set sphere_scale, 0.6, binding_pocket

orient nod2_full
zoom nod2_full, 3

scene full_domains, store

ray 2400, 1800
png C:/Users/vasud/nod2-screening-data/isef_figures/v5_winner/figA1_full_domains_white.png, dpi=300
print "Saved: figA1_full_domains_white.png"

# ============================================================
# FIGURE A2: FULL NOD2 (DARK BG - POSTER STYLE)
# ============================================================
bg_color black

ray 2400, 1800
png C:/Users/vasud/nod2-screening-data/isef_figures/v5_winner/figA2_full_domains_dark.png, dpi=300
print "Saved: figA2_full_domains_dark.png"

# ============================================================
# FIGURE A3: R702W CLOSEUP (NOW WITH ACTUAL RESIDUE!)
# ============================================================
hide everything
bg_color white

# Region around R702
select R702_region, nod2_full and byres (resi 702 around 12)

show cartoon, R702_region
color tv_orange, R702_region

# R702 as sticks - BRIGHT RED
show sticks, R702_site
color red, R702_site and elem C
color red, R702_site and elem O
color blue, R702_site and elem N
set stick_radius, 0.3, R702_site

# Nearby residues
select nearby_702, byres (R702_site around 5) and not R702_site
show sticks, nearby_702
color gray60, nearby_702 and elem C
util.cnc nearby_702
set stick_radius, 0.15, nearby_702

center R702_site
zoom R702_site, 8

scene r702_closeup, store

ray 2400, 1800
png C:/Users/vasud/nod2-screening-data/isef_figures/v5_winner/figA3_r702w_closeup.png, dpi=300
print "Saved: figA3_r702w_closeup.png"

# ============================================================
# PART B: BINDING SITE WITH FEBUXOSTAT
# ============================================================

# Load the complex (LRR + ligand)
load C:/Users/vasud/nod2-screening-data/MAC_ALL_FILES/md_simulation/structures/complex_febuxostat.pdb, complex

# Selections for complex
select protein_lrr, complex and chain A
select febuxostat, complex and resn FEB
select active_site, (byres chain A within 5 of febuxostat) and not febuxostat
select key_residues, complex and chain A and resi 1007+1008+1011+1034+1036+1037

# ============================================================
# FIGURE B1: BINDING SITE WITH H-BONDS (WHITE BG)
# ============================================================
hide everything, nod2_full
bg_color white

show cartoon, protein_lrr
color gray80, protein_lrr
set cartoon_transparency, 0.5

# Active site sticks
show sticks, active_site
color gray50, active_site and elem C
util.cnc active_site
set stick_radius, 0.18, active_site

# Key residues in GREEN
color tv_green, key_residues and elem C

# Febuxostat - ORANGE (pops!)
show sticks, febuxostat
color tv_orange, febuxostat and elem C
color red, febuxostat and elem O
color blue, febuxostat and elem N
color yellow, febuxostat and elem S
set stick_radius, 0.25, febuxostat

# H-bonds
distance hbonds, febuxostat, active_site, mode=2
color yellow, hbonds
set dash_width, 2.5
hide labels, hbonds

center febuxostat
zoom febuxostat, 6

scene binding_hbonds, store

ray 2400, 1800
png C:/Users/vasud/nod2-screening-data/isef_figures/v5_winner/figB1_binding_hbonds_white.png, dpi=300
print "Saved: figB1_binding_hbonds_white.png"

# ============================================================
# FIGURE B2: BINDING SITE (DARK BG)
# ============================================================
bg_color black

# Adjust colors for dark bg
color marine, protein_lrr
color tv_green, active_site and elem C
color tv_orange, febuxostat and elem C

ray 2400, 1800
png C:/Users/vasud/nod2-screening-data/isef_figures/v5_winner/figB2_binding_hbonds_dark.png, dpi=300
print "Saved: figB2_binding_hbonds_dark.png"

# ============================================================
# FIGURE B3: CLEAN BINDING (NO H-BONDS - FOR POWERPOINT LABELS)
# ============================================================
hide everything
delete hbonds
bg_color white

show cartoon, protein_lrr
color gray85, protein_lrr
set cartoon_transparency, 0.4

show sticks, key_residues
color tv_green, key_residues and elem C
util.cnc key_residues
set stick_radius, 0.2, key_residues

show sticks, febuxostat
color tv_orange, febuxostat and elem C
color red, febuxostat and elem O
color blue, febuxostat and elem N
set stick_radius, 0.25, febuxostat

center febuxostat
zoom febuxostat, 5

scene binding_clean, store

ray 2400, 1800
png C:/Users/vasud/nod2-screening-data/isef_figures/v5_winner/figB3_binding_clean.png, dpi=300
print "Saved: figB3_binding_clean.png"

# ============================================================
# FIGURE B4: SURFACE + LIGAND IN POCKET
# ============================================================
hide everything
bg_color white

show surface, protein_lrr
color white, protein_lrr
set transparency, 0.7, protein_lrr

show spheres, febuxostat
color tv_orange, febuxostat and elem C
color red, febuxostat and elem O
set sphere_scale, 0.45, febuxostat

center febuxostat
zoom febuxostat, 10

scene pocket_surface, store

ray 2400, 1800
png C:/Users/vasud/nod2-screening-data/isef_figures/v5_winner/figB4_pocket_surface.png, dpi=300
print "Saved: figB4_pocket_surface.png"

# ============================================================
# FIGURE B5: KEY RESIDUES ONLY (PLIP STYLE)
# ============================================================
hide everything
bg_color white

show sticks, key_residues
color tv_green, key_residues and elem C
util.cnc key_residues
set stick_radius, 0.22, key_residues

show sticks, febuxostat
color tv_orange, febuxostat and elem C
color red, febuxostat and elem O
color blue, febuxostat and elem N
set stick_radius, 0.28, febuxostat

# Show contacts
distance key_hbonds, febuxostat, key_residues, mode=2
color yellow, key_hbonds
set dash_width, 3
hide labels, key_hbonds

center febuxostat
zoom febuxostat, 4

scene key_residues_view, store

ray 2400, 1800
png C:/Users/vasud/nod2-screening-data/isef_figures/v5_winner/figB5_key_residues.png, dpi=300
print "Saved: figB5_key_residues.png"

# ============================================================
print ""
print "============================================"
print "COMPLETE POSTER FIGURES DONE!"
print ""
print "PART A - FULL NOD2 STRUCTURE:"
print "  figA1_full_domains_white.png - Domain colors (white bg)"
print "  figA2_full_domains_dark.png  - Domain colors (dark bg)"
print "  figA3_r702w_closeup.png      - R702W mutation site"
print ""
print "PART B - FEBUXOSTAT BINDING:"
print "  figB1_binding_hbonds_white.png - Binding + H-bonds (white)"
print "  figB2_binding_hbonds_dark.png  - Binding + H-bonds (dark)"
print "  figB3_binding_clean.png        - Clean (for PPT labels)"
print "  figB4_pocket_surface.png       - Surface view"
print "  figB5_key_residues.png         - Key PLIP residues"
print ""
print "DOMAIN COLOR KEY:"
print "  GRAY   = CARD1/CARD2 (N-terminal)"
print "  MARINE = NBD (nucleotide binding)"
print "  CYAN   = HD1"
print "  ORANGE = HD2 (contains R702W!)"
print "  GREEN  = LRR (ligand binding)"
print "  RED    = R702W mutation site"
print "  MAGENTA= Binding pocket residues"
print "============================================"
