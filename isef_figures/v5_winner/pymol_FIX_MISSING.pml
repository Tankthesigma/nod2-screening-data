# Fix script for missing figures B4 and B5
# Run this AFTER the main script (selections still loaded)

# ============================================================
# FIGURE B4: SURFACE + LIGAND (PROPER WAY)
# ============================================================
hide everything
bg_color white

# Surface - like the old working script
show surface, protein_lrr
color white, protein_lrr
set transparency, 0.7, protein_lrr
set surface_quality, 1

# Ligand as spheres inside
show spheres, febuxostat
color tv_orange, febuxostat and elem C
color red, febuxostat and elem O
color blue, febuxostat and elem N
set sphere_scale, 0.45, febuxostat

center febuxostat
zoom febuxostat, 10

ray 2400, 1800
png C:/Users/vasud/nod2-screening-data/isef_figures/v5_winner/figB4_pocket_surface.png, dpi=300
print "Saved: figB4_pocket_surface.png"

# ============================================================
# FIGURE B5: KEY RESIDUES ONLY
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

# H-bonds
distance key_hbonds, febuxostat, key_residues, mode=2
color yellow, key_hbonds
set dash_width, 3
hide labels, key_hbonds

center febuxostat
zoom febuxostat, 4

ray 2400, 1800
png C:/Users/vasud/nod2-screening-data/isef_figures/v5_winner/figB5_key_residues.png, dpi=300
print "Saved: figB5_key_residues.png"

print ""
print "DONE! Missing figures created."
