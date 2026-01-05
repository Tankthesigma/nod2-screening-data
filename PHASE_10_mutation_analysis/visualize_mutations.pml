# PyMOL Visualization Script for NOD2 Mutation Analysis
# Run in PyMOL: @visualize_mutations.pml

# Load structures
load C:/Users/vasud/nod2-screening-data/cross_validation/human_NOD2_LRR.pdb, NOD2_LRR
load C:/Users/vasud/nod2-screening-data/PHASE_6/structures/febuxostat_docked.sdf, febuxostat
load C:/Users/vasud/nod2-screening-data/PHASE_6/structures/natural_top_docked.sdf, natural_top

# Set up display
bg_color white
hide all
show cartoon, NOD2_LRR
color gray80, NOD2_LRR

# Show ligands
show sticks, febuxostat
show sticks, natural_top
color cyan, febuxostat
color magenta, natural_top

# Highlight binding pocket residues (1007-1011, 1034-1037)
select pocket, resi 1007-1011+1034-1037 and NOD2_LRR
show sticks, pocket
color yellow, pocket

# Highlight mutation sites (using LRR numbering if in structure)
# R702W - not in LRR (position 702 < 744)
# G908R - in LRR
select mut_G908R, resi 908 and NOD2_LRR
show spheres, mut_G908R and name CA
color red, mut_G908R

# N852S - in LRR
select mut_N852S, resi 852 and NOD2_LRR
show spheres, mut_N852S and name CA
color orange, mut_N852S

# M863V - in LRR
select mut_M863V, resi 863 and NOD2_LRR
show spheres, mut_M863V and name CA
color salmon, mut_M863V

# L1007fs truncation point
select mut_L1007fs, resi 1006-1007 and NOD2_LRR
show spheres, mut_L1007fs and name CA
color firebrick, mut_L1007fs

# Label mutation sites
label mut_G908R and name CA, "G908R"
label mut_N852S and name CA, "N852S"
label mut_M863V and name CA, "M863V"
label mut_L1007fs and name CA and resi 1007, "L1007fs"

# Set view
set_view (\
     0.9,    0.0,    0.4,\
     0.0,    1.0,    0.0,\
    -0.4,    0.0,    0.9,\
     0.0,    0.0, -100.0,\
     0.0,    0.0,    0.0,\
   50.0,  150.0,  -20.0 )

# Ray trace and save
ray 2400, 1800
png C:/Users/vasud/nod2-screening-data/PHASE_10_mutation_analysis/STRATIFICATION_MAP.png, dpi=300

# Also save session
save C:/Users/vasud/nod2-screening-data/PHASE_10_mutation_analysis/mutation_visualization.pse
