from pymol import cmd

# Create Gaussian density map around mutation
cmd.do("map_new mut_map, gaussian, 0.6, (mutation), 6.0")

# Create outer halo
cmd.do("isosurface mut_halo, mut_map, 0.5")
cmd.do("color magenta, mut_halo")
cmd.do("set surface_color, magenta, mut_halo")
cmd.do("set transparency, 0.6, mut_halo")

print("Glow created!")
