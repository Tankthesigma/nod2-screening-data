"""
BLENDER + MOLECULARNODES SCRIPT
Photorealistic NOD2 protein visualization
ISEF 2026

HOW TO RUN:
1. Open Blender
2. Go to Scripting workspace (top tabs)
3. Click "New" to create new script
4. Paste this entire script
5. Click "Run Script" (play button)

REQUIREMENTS:
- Blender 3.6+ (https://www.blender.org/download/)
- MolecularNodes addon (https://github.com/BradyAJohnston/MolecularNodes/releases)
  Install: Edit > Preferences > Add-ons > Install > select .zip file
"""

import bpy
import os

# ============================================================
# CONFIGURATION
# ============================================================
PDB_PATH = r"C:/Users/vasud/nod2-screening-data/PHASE_6/structures/NOD2_alphafold_full.pdb"
OUT_PATH = r"C:/Users/vasud/nod2-screening-data/isef_figures/v5_winner/ppt_photos/fig1_blender_nod2.png"

# Domain colors (RGBA)
DOMAIN_COLORS = {
    "CARD1": (0.6, 0.6, 0.6, 1.0),       # gray (1-93)
    "CARD2": (0.6, 0.6, 0.6, 1.0),       # gray (104-191)
    "NBD": (0.0, 0.4, 0.8, 1.0),         # blue #0066cc (273-577)
    "HD1": (0.0, 0.8, 0.8, 1.0),         # cyan #00cccc (578-628)
    "HD2": (1.0, 0.55, 0.0, 1.0),        # orange #ff8c00 (629-743)
    "LRR": (0.2, 0.8, 0.2, 1.0),         # green #32cd32 (744-1040)
}

# ============================================================
# SETUP FUNCTIONS
# ============================================================

def clear_scene():
    """Remove all objects from scene"""
    bpy.ops.object.select_all(action='SELECT')
    bpy.ops.object.delete(use_global=False)

    # Clear orphan data
    for block in bpy.data.meshes:
        if block.users == 0:
            bpy.data.meshes.remove(block)

def setup_render_engine():
    """Configure Cycles for high quality"""
    scene = bpy.context.scene
    scene.render.engine = 'CYCLES'

    # Try to use GPU
    prefs = bpy.context.preferences.addons.get('cycles')
    if prefs:
        try:
            bpy.context.preferences.addons['cycles'].preferences.compute_device_type = 'CUDA'
            scene.cycles.device = 'GPU'
        except:
            scene.cycles.device = 'CPU'

    # High quality settings
    scene.cycles.samples = 512  # Increase for final render
    scene.cycles.use_denoising = True
    scene.cycles.use_adaptive_sampling = True

def setup_render_output():
    """Configure 6K output"""
    scene = bpy.context.scene
    scene.render.resolution_x = 6000
    scene.render.resolution_y = 4500
    scene.render.resolution_percentage = 100
    scene.render.image_settings.file_format = 'PNG'
    scene.render.image_settings.color_mode = 'RGBA'
    scene.render.filepath = OUT_PATH

def setup_white_background():
    """Pure white background"""
    world = bpy.context.scene.world
    if world is None:
        world = bpy.data.worlds.new("World")
        bpy.context.scene.world = world

    world.use_nodes = True
    tree = world.node_tree
    tree.nodes.clear()

    bg = tree.nodes.new('ShaderNodeBackground')
    bg.inputs['Color'].default_value = (1, 1, 1, 1)
    bg.inputs['Strength'].default_value = 1.0

    output = tree.nodes.new('ShaderNodeOutputWorld')
    tree.links.new(bg.outputs['Background'], output.inputs['Surface'])

def add_studio_lighting():
    """Professional 3-point lighting setup"""
    # Key light (main)
    bpy.ops.object.light_add(type='AREA', location=(10, -8, 8))
    key = bpy.context.object
    key.name = "Key_Light"
    key.data.energy = 2000
    key.data.size = 6
    key.data.color = (1.0, 0.98, 0.95)

    # Fill light (softer, opposite side)
    bpy.ops.object.light_add(type='AREA', location=(-10, -5, 5))
    fill = bpy.context.object
    fill.name = "Fill_Light"
    fill.data.energy = 800
    fill.data.size = 8

    # Rim light (back)
    bpy.ops.object.light_add(type='AREA', location=(0, 10, 10))
    rim = bpy.context.object
    rim.name = "Rim_Light"
    rim.data.energy = 1200
    rim.data.size = 5

def add_camera():
    """Add camera positioned to view protein"""
    bpy.ops.object.camera_add(location=(0, -25, 8))
    cam = bpy.context.object
    cam.name = "Render_Camera"
    cam.rotation_euler = (1.3, 0, 0)
    bpy.context.scene.camera = cam

    # Depth of field for professional look
    cam.data.dof.use_dof = False  # Set True for blur effect
    cam.data.dof.focus_distance = 25
    cam.data.dof.aperture_fstop = 2.8

def create_material(name, color_rgba, metallic=0.0, roughness=0.3):
    """Create a Principled BSDF material"""
    mat = bpy.data.materials.new(name=name)
    mat.use_nodes = True

    bsdf = mat.node_tree.nodes.get("Principled BSDF")
    if bsdf:
        bsdf.inputs["Base Color"].default_value = color_rgba
        bsdf.inputs["Metallic"].default_value = metallic
        bsdf.inputs["Roughness"].default_value = roughness
        bsdf.inputs["Specular"].default_value = 0.5

    return mat

# ============================================================
# MOLECULARNODES IMPORT
# ============================================================

def import_with_molecular_nodes():
    """Import PDB using MolecularNodes addon"""
    try:
        # MolecularNodes operator (version 4.x)
        bpy.ops.mn.import_pdb(
            filepath=PDB_PATH,
            name="NOD2"
        )
        return bpy.context.object
    except AttributeError:
        print("MolecularNodes not found or different version.")
        print("Please install MolecularNodes addon.")
        print("Download: https://github.com/BradyAJohnston/MolecularNodes/releases")
        return None

# ============================================================
# MAIN EXECUTION
# ============================================================

def main():
    print("=" * 50)
    print("NOD2 BLENDER RENDER SCRIPT")
    print("=" * 50)

    # Setup scene
    clear_scene()
    setup_render_engine()
    setup_render_output()
    setup_white_background()
    add_studio_lighting()
    add_camera()

    print(f"Loading: {PDB_PATH}")

    # Check if file exists
    if not os.path.exists(PDB_PATH):
        print(f"ERROR: PDB file not found: {PDB_PATH}")
        return

    # Import protein
    protein = import_with_molecular_nodes()

    if protein is None:
        print("Failed to import protein. Check MolecularNodes installation.")
        return

    print("Protein imported successfully!")
    print(f"Output will be saved to: {OUT_PATH}")
    print("")
    print("NEXT STEPS:")
    print("1. In MolecularNodes panel, select 'Cartoon' style")
    print("2. Use 'Color' node to color by residue range")
    print("3. Add 'Sphere' style for R702 (residue 702) - RED")
    print("4. Add 'Sphere' style for pocket (1007,1008,1011,1034,1036,1037) - MAGENTA")
    print("5. Press F12 to render")
    print("")
    print("Or manually render:")
    print("  bpy.ops.render.render(write_still=True)")

# Run
if __name__ == "__main__":
    main()
