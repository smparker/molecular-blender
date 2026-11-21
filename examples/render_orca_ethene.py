"""Render the occupied ethene π orbital from the checked-in ORCA calculation.

blender -b --factory-startup --python-exit-code 1 --python examples/render_orca_ethene.py -- OUTPUT_DIR [ethene.molden.input]
"""
from pathlib import Path
import sys

import addon_utils
import bpy
from mathutils import Matrix

repo = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(repo))
addon_utils.enable('molecular_blender')
from molecular_blender.plotter import BlendMoleculeSetup, draw_surfaces
from molecular_blender.scene_state import ID, owned_objects, read_record

args = sys.argv[sys.argv.index('--') + 1:]
output = Path(args[0]).resolve()
source = Path(args[1]).resolve() if len(args) > 1 else repo / 'test/fixtures/orca_ethene/ethene.molden'
output.mkdir(parents=True, exist_ok=True)
scene = bpy.data.scenes.new('Ethene π orbital')
bpy.context.window.scene = scene
# ORCA counts from zero. MO("pi.cube",7,0) selects the eighth Molden orbital.
model, styler, settings = BlendMoleculeSetup(
    bpy.context, str(source), orbital='8', molden_convention='auto', isovalues='.11',
    cumulative=False, resolution=.055, remesh=0, plot_style='bs', ignore_hydrogen=False)
if read_record(model.scene_owner)['molden_convention'] != 'orca':
    raise ValueError('Expected an ORCA producer title in the ethene Molden file')
draw_surfaces(model, styler, bpy.context, settings)
# Use solid, smoothly shaded lobes for the side view.
entries = {s['id']: s for s in read_record(model.scene_owner)['surfaces']}
for obj in owned_objects(model.scene_owner, 'surface'):
    sign = 'plus' if entries[obj[ID]]['isovalue'] > 0 else 'minus'
    obj.active_material = styler.isosurface_material('Ethene π ' + sign)
    for polygon in obj.data.polygons:
        polygon.use_smooth = True

# Keep C=C horizontal. Tilt the view slightly out of the molecular plane to show
# the hydrogens as well as the orbital lobes above and below the carbons.
c1, c2, h1 = (atom.position for atom in model.atoms[:3])
center = (c1 + c2) / 2
right = (c2 - c1).normalized()
normal = right.cross(h1 - c1).normalized()
view = (right.cross(normal) + .35 * normal).normalized()
up = view.cross(right)
camera = bpy.data.cameras.new('Ethene camera')
camera.type, camera.ortho_scale = 'ORTHO', 5.1
scene.camera = bpy.data.objects.new(camera.name, camera)
scene.collection.objects.link(scene.camera)
scene.camera.location = center + 8 * view
scene.camera.rotation_euler = Matrix((right, up, view)).transposed().to_euler()

scene.world = bpy.data.worlds.new('Soft fill')
scene.world.use_nodes = True
scene.world.node_tree.nodes['Background'].inputs['Color'].default_value = (.65, .65, .65, 1)
lamp = bpy.data.lights.new('Key light', 'AREA')
lamp.energy, lamp.size = 240, 5
light = bpy.data.objects.new(lamp.name, lamp)
scene.collection.objects.link(light)
light.location = center + 4 * view + 3 * up - 2 * right
light.rotation_euler = (center - light.location).to_track_quat('-Z', 'Y').to_euler()
scene.render.engine = 'CYCLES'
scene.cycles.samples = 24
scene.render.resolution_x, scene.render.resolution_y = 480, 360
scene.render.resolution_percentage = 100
bpy.context.view_layer.objects.active = model.scene_owner
bpy.ops.molecular_blender.export_surface_series(directory=str(output))
