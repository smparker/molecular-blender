"""Render a water density surface colored by the electrostatic potential from the same RHF calculation.

blender -b --factory-startup --python-exit-code 1 --python examples/render_water_esp.py -- OUTPUT_DIR
"""
from pathlib import Path
import sys

import addon_utils
import bpy
from mathutils import Matrix

examples = Path(__file__).resolve().parent
sys.path.insert(0, str(examples.parent))
addon_utils.enable('molecular_blender')
from molecular_blender.plotter import BlendMoleculeSetup, draw_surfaces
from molecular_blender.scene_state import read_record

output = Path(sys.argv[sys.argv.index('--') + 1]).resolve()
output.mkdir(parents=True, exist_ok=True)
inputs = examples / 'water_esp'
scene = bpy.data.scenes.new('Water electrostatic potential')
bpy.context.window.scene = scene
model, styler, settings = BlendMoleculeSetup(
    bpy.context, str(inputs / 'water.molden'), volume='density', cumulative=False,
    isovalues='0.001', resolution=.055, remesh=0, ignore_hydrogen=False)
draw_surfaces(model, styler, bpy.context, settings)
owner = model.scene_owner
bpy.context.view_layer.objects.active = owner
group, = {surface['group'] for surface in read_record(owner)['surfaces']}
actions = bpy.ops.molecular_blender
actions.map_esp(group=group, filepath=str(inputs / 'potential.cube'),
                units='Hartree/e', minimum=-.035, maximum=.035)

# Look along the molecular plane normal with both hydrogens below the oxygen.
oxygen, h1, h2 = (atom.position for atom in model.atoms)
right = (h2 - h1).normalized()
normal = right.cross(oxygen - (h1 + h2) / 2).normalized()
up = normal.cross(right)
center = (oxygen + h1 + h2) / 3
camera = bpy.data.cameras.new('Water figure camera')
camera.type = 'ORTHO'
camera.ortho_scale = 4.7
scene.camera = bpy.data.objects.new(camera.name, camera)
scene.collection.objects.link(scene.camera)
scene.camera.rotation_euler = Matrix((right, up, normal)).transposed().to_euler()
scene.camera.location = center + 7 * normal

scene.world = bpy.data.worlds.new('Figure fill')
scene.world.use_nodes = True
scene.world.node_tree.nodes['Background'].inputs['Color'].default_value = (.65, .65, .65, 1)
lamp = bpy.data.lights.new('Softbox', 'AREA')
lamp.energy, lamp.size = 200, 5
light = bpy.data.objects.new(lamp.name, lamp)
scene.collection.objects.link(light)
light.location = center + 4 * normal + 3 * up - 2 * right
light.rotation_euler = (center - light.location).to_track_quat('-Z', 'Y').to_euler()
scene.render.engine = 'CYCLES'
scene.cycles.samples = 32
scene.view_settings.view_transform = 'Standard'
scene.render.resolution_x = scene.render.resolution_y = 480
scene.render.resolution_percentage = 100
actions.export_surface_series(directory=str(output), legends=True)
