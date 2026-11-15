"""Import bond_orders.sdf with its explicit orders in two rendering styles.

Run from any directory with Blender --python /path/to/examples/sdf_bonds.py.
Optional arguments after --: --save PATH and --render PATH. Clears the scene.
"""
import argparse
from pathlib import Path
import sys

import bpy
from mathutils import Vector

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from molecular_blender.plotter import BlendMoleculeSetup

parser = argparse.ArgumentParser()
parser.add_argument('--save')
parser.add_argument('--render')
args = parser.parse_args(sys.argv[sys.argv.index('--') + 1:] if '--' in sys.argv else [])

bpy.ops.object.select_all(action='SELECT')
bpy.ops.object.delete(use_global=False)
source = Path(__file__).with_name('bond_orders.sdf')
for style, y, label in [('sticks', 2, 'SDF import: sticks'),
                         ('fixedbs', -2.2, 'SDF import: ball and stick')]:
    molecule, _, _ = BlendMoleculeSetup(
        bpy.context, str(source), plot_style=style, bond_thickness=0.14,
        ignore_hydrogen=False, hook_atoms='on')
    bpy.data.objects[molecule.name].location += Vector((0.7, y - 0.6, 0))
    font = bpy.data.curves.new(label, 'FONT')
    font.body = label
    font.align_x = 'CENTER'
    font.size = 0.3
    obj = bpy.data.objects.new(label, font)
    bpy.context.collection.objects.link(obj)
    obj.location = (0, y - 2.1, 0)

scene = bpy.context.scene
scene.render.engine = 'CYCLES'
scene.cycles.samples = 48
scene.cycles.use_denoising = True
scene.world.color = (0.3, 0.3, 0.3)
bpy.ops.object.camera_add(location=(0, -0.4, 18))
camera = bpy.context.object
camera.rotation_euler = (Vector((0, -0.4, 0)) - camera.location).to_track_quat('-Z', 'Y').to_euler()
camera.data.type = 'ORTHO'
camera.data.ortho_scale = 10
scene.camera = camera
bpy.ops.object.light_add(type='AREA', location=(1, 1, 8))
bpy.context.object.data.energy = 1800
bpy.context.object.data.size = 8
scene.render.resolution_x = 1200
scene.render.resolution_y = 1100
scene.render.resolution_percentage = 100
if args.save:
    bpy.ops.wm.save_as_mainfile(filepath=str(Path(args.save).resolve()))
if args.render:
    scene.render.filepath = str(Path(args.render).resolve())
    bpy.ops.render.render(write_still=True)
