"""Manual multiple bonds; run from the repo root with Blender's --python option.

No bond orders are inferred. Every strand is extruded from a two-point curve.
The optional --render PATH argument writes a preview; --save PATH saves a .blend.
"""
import argparse
from pathlib import Path
import sys

import bpy
from mathutils import Vector

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from molecular_blender.molecule import Atom, Bond, Molecule
from molecular_blender.plotter import PlotMolecule, default_options, make_materials, plot_prep
from molecular_blender.stylers import get_styler


def draw(name, symbols, positions, bonds, xoffset, y, style):
    atoms = [Atom(s, (x + xoffset, py + y, z), i)
             for i, (s, (x, py, z)) in enumerate(zip(symbols, positions))]
    mol = Molecule(name, atoms)
    mol.bonds = [Bond(atoms[i], atoms[j], order=order, direction=(0, 1, 0))
                 for i, j, order in bonds]
    options = default_options({'plot_style': style, 'bond_thickness': 0.14, 'draw_bond_orders': True,
                               'hook_atoms': 'on'})
    make_materials(mol, get_styler(options), options)
    plot_prep(mol, options)
    PlotMolecule(bpy.context, mol, options)
    label = bpy.data.curves.new(name + '_label', 'FONT')
    label.body = name.split(' / ')[0]
    label.align_x = 'CENTER'
    label.size = 0.25
    obj = bpy.data.objects.new(name + '_label', label)
    bpy.context.collection.objects.link(obj)
    obj.location = (xoffset, y - 1.4, 0)


# Use a fresh scene when running this example.
bpy.ops.object.select_all(action='SELECT')
bpy.ops.object.delete(use_global=False)
examples = [
    ('Ethane: single bond', ['c', 'c', 'h', 'h', 'h', 'h', 'h', 'h'],
     [(-0.77, 0, 0), (0.77, 0, 0), (-1.13, 1.02, 0), (-1.13, -0.51, 0.88),
      (-1.13, -0.51, -0.88), (1.13, -1.02, 0), (1.13, 0.51, 0.88), (1.13, 0.51, -0.88)],
     [(0, 1, 1), (0, 2, 1), (0, 3, 1), (0, 4, 1), (1, 5, 1), (1, 6, 1), (1, 7, 1)]),
    ('Ethene: double bond', ['c', 'c', 'h', 'h', 'h', 'h'],
     [(-0.67, 0, 0), (0.67, 0, 0), (-1.23, 0.93, 0), (-1.23, -0.93, 0),
      (1.23, 0.93, 0), (1.23, -0.93, 0)],
     [(0, 1, 2), (0, 2, 1), (0, 3, 1), (1, 4, 1), (1, 5, 1)]),
    ('Formaldehyde: double bond', ['c', 'o', 'h', 'h'],
     [(-0.6, 0, 0), (0.6, 0, 0), (-1.15, 0.94, 0), (-1.15, -0.94, 0)],
     [(0, 1, 2), (0, 2, 1), (0, 3, 1)]),
    ('Acetylene: triple bond', ['c', 'c', 'h', 'h'],
     [(-0.6, 0, 0), (0.6, 0, 0), (-1.66, 0, 0), (1.66, 0, 0)],
     [(0, 1, 3), (0, 2, 1), (1, 3, 1)]),
]
for xoffset, style, heading in [(-3.1, 'sticks', 'Sticks'),
                                 (3.1, 'fixedbs', 'Ball and stick')]:
    for row, (name, symbols, positions, bonds) in enumerate(examples):
        draw(name + ' / ' + style, symbols, positions, bonds,
             xoffset, 4.8 - row * 3.3, style)
    text = bpy.data.curves.new(heading, 'FONT')
    text.body = heading
    text.align_x = 'CENTER'
    text.size = 0.4
    obj = bpy.data.objects.new(heading, text)
    bpy.context.collection.objects.link(obj)
    obj.location = (xoffset, 6.5, 0)

scene = bpy.context.scene
scene.render.engine = 'CYCLES'
scene.cycles.samples = 48
scene.cycles.use_denoising = True
scene.world.color = (0.3, 0.3, 0.3)
bpy.ops.object.camera_add(location=(0, -0.3, 20))
camera = bpy.context.object
camera.rotation_euler = (Vector((0, -0.3, 0)) - camera.location).to_track_quat('-Z', 'Y').to_euler()
camera.data.type = 'ORTHO'
camera.data.ortho_scale = 16
scene.camera = camera
bpy.ops.object.light_add(type='AREA', location=(1, 1, 8))
bpy.context.object.data.energy = 1800
bpy.context.object.data.shape = 'DISK'
bpy.context.object.data.size = 8
scene.render.resolution_x = 1400
scene.render.resolution_y = 1600
scene.render.resolution_percentage = 100
args = argparse.ArgumentParser()
args.add_argument('--render')
args.add_argument('--save')
opts = args.parse_args(sys.argv[sys.argv.index('--') + 1:] if '--' in sys.argv else [])
if opts.save:
    bpy.ops.wm.save_as_mainfile(filepath=str(Path(opts.save).resolve()))
if opts.render:
    scene.render.filepath = str(Path(opts.render).resolve())
    bpy.ops.render.render(write_still=True)
