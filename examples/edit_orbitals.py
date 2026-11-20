"""Style a water molecule and export two orbitals with the startup scene camera.

blender --background --factory-startup --python examples/edit_orbitals.py -- OUTPUT_DIR
"""
from pathlib import Path
import sys

import addon_utils
import bpy

examples = Path(__file__).resolve().parent
sys.path.insert(0, str(examples.parent))
addon_utils.enable('molecular_blender')
from molecular_blender.plotter import BlendMoleculeSetup, draw_surfaces

output = Path(sys.argv[sys.argv.index('--') + 1]).resolve()
output.mkdir(parents=True, exist_ok=True)
scene = bpy.context.scene
if scene.camera is None:
    raise ValueError('This example uses the active scene camera')
if '--factory-startup' in sys.argv and scene.objects.get('Cube'):
    scene.objects['Cube'].hide_render = True

model, styler, settings = BlendMoleculeSetup(
    bpy.context, str(examples / 'water-sto-3g.molden'), orbital='homo,lumo',
    resolution=.12, isovalues='.5', remesh=0, ignore_hydrogen=False, plot_style='fixedbs')
draw_surfaces(model, styler, bpy.context, settings)
model.scene_owner.location = (0, 0, 0)
actions = bpy.ops.molecular_blender
actions.select_atom_indices(indices='1,3')
actions.style_selected(color=(.15, .5, .8))

scene.render.engine = 'CYCLES'
scene.cycles.samples = 8
scene.render.resolution_x = scene.render.resolution_y = 480
scene.render.resolution_percentage = 100
actions.export_surface_series(directory=str(output))
