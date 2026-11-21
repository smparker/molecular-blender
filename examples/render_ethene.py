"""Render an ethene HOMO (TURBOMOLE or ORCA) or valence density (xTB).

blender -b --factory-startup --python-exit-code 1 --python examples/render_ethene.py -- WRITER OUTPUT_DIR
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


def render(writer, output, source=None):
    fields = repo / 'examples/ethene_fields'
    inputs = {
        'orca': (repo / 'test/fixtures/orca_ethene/ethene.molden', dict(orbital='8', molden_convention='auto')),
        'turbomole': (fields / 'turbomole.molden', dict(orbital='homo', molden_convention='turbomole_legacy')),
        'xtb': (fields / 'xtb-density.cub', dict(cube_field_type='density', cube_both_signs=False)),
    }
    default_path, options = inputs[writer]
    output = Path(output).resolve()
    output.mkdir(parents=True, exist_ok=True)
    scene = bpy.data.scenes.new(f'{writer} ethene')
    bpy.context.window.scene = scene
    density = writer == 'xtb'
    model, styler, settings = BlendMoleculeSetup(
        bpy.context, str(source or default_path), **options, isovalues='.035' if density else '.11',
        cumulative=False, resolution=.055, remesh=0, plot_style='bs', ignore_hydrogen=False)
    if writer == 'orca' and read_record(model.scene_owner)['molden_convention'] != 'orca':
        raise ValueError('Expected an ORCA producer title in the ethene Molden file')
    draw_surfaces(model, styler, bpy.context, settings)
    entries = {s['id']: s for s in read_record(model.scene_owner)['surfaces']}
    lobes = owned_objects(model.scene_owner, 'surface')
    if len(lobes) != (1 if density else 2):
        raise ValueError(f'Unexpected number of surfaces for {writer} ethene')
    for obj in lobes:
        sign = 'plus' if entries[obj[ID]]['isovalue'] > 0 else 'minus'
        obj.active_material = (styler.outer_isosurface_material('Valence density plus') if density
                               else styler.isosurface_material('Ethene π ' + sign))
        for polygon in obj.data.polygons:
            polygon.use_smooth = True

    # Keep C=C horizontal. Tilt the view slightly out of the molecular plane to show
    # the hydrogens as well as the orbital lobes above and below the carbons.
    c1, c2, h1 = (atom.position for atom in model.atoms[:3])
    center = (c1 + c2) / 2
    right = (c2 - c1).normalized()
    normal = right.cross(h1 - c1).normalized()
    view = normal if density else (right.cross(normal) + .35 * normal).normalized()
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
    if bpy.ops.molecular_blender.export_surface_series(directory=str(output)) != {'FINISHED'}:
        raise RuntimeError('Ethene figure export failed')


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('writer', choices=('turbomole', 'xtb', 'orca'))
    parser.add_argument('output_dir', type=Path)
    args = parser.parse_args(sys.argv[sys.argv.index('--') + 1:] if '--' in sys.argv else [])
    render(args.writer, args.output_dir)
