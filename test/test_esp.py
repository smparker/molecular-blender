"""ESP mapping in source coordinates on the final displayed density geometry."""
from pathlib import Path
import tempfile
import unittest

import numpy as np
from molecular_blender.esp import sample_cube

try:
    import bpy
except ImportError:
    bpy = None


def write_cube(path, values, atom_x=0, origin=(-2., -2., -2.)):
    axes = ''.join(' '.join(map(str, [-n, *axis])) + '\n'
                   for n, axis in zip(values.shape, np.eye(3) * .5))
    path.write_text('analytic test field\nexplicit scalar units in test\n1 ' + ' '.join(map(str, origin))
                    + '\n' + axes + f'2 0 {atom_x} 0 0\n' + ' '.join(map(str, values.ravel())) + '\n')


class TestInterpolation(unittest.TestCase):
    def test_rotated_left_handed_grid_and_closed_boundary(self):
        axes = np.array([[.2, .2, 0], [-.4, .4, 0], [0, 0, -.8]])
        origin = np.array([2., -3., 5.])
        i, j, k = np.indices((3, 4, 5))
        grid = dict(nres=(3, 4, 5), axes=axes, origin=origin, data=2*i-j+3*k+i*j*k)
        positions = np.array([[.25, 2., 1.5], [2., 3., 4.], [0., 0., 0.]])
        x, y, z = positions.T
        np.testing.assert_allclose(sample_cube(grid, positions @ axes + origin), 2*x-y+3*z+x*y*z)
        with self.assertRaisesRegex(ValueError, 'cover'):
            sample_cube(grid, np.array([[2.01, 1., 2.]]) @ axes + origin)


@unittest.skipUnless(bpy, 'requires Blender')
class TestMappedDensity(unittest.TestCase):
    def test_remeshed_potential_renders_and_restores(self):
        from molecular_blender import scene_panel
        from molecular_blender.plotter import BlendMoleculeSetup, draw_surfaces
        from molecular_blender.scene_state import owned_objects, read_record
        from molecular_blender.scene_surfaces import surface_groups
        from molecular_blender.scene_esp import ATTRIBUTE, field_map_settings
        old = bpy.context.window.scene
        scene = bpy.data.scenes.new('ESP regression')
        bpy.context.window.scene = scene
        objs, colls = set(bpy.data.objects), set(bpy.data.collections)
        cam = bpy.data.cameras.new('Ortho')
        world = bpy.data.worlds.new('Uniform light')
        scene.world = world
        world.use_nodes = True
        world.node_tree.nodes.get('Background').inputs['Color'].default_value = (1, 1, 1, 1)
        scene_panel.register()
        try:
            with tempfile.TemporaryDirectory() as tmp:
                root = Path(tmp)
                x, y, z = .5 * np.indices((9, 9, 9)) - 2
                density, esp = root / 'density.cube', root / 'esp.cube'
                write_cube(density, np.exp(-(x*x+y*y+z*z)))
                write_cube(esp, x)
                model, styler, options = BlendMoleculeSetup(bpy.context, str(density), volume='density',
                    cube_field_type='density', isovalues='.35', cumulative=False, remesh=4)
                draw_surfaces(model, styler, bpy.context, options)
                owner = model.scene_owner
                owner.location.x += 3.
                key, = surface_groups(owner, read_record(owner))
                surf, = owned_objects(owner, 'surface')
                mesh, mat = surf.data, surf.active_material
                for atom in owned_objects(owner, 'atom'):
                    atom.hide_render = True
                def apply():
                    return bpy.ops.molecular_blender.map_esp(group=key, filepath=str(esp), minimum=-.6, maximum=.6)
                self.assertEqual(apply(), {'FINISHED'})
                self.assertFalse(surf.modifiers[0].show_viewport)
                self.assertNotEqual(len(surf.data.vertices), len(mesh.vertices))
                def check_values():
                    values = [v.value for v in surf.data.attributes[ATTRIBUTE].data]
                    np.testing.assert_allclose(values, [v.co.x for v in surf.data.vertices], atol=1e-6)
                check_values()
                write_cube(esp, x, atom_x=.1)
                with self.assertRaisesRegex(RuntimeError, 'atom order and source coordinates'):
                    apply()
                write_cube(esp, x, origin=(2., 2., 2.))
                with self.assertRaisesRegex(RuntimeError, 'cover'):
                    apply()
                check_values()
                write_cube(esp, x)
                self.assertEqual(apply(), {'FINISHED'})
                scene.camera = bpy.data.objects.new('ESP camera', cam)
                scene.collection.objects.link(scene.camera)
                scene.camera.location = (3., 0., 5.)
                cam.type = 'ORTHO'
                cam.ortho_scale = 2.8
                scene.render.engine = 'CYCLES'
                scene.cycles.device = 'CPU'
                scene.cycles.samples = 4
                scene.render.film_transparent = True
                scene.render.resolution_x = scene.render.resolution_y = 48
                scene.render.resolution_percentage = 100
                scene.render.filepath = str(root / 'mapped.png')
                bpy.ops.render.render(write_still=True)
                image = bpy.data.images.load(scene.render.filepath)
                try:
                    pixels = np.asarray(image.pixels[:]).reshape(48, 48, 4)
                    red, blue = pixels[22:26, 11:15].mean((0, 1)), pixels[22:26, 33:37].mean((0, 1))
                    self.assertGreater(red[0], 2 * red[2])
                    self.assertGreater(blue[2], 2 * blue[0])
                finally:
                    bpy.data.images.remove(image)
                # Materials saved by the original ESP-only mapper have no kind or label.
                del surf.active_material['mb_field_kind'], surf.active_material['mb_field_label']
                self.assertEqual(field_map_settings(surf.active_material)['kind'], 'esp')
                for kind, low, high, stops in [('signed', -2., 1., 3), ('unsigned', .25, 2., 2)]:
                    values = x if kind == 'signed' else x+2.
                    write_cube(esp, values)
                    label = 'Density change' if kind == 'signed' else 'Density'
                    args = dict(group=key, filepath=str(esp), field_kind=kind, field_label=label,
                                field_units='e/bohr^3', minimum=low, maximum=high)
                    self.assertEqual(bpy.ops.molecular_blender.map_esp(**args), {'FINISHED'})
                    self.assertNotIn(ATTRIBUTE, surf.data.attributes)
                    np.testing.assert_allclose([p.value for p in surf.data.attributes['mb_field'].data],
                                               [v.co.x + (2. if kind == 'unsigned' else 0.)
                                                for v in surf.data.vertices], atol=1e-6)
                    saved = field_map_settings(surf.active_material)
                    self.assertEqual((saved['label'], saved['units'], saved['kind']), (label, 'e/bohr^3', kind))
                    np.testing.assert_allclose(saved['limits'], [low, high])
                    nodes = surf.active_material.node_tree.nodes
                    self.assertEqual(nodes.get('Attribute').attribute_name, 'mb_field')
                    colors = nodes.get('Color Ramp').color_ramp
                    self.assertEqual(len(colors.elements), stops)
                    white = colors.evaluate(0. if kind == 'unsigned' else -low / (high-low))
                    np.testing.assert_allclose(white[:3], [1., 1., 1.], atol=1e-5)
                    old_data = surf.data
                    with self.assertRaisesRegex(RuntimeError, 'units'):
                        bpy.ops.molecular_blender.map_esp(**dict(args, field_units=''))
                    self.assertEqual(surf.data, old_data)
                write_cube(esp, x)
                self.assertEqual(apply(), {'FINISHED'})
                self.assertNotIn('mb_field', surf.data.attributes)
                check_values()
                self.assertEqual(bpy.ops.molecular_blender.clear_esp(group=key), {'FINISHED'})
                self.assertEqual((surf.data, surf.active_material), (mesh, mat))
                self.assertTrue(surf.modifiers[0].show_render)
        finally:
            scene_panel.unregister()
            bpy.context.window.scene = old
            for obj in set(bpy.data.objects) - objs:
                bpy.data.objects.remove(obj, do_unlink=True)
            for coll in set(bpy.data.collections) - colls:
                bpy.data.collections.remove(coll)
            bpy.data.scenes.remove(scene)
            bpy.data.cameras.remove(cam)
            bpy.data.worlds.remove(world)
