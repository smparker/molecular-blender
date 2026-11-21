"""ESP mapping in source coordinates on the final displayed density geometry."""
import json
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
    def test_replacing_field_maps_keeps_compatible_ramps_per_contour(self):
        from molecular_blender import scene_panel
        from molecular_blender.field_scale import scale_values
        from molecular_blender.plotter import BlendMoleculeSetup, draw_surfaces
        from molecular_blender.scene_esp import scale_nodes
        from molecular_blender.scene_state import read_record
        from molecular_blender.scene_surfaces import surface_groups
        scene_panel.register()
        self.addCleanup(scene_panel.unregister)
        blocks = set(bpy.data.objects) | set(bpy.data.collections)
        self.addCleanup(lambda: bpy.data.batch_remove((set(bpy.data.objects) | set(bpy.data.collections)) - blocks))
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            density, first, second = (directory / name for name in ('density.cube', 'x.cube', 'y.cube'))
            x, y, z = .5 * np.indices((9, 9, 9)) - 2
            write_cube(density, np.exp(-(x*x + y*y + z*z)))
            write_cube(first, x)
            write_cube(second, y)
            model, styler, options = BlendMoleculeSetup(bpy.context, str(density), volume='density',
                isovalues='.35,.6', cumulative=False, remesh=0)
            draw_surfaces(model, styler, bpy.context, options)
            owner = model.scene_owner
            bpy.context.view_layer.objects.active = owner
            key, group = next(iter(surface_groups(owner, read_record(owner)).items()))
            members = group['objects']
            a, b = members
            op = bpy.ops.molecular_blender.map_esp
            args = dict(group=key, field_kind='signed', field_units='e/bohr³', minimum=-2., maximum=2.)
            self.assertEqual(op(filepath=str(first), palette='magma', **args), {'FINISHED'})
            b.active_material = a.active_material.copy()
            ramp = scale_nodes(b.active_material)[0].color_ramp
            ramp.elements[5].color = (.2, .7, .3, 1.)
            ramp.elements[5].position = .17
            ramp.interpolation = 'EASE'
            edited = [scale_values(obj.active_material) for obj in members]
            self.assertNotEqual(*edited)
            self.assertEqual(op(filepath=str(second), **args), {'FINISHED'})
            for obj, expected in zip(members, edited):
                self.assertEqual(scale_values(obj.active_material), expected)
                np.testing.assert_allclose([v.value for v in obj.data.attributes['mb_field'].data],
                                           [v.co.y for v in obj.data.vertices], atol=1e-6)
            manifest = directory / 'sum.json'
            entry = dict(units='e/bohr³', linear_combination=[
                dict(filename=first.name, coefficient=1.), dict(filename=second.name, coefficient=-1.)])
            manifest.write_text(json.dumps(dict(molecules=[entry])))
            self.assertEqual(op(filepath=str(manifest), **dict(args, field_units='ignored')), {'FINISHED'})
            self.assertEqual([scale_values(o.active_material) for o in members], edited)
            self.assertEqual(op(filepath=str(second), palette='viridis', **args), {'FINISHED'})
            self.assertEqual(scale_nodes(b.active_material)[0].label, 'Viridis')
            # A custom shader without a usable ramp can still receive another map.
            node, _ = scale_nodes(b.active_material)
            b.active_material.node_tree.links.remove(node.inputs['Fac'].links[0])
            self.assertEqual(op(filepath=str(second), **args), {'FINISHED'})
            self.assertEqual(len(scale_nodes(b.active_material)[0].color_ramp.elements), 3)
            self.assertEqual(op(filepath=str(second), palette='magma', **args), {'FINISHED'})
            self.assertEqual(op(filepath=str(second), **dict(args, field_units='V')), {'FINISHED'})
            self.assertEqual(len(scale_nodes(b.active_material)[0].color_ramp.elements), 3)
            self.assertEqual(op(filepath=str(second), **dict(args, field_kind='unsigned', field_units='V',
                                                              minimum=0.)), {'FINISHED'})
            self.assertEqual(len(scale_nodes(b.active_material)[0].color_ramp.elements), 2)

    def test_cube_contour_resamples_edited_map_and_restores_new_base(self):
        from molecular_blender.plotter import BlendMoleculeSetup, compute_surface_vertex_sets, draw_surfaces
        from molecular_blender.scene_state import read_record
        from molecular_blender.scene_surfaces import replace_surface_geometry, surface_groups
        from molecular_blender.scene_esp import ATTRIBUTE, BASE_MESH, map_esp, clear_esp
        with tempfile.TemporaryDirectory() as tmp:
            x, y, z = .5 * np.indices((9, 9, 9)) - 2
            root = Path(tmp)
            density, field = root / 'density.cube', root / 'field.cube'
            write_cube(density, np.exp(-(x*x + y*y + z*z)))
            write_cube(field, x)
            blocks = set(bpy.data.objects) | set(bpy.data.collections)
            self.addCleanup(lambda: bpy.data.batch_remove((set(bpy.data.objects) | set(bpy.data.collections)) - blocks))
            model, styler, options = BlendMoleculeSetup(bpy.context, str(density), volume='density',
                isovalues='.35', cumulative=False, remesh=4)
            draw_surfaces(model, styler, bpy.context, options)
            owner = model.scene_owner
            record = read_record(owner)
            key, group = next(iter(surface_groups(owner, record).items()))
            map_esp(bpy.context, record, group, str(field), 0, 'V', -1., 1.)
            obj, = group['objects']
            old_mesh, mat = obj.data, obj.active_material
            old_base = obj[BASE_MESH].as_pointer()
            options['isovalues'] = [.55]
            surfaces, _ = compute_surface_vertex_sets(model, options, indexed=True)
            field.rename(root / 'missing')
            with self.assertRaises(FileNotFoundError):
                replace_surface_geometry(bpy.context, owner, record, key, group, surfaces, options, styler)
            self.assertEqual(obj.data, old_mesh)
            (root / 'missing').rename(field)
            replace_surface_geometry(bpy.context, owner, record, key, group, surfaces, options, styler)
            self.assertEqual(obj.active_material, mat)
            self.assertNotEqual(obj[BASE_MESH].as_pointer(), old_base)
            np.testing.assert_allclose([p.value for p in obj.data.attributes[ATTRIBUTE].data],
                                       [p.co.x for p in obj.data.vertices], atol=1e-5)
            base = obj[BASE_MESH]
            clear_esp(group)
            self.assertEqual(obj.data, base)
            self.assertTrue(obj.modifiers[0].show_render)

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
                    # Edit the actual shader ramp. Check its compositing legend, including
                    # asymmetric zero placement and XML-sensitive characters in the title.
                    import xml.etree.ElementTree as ET
                    surf.active_material['mb_field_label'] = 'Field A & B'
                    colors.elements[0].color = (.21404114, 0., 0., 1.)
                    path = root / f'{kind}.svg'
                    self.assertEqual(bpy.ops.molecular_blender.export_field_legend(group=key, filepath=str(path)),
                                     {'FINISHED'})
                    svg = ET.parse(path)
                    self.assertEqual(svg.find('.//{*}title').text, 'Field A & B (e/bohr^3)')
                    self.assertEqual(svg.find('.//{*}stop').get('stop-color'), '#800000')
                    ticks = {el.text: float(el.get('x')) for el in svg.findall('.//{*}text')}
                    if kind == 'signed':
                        self.assertAlmostEqual((ticks['0']-ticks['-2']) / (ticks['1']-ticks['-2']), 2/3, places=4)
                    else:
                        self.assertNotIn('0', ticks)
                    old_data = surf.data
                    with self.assertRaisesRegex(RuntimeError, 'units'):
                        bpy.ops.molecular_blender.map_esp(**dict(args, field_units=''))
                    self.assertEqual(surf.data, old_data)
                manifest = root / 'induced.json'
                shifted = root / 'polarized.cube'
                write_cube(shifted, x + 2.)
                write_cube(esp, np.full(x.shape, 2.))
                entry = dict(units='e/bohr³', linear_combination=[
                    dict(filename=shifted.name, coefficient=1.), dict(filename=esp.name, coefficient=-1.)])
                manifest.write_text(json.dumps(dict(molecules=[entry])))
                def map_difference():
                    return bpy.ops.molecular_blender.map_esp(group=key, filepath=str(manifest),
                        field_label='Induced density', minimum=-2., maximum=1.)
                self.assertEqual(map_difference(), {'FINISHED'})
                np.testing.assert_allclose([p.value for p in surf.data.attributes['mb_field'].data],
                                           [v.co.x for v in surf.data.vertices], atol=1.e-6)
                info = field_map_settings(surf.active_material)
                self.assertEqual((info['kind'], info['units'], info['dataset'], info['path']),
                                 ('signed', 'e/bohr³', 0, str(manifest)))
                legend = root / 'induced.svg'
                bpy.ops.molecular_blender.export_field_legend(group=key, filepath=str(legend))
                self.assertIn('Induced density (e/bohr³)', legend.read_text())
                mapped = surf.data
                for invalid in ([entry, entry], [dict(filename=shifted.name)]):
                    manifest.write_text(json.dumps(dict(molecules=invalid)))
                    with self.assertRaisesRegex(RuntimeError, 'one linear_combination'):
                        map_difference()
                    self.assertEqual(surf.data, mapped)

                # A second panel can take its range and edited ramp from this
                # map after both scalar sources move. Mesh values stay intact.
                from molecular_blender.field_scale import scale_values
                from molecular_blender.scene_esp import scale_nodes
                peer, style, opts = BlendMoleculeSetup(bpy.context, str(density), volume='density',
                    cube_field_type='density', isovalues='.35', cumulative=False, remesh=0)
                draw_surfaces(peer, style, bpy.context, opts)
                peer_key, = surface_groups(peer.scene_owner, read_record(peer.scene_owner))
                target, = owned_objects(peer.scene_owner, 'surface')
                bpy.context.view_layer.objects.active = peer.scene_owner
                opts = dict(group=peer_key, filepath=str(esp), field_kind='signed', field_units='e/bohr³',
                            field_label='Second panel', minimum=-1., maximum=1.)
                self.assertEqual(bpy.ops.molecular_blender.map_esp(**opts), {'FINISHED'})
                mesh_before = target.data
                duplicate = target.copy()
                duplicate.parent = None
                scene.collection.objects.link(duplicate)
                original_material = duplicate.active_material
                edit = bpy.ops.molecular_blender.edit_field_scale
                node, _ = scale_nodes(surf.active_material)
                node.color_ramp.elements.new(.2).color = (.1, .3, .2, 1.)
                node.color_ramp.interpolation = 'EASE'
                source_scale = scale_values(surf.active_material)
                esp.unlink()
                manifest.unlink()
                with self.assertRaisesRegex(RuntimeError, 'matching units'):
                    original_material['mb_field_units'] = 'V'
                    edit(group=peer_key, copy_from=surf.name)
                original_material['mb_field_units'] = 'e/bohr³'
                self.assertEqual(edit(group=peer_key, copy_from=surf.name), {'FINISHED'})
                self.assertEqual(target.data, mesh_before)
                self.assertEqual(duplicate.active_material, original_material)
                self.assertEqual(scale_values(target.active_material), source_scale)
                self.assertEqual(edit(group=peer_key, minimum=-1., maximum=3.), {'FINISHED'})
                self.assertEqual(scale_values(surf.active_material), source_scale)
                self.assertEqual(field_map_settings(target.active_material)['limits'], (-1., 3.))
                ramp = scale_nodes(target.active_material)[0].color_ramp
                np.testing.assert_allclose(ramp.evaluate(.25)[:3], [1., 1., 1.], atol=1e-6)

                output = root / 'panels'
                output.mkdir()
                stem = bpy.path.clean_name(peer.scene_owner.name + '-Cube', replace='-')
                occupied = output / (stem + '.svg')
                occupied.write_text('existing legend')
                series = bpy.ops.molecular_blender.export_surface_series
                self.assertEqual(series(directory=str(output), legends=True), {'FINISHED'})
                png, = output.glob('*.png')
                self.assertEqual(png.stem, stem + '-2')
                self.assertGreater(png.stat().st_size, 0)
                self.assertIn('Second panel (e/bohr³)', png.with_suffix('.svg').read_text())
                self.assertEqual(occupied.read_text(), 'existing legend')
                bpy.context.view_layer.objects.active = owner

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
