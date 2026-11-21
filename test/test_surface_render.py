"""Check rendered coverage on a fixed transparent scene using Blender's CPU renderer."""
import math
from pathlib import Path
import subprocess
import tempfile
import unittest

import numpy as np

try:
    import bpy
except ImportError:
    bpy = None


@unittest.skipUnless(bpy, 'requires Blender')
class TestSurfacePNG(unittest.TestCase):
    def test_carbonyl_strands_reach_atoms_and_switch_color_at_midpoint(self):
        from mathutils import Vector
        from molecular_blender.plotter import BlendMoleculeSetup
        from molecular_blender.scene_state import owned_objects

        previous = bpy.context.window.scene
        kinds = ('scenes', 'worlds', 'cameras', 'collections', 'objects', 'meshes', 'materials', 'curves')
        def blocks():
            return {block for kind in kinds for block in getattr(bpy.data, kind)}
        before = blocks()
        try:
            scene = bpy.data.scenes.new('Carbonyl render test')
            bpy.context.window.scene = scene
            model, _, _ = BlendMoleculeSetup(bpy.context, 'examples/bond_orders.sdf',
                draw_bond_orders=True, hook_atoms='on', colors='vmd', bond_thickness=.12, isovalues='')
            carbonyl = model.bonds[0]
            self.assertEqual(carbonyl.order, 2)
            for obj in owned_objects(model.scene_owner):
                # Hide the junction spheres so they cannot conceal shortened strands.
                obj.hide_render = obj.name not in carbonyl.name
            o, c = carbonyl.iatom.position, carbonyl.jatom.position
            camera = bpy.data.cameras.new('Along molecular plane normal')
            camera.type = 'ORTHO'
            camera.ortho_scale = 1.5 * (c - o).length
            scene.camera = bpy.data.objects.new(camera.name, camera)
            scene.collection.objects.link(scene.camera)
            scene.camera.location = (o + c) / 2 + Vector((0, 0, 4))
            scene.camera.rotation_euler.z = math.atan2(c.y - o.y, c.x - o.x)
            scene.world = bpy.data.worlds.new('Uniform light')
            scene.world.use_nodes = True
            scene.world.node_tree.nodes['Background'].inputs['Color'].default_value = (1, 1, 1, 1)
            scene.render.engine = 'CYCLES'
            scene.cycles.device = 'CPU'
            scene.cycles.samples = 8
            scene.cycles.use_denoising = False
            scene.render.film_transparent = True
            scene.view_settings.view_transform = 'Standard'
            scene.render.resolution_x, scene.render.resolution_y = 192, 96
            scene.render.resolution_percentage = 100
            with tempfile.TemporaryDirectory() as tmp:
                scene.render.filepath = str(Path(tmp) / 'carbonyl.png')
                bpy.ops.render.render(write_still=True)
                image = bpy.data.images.load(scene.render.filepath)
                try:
                    pixels = np.array(image.pixels[:]).reshape(96, 192, 4)
                finally:
                    bpy.data.images.remove(image)

            # The source O and C project to x=32 and 160. Find the two strands
            # independently from the PNG alpha at a point inside the oxygen half.
            rows = np.flatnonzero(pixels[:, 60, 3] > .75)
            strands = np.split(rows, np.flatnonzero(np.diff(rows) > 1) + 1)
            self.assertEqual(len(strands), 2)
            for strand in strands:
                self.assertGreater(len(strand), 6)
                rgb = pixels[strand[len(strand) // 2]]
                filled = np.flatnonzero(rgb[:, 3] > .75)
                self.assertAlmostEqual(filled[0], 32, delta=1)
                self.assertAlmostEqual(filled[-1], 159, delta=1)
                self.assertTrue(np.all(np.diff(filled) == 1))
                self.assertTrue(np.all(rgb[34:95, 0] > 1.15 * rgb[34:95, 1]))   # oxygen red
                self.assertTrue(np.all(rgb[97:158, 1] > 1.15 * rgb[97:158, 0]))  # carbon cyan
            gap = (strands[0][-1] + strands[1][0]) // 2
            self.assertTrue(np.all(pixels[gap, 34:158, 3] < .25))
        finally:
            bpy.context.window.scene = previous
            bpy.data.batch_remove(blocks() - before)

    def test_gradient_arrow_projects_from_atom_to_tip(self):
        from mathutils import Vector
        from molecular_blender.plotter import BlendMoleculeSetup
        previous = bpy.context.window.scene
        kinds = ('scenes', 'cameras', 'collections', 'objects', 'meshes', 'materials', 'actions')
        def blocks():
            return {item for kind in kinds for item in getattr(bpy.data, kind)}
        before = blocks()
        try:
            scene = bpy.data.scenes.new('Gradient render test')
            bpy.context.window.scene = scene
            with tempfile.TemporaryDirectory() as directory:
                xyz = Path(directory) / 'gradient.xyz'
                xyz.write_text(''.join(f'3\nvector\nC {line}\nO 3 0 0 0 0 0 0\nH 0 1 0 0 2 0 0\n'
                                       for line in ('1 2 0 0 2 0 0', '2 3 0 0 0 -1 0', '3 4 0 0 0 0 0')))
                model, _, _ = BlendMoleculeSetup(bpy.context, str(xyz), gradient=True, bond_thickness=0,
                                                  trajectory_interpolation='hold', keystride=2, isovalues='')
                arrow, = [o for o in model.scene_collection.objects if o.name.endswith('_gradient')]
                for atom in model.atoms:
                    if not atom.hidden:
                        bpy.data.objects[atom.name].hide_render = True
                model.scene_owner.location.x += .5
                camera = bpy.data.cameras.new('Gradient camera')
                scene.camera = bpy.data.objects.new(camera.name, camera)
                scene.collection.objects.link(scene.camera)
                scene.camera.location = (2.5, 2, 5)
                camera.type = 'ORTHO'
                camera.ortho_scale = 3
                scene.render.engine = 'CYCLES'
                scene.cycles.samples = 4
                scene.cycles.device = 'CPU'
                scene.cycles.use_denoising = False
                scene.render.film_transparent = True
                scene.render.resolution_x, scene.render.resolution_y = 192, 64
                scene.render.resolution_percentage = 100
                scene.render.filepath = str(Path(directory) / 'gradient.png')
                scene.frame_set(1)
                bpy.ops.render.render(write_still=True)
                image = bpy.data.images.load(scene.render.filepath)
                try:
                    mask = np.array(image.pixels[:]).reshape(64, 192, 4)[:, :, 3] > .5
                finally:
                    bpy.data.images.remove(image)
                # Source and tip project to x=32 and x=160. The head widens, then tapers.
                filled = np.flatnonzero(mask.any(axis=0))
                self.assertAlmostEqual(filled[0], 32, delta=1)
                self.assertAlmostEqual(filled[-1], 159, delta=3)
                shaft, head, near_tip = mask[:, [70, 127, 151]].sum(axis=0)
                self.assertGreater(shaft, 2)
                self.assertGreater(head, 1.7 * shaft)
                self.assertLess(near_tip, shaft)
                scene.frame_set(3)
                np.testing.assert_allclose(arrow.matrix_world.translation, (2.5, 3, 0), atol=1e-5)
                np.testing.assert_allclose(arrow.matrix_world @ Vector((0, 0, 1)), (2.5, 2, 0), atol=1e-5)
                scene.frame_set(5)
                self.assertEqual(tuple(arrow.scale), (0, 0, 0))
        finally:
            bpy.context.window.scene = previous
            bpy.data.batch_remove(blocks() - before)

    def test_translated_parent_and_two_distinct_alpha_masks(self):
        from molecular_blender import scene_panel
        from molecular_blender.plotter import BlendMoleculeSetup, create_surface_meshes
        from molecular_blender.scene_state import owned_objects
        window = bpy.context.window
        previous_scene = window.scene
        scene = bpy.data.scenes.new('Surface render test')
        window.scene = scene
        objects = set(bpy.data.objects)
        collections = set(bpy.data.collections)
        camera = bpy.data.cameras.new('Test orthographic projection')
        scene_panel.register()
        try:
            model, styler, options = BlendMoleculeSetup(bpy.context, 'examples/water-sto-3g.molden',
                orbital='4,5', cumulative=False, isovalues='.05', remesh=0)
            owner = model.scene_owner
            for obj in owned_objects(owner):
                obj.hide_render = True
            # Source triangles centered left and right. Move the parent up by one unit.
            definitions = [dict(name=name, isovalue=.05, triangles=[[(x-.6,-.5,0), (x+.6,-.5,0), (x,.6,0)]])
                           for name, x in (('4', -1.), ('5', 1.))]
            create_surface_meshes(model, styler, bpy.context, options, definitions, {})
            # The normal outer-lobe shader is transparent face-on. Give the
            # fixture opaque faces so alpha measures geometry independently of lighting.
            for item in definitions:
                item['material'].use_nodes = False
            owner.location.y += 1.
            scene.camera = bpy.data.objects.new('Render camera', camera)
            scene.collection.objects.link(scene.camera)
            scene.camera.location = (0, 0, 6)
            camera.type = 'ORTHO'
            camera.ortho_scale = 5
            scene.render.engine = 'CYCLES'
            scene.cycles.device = 'CPU'
            scene.cycles.samples = 1
            scene.cycles.use_denoising = False
            scene.render.resolution_x = scene.render.resolution_y = 48
            scene.render.resolution_percentage = 100
            # Exercise the ordinary-render preset and the temporary batch settings independently.
            self.assertEqual(bpy.ops.molecular_blender.transparent_png(), {'FINISHED'})
            self.assertTrue(scene.render.film_transparent)
            self.assertEqual(scene.render.image_settings.color_mode, 'RGBA')
            scene.render.film_transparent = False
            scene.render.image_settings.file_format = 'JPEG'
            with tempfile.TemporaryDirectory() as directory:
                self.assertEqual(bpy.ops.molecular_blender.export_surface_series(directory=directory), {'FINISHED'})
                files = sorted(Path(directory).glob('*.png'))
                self.assertEqual(len(files), 2)
                for i, path in enumerate(files):
                    image = bpy.data.images.load(str(path))
                    try:
                        self.assertEqual(list(image.size), [48, 48])
                        alpha = np.array(image.pixels[:]).reshape(48, 48, 4)[:, :, 3]
                        y, x = np.where(alpha > .5)
                        self.assertGreater(len(x), 30)
                        self.assertLess(len(x), 100)
                        self.assertAlmostEqual(x.mean(), 14 if i == 0 else 33, delta=2)
                        self.assertGreater(y.mean(), 29)
                        self.assertEqual(alpha[0, 0], 0)
                    finally:
                        bpy.data.images.remove(image)
                if hasattr(bpy.utils, 'register_cli_command'):
                    # Start another Blender on the saved scene and use the checkout's CLI.
                    blend = Path(directory) / 'styled figure.blend'
                    bpy.data.libraries.write(str(blend), {scene})
                    checkout = str(Path(__file__).resolve().parents[1])
                    register = (f'import sys; sys.path.insert(0, {checkout!r}); '
                                'from molecular_blender import render_cli; render_cli.register()')
                    output = Path(directory) / 'orbital panels'
                    cmd = [bpy.app.binary_path, '-b', '--factory-startup', str(blend), '--python-expr', register,
                           '-c', 'molecular_render']
                    listing = subprocess.run(cmd + ['--list'], capture_output=True, text=True, timeout=90)
                    self.assertEqual(listing.returncode, 0, listing.stdout + listing.stderr)
                    self.assertIn('Camera: ' + scene.camera.name, listing.stdout)
                    self.assertIn('Molecule: ' + owner.name, listing.stdout)
                    from molecular_blender.scene_state import read_record
                    for surface in read_record(owner)['surfaces']:
                        self.assertIn(f"Orbital {surface['field']}  ({surface['group']})", listing.stdout)
                    self.assertFalse(output.exists())
                    cmd += ['--molecule', owner.name, '--output', str(output)]
                    def run(label):
                        return subprocess.run(cmd + ['--surface', label], capture_output=True, text=True,
                                              timeout=90)
                    result = run('Orbital 5')
                    self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
                    self.assertEqual([p.name for p in output.glob('*.png')], [files[1].name])
                    result = run('Missing group')
                    self.assertNotEqual(result.returncode, 0)
                    self.assertIn('Available surfaces', result.stderr + result.stdout)
            self.assertFalse(scene.render.film_transparent)
            self.assertEqual(scene.render.image_settings.file_format, 'JPEG')
        finally:
            scene_panel.unregister()
            window.scene = previous_scene
            for obj in set(bpy.data.objects) - objects:
                bpy.data.objects.remove(obj, do_unlink=True)
            for coll in set(bpy.data.collections) - collections:
                bpy.data.collections.remove(coll)
            bpy.data.scenes.remove(scene)
            bpy.data.cameras.remove(camera)
