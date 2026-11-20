"""Exercise a mode through the same operators used by the molecule sidebar."""
from pathlib import Path
import tempfile
import unittest

import numpy as np

try:
    import bpy
except ImportError:
    bpy = None


@unittest.skipUnless(bpy, 'requires Blender')
class TestModePreview(unittest.TestCase):
    def test_signed_mode_motion_preserves_authored_positions_and_offsets(self):
        from molecular_blender import scene_panel
        from molecular_blender.plotter import BlendMoleculeSetup
        from molecular_blender.scene_state import owned_objects
        from molecular_blender.vibrations import BASE
        objects, collections = set(bpy.data.objects), set(bpy.data.collections)
        scene = bpy.context.scene
        frame = scene.frame_current
        scene_panel.register()
        try:
            with tempfile.TemporaryDirectory() as tmp:
                path = Path(tmp) / 'stretch.molden'
                path.write_text('[Molden Format]\n[FR-COORD]\nC 0 0 0\nC 2.5 0 0\n'
                                '[FREQ]\n0\n1200\n[FR-NORM-COORD]\nvibration 1\n0 0 0\n0 0 0\n'
                                'vibration 2\n0 0 -2\n0 0 1\n')
                model, _, _ = BlendMoleculeSetup(bpy.context, str(path), isovalues='')
            owner = model.scene_owner
            atoms = sorted(owned_objects(owner, 'atom'), key=lambda obj: obj['mb_atom_index'])
            bond, = owned_objects(owner, 'bond')
            self.assertEqual({m.object for m in bond.modifiers if m.type == 'HOOK'}, set(atoms))
            bpy.context.view_layer.objects.active = atoms[0]
            atoms[0].location.y += .4
            atoms[1].delta_location.x = .1
            locations = [tuple(atom.location) for atom in atoms]
            self.assertEqual(owner.mb_vibration_mode, '2')
            owner.mb_mode_period = 16
            owner.mb_mode_amplitude = .4
            self.assertEqual(bpy.ops.molecular_blender.preview_mode(), {'FINISHED'})
            for frame_num, multiplier in [(1, 0), (5, 1), (13, -1), (21, 1)]:
                scene.frame_set(frame_num)
                np.testing.assert_allclose([a.delta_location.z for a in atoms], [-.4*multiplier, .2*multiplier],
                                           atol=1.e-6)
                self.assertAlmostEqual(atoms[1].delta_location.x, .1)
            start, end = [a.matrix_world.translation for a in atoms]
            axis = (end - start).normalized()
            evaluated = bond.evaluated_get(bpy.context.evaluated_depsgraph_get())
            try:
                mesh = evaluated.to_mesh()
                projection = [(evaluated.matrix_world @ v.co - start).dot(axis) for v in mesh.vertices]
                np.testing.assert_allclose([min(projection), max(projection)], [0., (end-start).length], atol=1e-4)
            finally:
                evaluated.to_mesh_clear()
            owner.mb_mode_phase = np.pi/2
            bpy.ops.molecular_blender.preview_mode()
            scene.frame_set(1)
            self.assertAlmostEqual(atoms[0].delta_location.z, -.4)
            self.assertEqual(bpy.ops.molecular_blender.clear_mode(), {'FINISHED'})
            for atom, location, offset in zip(atoms, locations, [(0., 0., 0.), (.1, 0., 0.)]):
                np.testing.assert_allclose(atom.delta_location, offset)
                np.testing.assert_allclose(atom.location, location)
                self.assertNotIn(BASE, atom)
            atoms[1].keyframe_insert(data_path='delta_location', frame=1)
            with self.assertRaisesRegex(RuntimeError, 'already has animated'):
                bpy.ops.molecular_blender.preview_mode()
            self.assertNotIn(BASE, atoms[0])
        finally:
            scene_panel.unregister()
            for obj in set(bpy.data.objects) - objects:
                bpy.data.objects.remove(obj, do_unlink=True)
            for coll in set(bpy.data.collections) - collections:
                bpy.data.collections.remove(coll)
            scene.frame_set(frame)
