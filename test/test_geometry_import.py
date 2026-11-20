"""Failed geometry imports must leave a styled scene ready for another attempt."""
import unittest
from unittest.mock import patch

try:
    import bpy
except ImportError:
    bpy = None


@unittest.skipUnless(bpy, 'requires Blender')
class TestGeometryRollback(unittest.TestCase):
    def setUp(self):
        self.scene = bpy.data.scenes.new('Rollback figure')
        self.original_scene = bpy.context.window.scene
        bpy.context.window.scene = self.scene
        self.before = self.ids()
        mesh = bpy.data.meshes.new('Figure geometry')
        mesh.from_pydata([(1., 2., 3.)], [], [])
        self.figure = bpy.data.objects.new('Figure', mesh)
        self.scene.collection.objects.link(self.figure)
        self.figure.select_set(True)
        bpy.context.view_layer.objects.active = self.figure
        self.figure.location = (4., 5., 6.)

    @staticmethod
    def ids():
        return {item for kind in ('objects', 'collections', 'meshes', 'curves', 'metaballs', 'materials', 'actions')
                for item in getattr(bpy.data, kind)}

    def tearDown(self):
        if bpy.context.mode != 'OBJECT':
            bpy.ops.object.mode_set(mode='OBJECT')
        bpy.context.window.scene = self.original_scene
        bpy.data.batch_remove(self.ids() - self.before)
        bpy.data.scenes.remove(self.scene)

    def test_failed_bonds_restore_edit_mode_and_allow_retry(self):
        from molecular_blender import plotter
        bpy.ops.object.mode_set(mode='EDIT')
        before = self.ids()
        bind = plotter.bind_bond_endpoints
        curves = []
        other = []
        def fail(bond, first, second):
            curves.append(bond.name)
            bind(bond, first, second)
            if len(curves) == 3:
                # An unrelated datablock created during the same call must survive.
                other.append(bpy.data.materials.new('Separate figure material'))
                raise ValueError('injected bond failure')

        args = dict(draw_bond_orders=True, plot_style='sticks', hook_atoms='on',
                    ignore_hydrogen=False, isovalues='')
        with patch.object(plotter, 'bind_bond_endpoints', fail), self.assertRaisesRegex(ValueError, 'bond failure'):
            plotter.BlendMoleculeSetup(bpy.context, 'examples/bond_orders.sdf', **args)
        self.assertEqual(len(curves), 3)
        self.assertEqual(self.ids(), before | set(other))
        self.assertEqual(bpy.context.mode, 'EDIT_MESH')
        self.assertEqual(self.figure, bpy.context.view_layer.objects.active)
        self.assertEqual(set(bpy.context.selected_objects), {self.figure})
        self.assertEqual(tuple(self.figure.location), (4., 5., 6.))
        self.assertEqual(len(self.figure.data.vertices), 1)
        model, _, _ = plotter.BlendMoleculeSetup(bpy.context, 'examples/bond_orders.sdf', **args)
        self.assertTrue(model.scene_owner.get('mb_molecule'))
        self.assertEqual(len(model.bonds), 8)

    def test_failure_after_keyframing_or_wire_conversion_restores_selection(self):
        from molecular_blender.plotter import BlendMoleculeSetup
        second = bpy.data.objects.new('Other selected object', None)
        self.scene.collection.objects.link(second)
        second.select_set(True)
        for file, settings in [('examples/benzopentyl_twist.xyz', dict(charges='scale')),
                               ('examples/water-sto-3g.molden', dict(object_type='wireframe'))]:
            with self.subTest(file=file):
                before = self.ids()
                with patch('molecular_blender.scene_state.persist_molecule', side_effect=ValueError('late failure')):
                    with self.assertRaisesRegex(ValueError, 'late failure'):
                        BlendMoleculeSetup(bpy.context, file, isovalues='', **settings)
                self.assertEqual(before, self.ids())
                self.assertEqual(set(bpy.context.selected_objects), {self.figure, second})
                self.assertEqual(bpy.context.view_layer.objects.active, self.figure)

    def test_second_molecule_uses_its_own_charges_and_recycled_ring(self):
        from molecular_blender.plotter import BlendMoleculeSetup
        args = dict(find_aromatic=True, charges='scale', recycle_materials=True, isovalues='')
        ring = bpy.data.materials.get('ring_mat') or bpy.data.materials.new('ring_mat')
        models = [BlendMoleculeSetup(bpy.context, 'examples/tetracene_dimer.xyz', **args)[0] for _ in range(2)]
        for model in models:
            rings = [obj for obj in model.scene_owner.children_recursive if obj.name.startswith(model.name + '_ring')]
            self.assertTrue(rings)
            self.assertTrue(all(obj.active_material == ring for obj in rings))
            bpy.context.view_layer.update()
            for obj in rings:
                for vertex in obj.data.vertices:
                    self.assertLess((obj.matrix_world @ vertex.co - vertex.co).length, 1e-5)
            for atom in model.atoms:
                if atom.hidden:
                    continue
                for name in (atom.name, atom.plus_charge, atom.neg_charge):
                    self.assertLess((bpy.data.objects[name].matrix_world.translation - atom.position).length, 1e-5)
            atom = model.atoms[0]
            self.assertEqual(bpy.data.objects[atom.plus_charge].active_material.name, model.pluscharge_mat)
            self.assertEqual(bpy.data.objects[atom.neg_charge].active_material.name, model.negcharge_mat)
        self.assertNotEqual(models[0].pluscharge_mat, models[1].pluscharge_mat)
        self.assertNotEqual(models[0].negcharge_mat, models[1].negcharge_mat)
