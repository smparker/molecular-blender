"""Manual topology, locked perception, parent transforms and scene lifecycle."""
from pathlib import Path
import subprocess
import tempfile
from types import SimpleNamespace
import unittest

from molecular_blender.bond_perception import perceive_bond_orders, BondPerceptionError

try:
    import bpy
except ImportError:
    bpy = None


class TestLockedPerception(unittest.TestCase):
    def test_resonance_assignment_respects_fixed_orders(self):
        symbols = ['c'] * 6 + ['h'] * 6
        edges = [(i, (i + 1) % 6) for i in range(6)] + [(i, i + 6) for i in range(6)]
        orders, charges = perceive_bond_orders(symbols, edges, fixed_orders={0: 2, 1: 1})
        self.assertEqual(orders[:6], [2, 1, 2, 1, 2, 1])
        self.assertEqual(sum(charges), 0)
        with self.assertRaises(BondPerceptionError):
            perceive_bond_orders(symbols, edges, fixed_orders={0: 3})

    def test_invalid_fixed_orders(self):
        for fixed in ({2: 1}, {True: 1}, {0: True}, {0: 4}, {0: 1.5}):
            with self.subTest(fixed=fixed), self.assertRaises(ValueError):
                perceive_bond_orders(['h', 'h'], [(0, 1)], fixed_orders=fixed)


@unittest.skipUnless(bpy, 'requires Blender')
class TestBondEditing(unittest.TestCase):
    def setUp(self):
        self.before = set(bpy.data.objects)
        self.collections = set(bpy.data.collections)
        self.temporary = tempfile.TemporaryDirectory()

    def tearDown(self):
        for obj in set(bpy.data.objects) - self.before:
            bpy.data.objects.remove(obj, do_unlink=True)
        for coll in set(bpy.data.collections) - self.collections:
            bpy.data.collections.remove(coll)
        self.temporary.cleanup()

    def molecule(self, text='2\nCO\nC 0 0 0\nO 1.3 0 0\n', **options):
        from molecular_blender.plotter import BlendMoleculeSetup
        path = Path(self.temporary.name) / 'atoms.xyz'
        path.write_text(text)
        opts = dict(isovalues='', plot_style='sticks', ignore_hydrogen=False)
        opts.update(options)
        return BlendMoleculeSetup(bpy.context, str(path), **opts)[0]

    def test_add_change_delete_split_bond_and_junctions(self):
        from molecular_blender.scene_state import read_record, owned_objects, ID
        from molecular_blender.bond_edit import set_manual_bond, delete_bond
        molecule = self.molecule()
        owner = molecule.scene_owner
        record = read_record(owner)
        ids = [a['id'] for a in record['atoms']]
        delete_bond(owner, atom_ids=ids)
        self.assertEqual(read_record(owner)['bonds'], [])
        self.assertEqual(owned_objects(owner, 'bond'), [])
        bond_id = set_manual_bond(owner, 2, atom_ids=ids)
        curves = owned_objects(owner, 'bond')
        self.assertEqual(len(curves), 2)
        self.assertEqual({obj[ID] for obj in curves}, {bond_id})
        self.assertEqual({(obj.data.bevel_factor_start, obj.data.bevel_factor_end) for obj in curves}, {(0, 0.5), (0.5, 1)})
        self.assertTrue(all(len(obj.data.bevel_object.data.splines) == 2 for obj in curves))
        double_scale = owned_objects(owner, 'atom')[0].scale.x
        set_manual_bond(owner, 3, bond_id=bond_id)
        self.assertTrue(all(obj['mb_order'] == 3 for obj in curves))
        self.assertTrue(all(len(obj.data.bevel_object.data.splines) == 3 for obj in curves))
        self.assertGreater(owned_objects(owner, 'atom')[0].scale.x, double_scale)
        record = read_record(owner)
        self.assertEqual(len(record['bonds']), 1)
        self.assertEqual(record['bonds'][0]['provenance'], 'manual')
        self.assertTrue(record['bonds'][0]['locked'])
        set_manual_bond(owner, 1, atom_ids=list(reversed(ids)))
        self.assertEqual(len(read_record(owner)['bonds']), 1)
        self.assertTrue(all(obj.data.bevel_object is None for obj in curves))
        delete_bond(owner, bond_id=bond_id)
        self.assertAlmostEqual(owned_objects(owner, 'atom')[0].scale.x, 1.0)

    def test_manual_order_survives_perception_and_conflict_is_atomic(self):
        from molecular_blender.scene_state import read_record, owned_objects
        from molecular_blender.bond_edit import set_manual_bond, reperceive_bonds
        molecule = self.molecule('4\nethene fragment\nC 0 0 0\nC 1.2 0 0\nH -1.05 0 0\nH 2.25 0 0\n')
        owner = molecule.scene_owner
        record = read_record(owner)
        carbons = [a['id'] for a in record['atoms'] if a['element'] == 'c']
        bond_id = set_manual_bond(owner, 3, atom_ids=carbons)
        reperceive_bonds(owner)
        bond = next(b for b in read_record(owner)['bonds'] if b['id'] == bond_id)
        self.assertEqual((bond['order'], bond['provenance'], bond['locked']), (3, 'manual', True))
        set_manual_bond(owner, 1, bond_id=bond_id)
        before = owner['mb_molecule']
        geometry = [(o.name, o['mb_order']) for o in owned_objects(owner, 'bond')]
        with self.assertRaises(BondPerceptionError):
            reperceive_bonds(owner)
        self.assertEqual(owner['mb_molecule'], before)
        self.assertEqual([(o.name, o['mb_order']) for o in owned_objects(owner, 'bond')], geometry)

    def test_cross_molecule_selection_and_invalid_orders_are_rejected(self):
        from molecular_blender.scene_state import owned_objects, read_record
        from molecular_blender.bond_edit import selected_bond, set_manual_bond
        first = self.molecule().scene_owner
        second = self.molecule().scene_owner
        selected = [owned_objects(first, 'atom')[0], owned_objects(second, 'atom')[0]]
        with self.assertRaisesRegex(ValueError, 'one molecule'):
            selected_bond(SimpleNamespace(selected_objects=selected))
        before = first['mb_molecule']
        with self.assertRaises(ValueError):
            set_manual_bond(first, 4, atom_ids=[a['id'] for a in read_record(first)['atoms']])
        self.assertEqual(first['mb_molecule'], before)
        parts = owned_objects(first, 'bond')
        owner, atoms, bond = selected_bond(SimpleNamespace(selected_objects=parts))
        self.assertEqual(owner, first)
        self.assertIsNone(atoms)
        self.assertIsNotNone(bond)

    def test_new_bond_hooks_follow_atoms_under_parent_transform(self):
        from mathutils import Vector
        from molecular_blender.scene_state import read_record, owned_objects
        from molecular_blender.bond_edit import set_manual_bond, delete_bond
        molecule = self.molecule('2\nCC\nC 0 0 0\nC 1.4 0 0\n')
        owner = molecule.scene_owner
        ids = [a['id'] for a in read_record(owner)['atoms']]
        delete_bond(owner, atom_ids=ids)
        owner.location = (3, -4, 2)
        owner.rotation_euler = (0.2, 0.3, 0.4)
        owner.scale = (2, 1.5, 0.8)
        set_manual_bond(owner, 1, atom_ids=ids)
        atoms = owned_objects(owner, 'atom')
        atoms[1].location += Vector((0.3, 0.4, 0.2))
        bpy.context.view_layer.update()
        from molecular_blender.scene_edit import restyle_molecule
        restyle_molecule(owner, atom_scale=4.0)
        bpy.context.view_layer.update()
        curve = owned_objects(owner, 'bond')[0]
        self.assertEqual(len(curve.modifiers), 2)
        inverse = owner.matrix_world.inverted()
        start, end = [inverse @ atom.matrix_world.translation for atom in atoms]
        axis = (end - start).normalized()
        evaluated = curve.evaluated_get(bpy.context.evaluated_depsgraph_get())
        mesh = evaluated.to_mesh()
        try:
            points = [inverse @ evaluated.matrix_world @ v.co for v in mesh.vertices]
            projection = [(p - start).dot(axis) for p in points]
            self.assertAlmostEqual(min(projection), 0, places=4)
            self.assertAlmostEqual(max(projection), (end - start).length, places=4)
        finally:
            evaluated.to_mesh_clear()

    def test_duplicate_edits_and_deletion_do_not_touch_original(self):
        from molecular_blender.scene_state import read_record, owned_objects, normalize_duplicate_ids
        from molecular_blender.bond_edit import set_manual_bond, delete_bond
        molecule = self.molecule()
        original = molecule.scene_owner
        ids = [a['id'] for a in read_record(original)['atoms']]
        set_manual_bond(original, 2, atom_ids=ids)
        copied = original.copy()
        molecule.scene_collection.objects.link(copied)
        for obj in owned_objects(original):
            child = obj.copy()
            child.parent = copied
            molecule.scene_collection.objects.link(child)
        normalize_duplicate_ids()
        profiles = [o.data.bevel_object for o in owned_objects(original, 'bond')]
        set_manual_bond(copied, 3, atom_ids=ids)
        delete_bond(copied, atom_ids=ids)
        self.assertEqual(len(owned_objects(original, 'bond')), 2)
        self.assertEqual([o.data.bevel_object for o in owned_objects(original, 'bond')], profiles)
        self.assertEqual(read_record(original)['bonds'][0]['order'], 2)

    def test_edit_undo_redo_and_save_reopen_in_fresh_blender(self):
        root = Path(__file__).resolve().parents[1]
        script = '''
import bpy, sys
sys.path.insert(0, ROOT)
import molecular_blender as mb
from molecular_blender.scene_state import read_record, owned_objects
mb.register()
owner = mb.BlendMoleculeSetup(bpy.context, SOURCE, isovalues='', ignore_hydrogen=False)[0].scene_owner
name = owner.name
atoms = owned_objects(owner, 'atom')
bpy.ops.object.select_all(action='DESELECT')
for obj in atoms[:2]: obj.select_set(True)
bpy.context.view_layer.objects.active = atoms[0]
bpy.context.preferences.edit.use_global_undo = True
before = owner['mb_molecule']
bpy.ops.ed.undo_push(message='Before bond edit')
assert bpy.ops.molecular_blender.set_bond(order='3') == {'FINISHED'}
after = owner['mb_molecule']
assert before != after
bpy.ops.ed.undo_push(message='After bond edit')
bpy.ops.ed.undo()
assert bpy.data.objects[name]['mb_molecule'] == before
bpy.ops.ed.redo()
assert bpy.data.objects[name]['mb_molecule'] == after
bpy.ops.wm.save_as_mainfile(filepath=BLEND)
bpy.ops.wm.open_mainfile(filepath=BLEND)
owner = bpy.data.objects[name]
assert owner['mb_molecule'] == after
assert any(b['order'] == 3 and b['locked'] for b in read_record(owner)['bonds'])
bpy.ops.object.select_all(action='DESELECT')
parts = owned_objects(owner, 'bond')
target = next(o for o in parts if o['mb_order'] == 3)
target.select_set(True)
target_id = target['mb_id']
bpy.context.view_layer.objects.active = target
bpy.ops.ed.undo_push(message='Before delete')
assert bpy.ops.molecular_blender.delete_bond() == {'FINISHED'}
assert not any(o.get('mb_id') == target_id for o in owned_objects(owner, 'bond'))
bpy.ops.ed.undo_push(message='After delete')
bpy.ops.ed.undo()
assert bpy.data.objects[name]['mb_molecule'] == after
'''.replace('ROOT', repr(str(root))).replace('SOURCE', repr(str(root / 'examples/bond_orders.sdf'))).replace('BLEND', repr(str(Path(self.temporary.name) / 'edited.blend')))
        path = Path(self.temporary.name) / 'bond_undo.py'
        path.write_text(script)
        result = subprocess.run([bpy.app.binary_path, '--background', '--factory-startup',
                                 '--python-exit-code', '1', '--python', str(path)],
                                capture_output=True, text=True, timeout=60)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
