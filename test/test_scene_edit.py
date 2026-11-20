"""Independent styles, geometry preservation and actual Blender undo/redo."""
from pathlib import Path
import subprocess
import tempfile
import unittest
from unittest.mock import patch

try:
    import bpy
except ImportError:
    bpy = None


@unittest.skipUnless(bpy, 'requires Blender')
class TestRestyle(unittest.TestCase):
    def test_make_independent_preserves_internal_sharing_and_detaches_bevel_sources(self):
        from molecular_blender.plotter import BlendMoleculeSetup
        from molecular_blender.scene_edit import make_molecule_independent
        from molecular_blender.scene_state import owned_objects
        model, _, _ = BlendMoleculeSetup(bpy.context, 'examples/bond_orders.sdf', draw_bond_orders=True,
                                        ignore_hydrogen=False, isovalues='')
        root = model.scene_owner
        other = root.copy()
        model.scene_collection.objects.link(other)
        pairs = []
        for obj in owned_objects(root):
            copy = obj.copy()
            copy.parent = other
            model.scene_collection.objects.link(copy)
            pairs.append((obj, copy))
        make_molecule_independent(other)
        groups = {}
        for original, copy in pairs:
            self.assertNotEqual(original.data, copy.data)
            self.assertEqual(groups.setdefault(original.data, copy.data), copy.data)
            self.assertNotEqual(original.active_material, copy.active_material)
            if original.get('mb_order', 1) > 1:
                self.assertNotEqual(original.data.bevel_object.data, copy.data.bevel_object.data)
        baseline = [(b.data, b.active_material) for _, b in pairs]
        make_molecule_independent(other)
        self.assertEqual(baseline, [(b.data, b.active_material) for _, b in pairs])

    def setUp(self):
        self.before = set(bpy.data.objects)
        self.collections = set(bpy.data.collections)

    def tearDown(self):
        for obj in set(bpy.data.objects) - self.before:
            bpy.data.objects.remove(obj, do_unlink=True)
        for coll in set(bpy.data.collections) - self.collections:
            bpy.data.collections.remove(coll)

    def test_linked_duplicate_is_independent_and_repeatable(self):
        from molecular_blender.plotter import BlendMoleculeSetup
        from molecular_blender.scene_state import owned_objects, normalize_duplicate_ids, read_record
        from molecular_blender.scene_edit import restyle_molecule
        molecule, _, _ = BlendMoleculeSetup(bpy.context, 'examples/bond_orders.sdf',
                                           plot_style='sticks', draw_bond_orders=True,
                                           ignore_hydrogen=False, isovalues='')
        original = molecule.scene_owner
        copied = original.copy()
        molecule.scene_collection.objects.link(copied)
        pairs = []
        for obj in owned_objects(original):
            duplicate = obj.copy()
            duplicate.parent = copied
            molecule.scene_collection.objects.link(duplicate)
            pairs.append((obj, duplicate))
        normalize_duplicate_ids()
        snapshot = {a.name: (tuple(a.scale), a.data, a.active_material) for a, _ in pairs}
        old_profiles = {a.name: a.data.bevel_object for a, _ in pairs if a.get('mb_order', 1) > 1}
        self.assertTrue(old_profiles)
        restyle_molecule(copied, atom_scale=1.5, bond_thickness=0.4, colors='vmd')
        scales = {b.name: tuple(b.scale) for _, b in pairs}
        restyle_molecule(copied, atom_scale=1.5, bond_thickness=0.4, colors='vmd')
        for original_obj, edited in pairs:
            self.assertEqual((tuple(original_obj.scale), original_obj.data, original_obj.active_material), snapshot[original_obj.name])
            self.assertEqual(tuple(edited.scale), scales[edited.name])
            if edited.get('mb_role') == 'atom':
                self.assertEqual(original_obj.data, edited.data)
                self.assertNotEqual(original_obj.active_material, edited.active_material)
            elif edited.get('mb_role') == 'bond':
                self.assertNotEqual(original_obj.data, edited.data)
                self.assertAlmostEqual(edited.data.bevel_depth, 0.4)
                if edited.get('mb_order', 1) > 1:
                    self.assertEqual(original_obj.data.bevel_object, old_profiles[original_obj.name])
                    self.assertNotEqual(edited.data.bevel_object, original_obj.data.bevel_object)
                    self.assertEqual(len(edited.data.bevel_object.data.splines), edited['mb_order'])
        self.assertEqual(read_record(copied)['settings']['colors'], 'vmd')
        self.assertEqual(read_record(original)['settings']['colors'], 'default')
        restyle_molecule(copied, surface_opacity=0.6)
        self.assertEqual(read_record(copied)['settings']['colors'], 'vmd')
        self.assertEqual(read_record(copied)['settings']['bond_thickness'], 0.4)
        self.assertEqual(read_record(copied)['settings']['atom_scale'], 1.5)

    def test_surface_opacity_preserves_geometry_and_unrelated_materials(self):
        from molecular_blender.plotter import BlendMoleculeSetup, draw_surfaces
        from molecular_blender.scene_state import owned_objects
        from molecular_blender.scene_edit import restyle_molecule
        molecule, styler, options = BlendMoleculeSetup(
            bpy.context, 'examples/fake-atom.molden', isovalues='0.05',
            cumulative=False, orbital='1', resolution=0.3, remesh=0)
        draw_surfaces(molecule, styler, bpy.context, options)
        surfaces = owned_objects(molecule.scene_owner, 'surface')
        self.assertTrue(surfaces)
        surface = surfaces[0]
        other = surface.copy()
        other.parent = None
        bpy.context.collection.objects.link(other)
        original_material = other.active_material
        vertices = [tuple(v.co) for v in surface.data.vertices]
        with patch('molecular_blender.plotter.compute_surface_vertex_sets', side_effect=AssertionError('must not recompute')):
            restyle_molecule(molecule.scene_owner, surface_opacity=0.35)
        self.assertEqual([tuple(v.co) for v in surface.data.vertices], vertices)
        self.assertEqual(other.active_material, original_material)
        self.assertNotEqual(surface.active_material, original_material)
        self.assertAlmostEqual(surface.active_material.node_tree.nodes['MB Surface Opacity'].inputs[0].default_value, 0.35)
        self.assertIsNone(original_material.node_tree.nodes.get('MB Surface Opacity'))

    def test_selected_accents_and_restore_handle_split_bonds_and_shared_atoms(self):
        from molecular_blender.plotter import BlendMoleculeSetup
        from molecular_blender import scene_edit
        from molecular_blender.scene_state import ID, owned_objects
        from molecular_blender.stylers import get_styler
        from molecular_blender.periodictable import elements
        scene_edit.register()
        self.addCleanup(scene_edit.unregister)
        model, _, _ = BlendMoleculeSetup(bpy.context, 'examples/water-sto-3g.molden',
                                         ignore_hydrogen=False, isovalues='')
        owner = model.scene_owner
        atoms = owned_objects(owner, 'atom')
        h, other_h = [o for o in atoms if o['mb_element'] == 'h']
        oxygen, = [o for o in atoms if o['mb_element'] == 'o']
        bonds = owned_objects(owner, 'bond')
        halves = [o for o in bonds if h[ID] in o['mb_endpoints']]
        self.assertEqual(len(halves), 2)
        bpy.ops.object.select_all(action='DESELECT')
        bpy.context.view_layer.objects.active = h
        h.select_set(True)
        oxygen.select_set(True)
        shared = other_h.active_material
        for value in (.1, .3):
            self.assertEqual(bpy.ops.molecular_blender.style_selected(color=(value, .2, .8)), {'FINISHED'})
            for obj in (h, oxygen, *halves):
                self.assertAlmostEqual(obj.active_material.diffuse_color[0], value)
            self.assertEqual(other_h.active_material, shared)
            self.assertTrue(all('mb_accent_color' not in o for o in bonds if o not in halves))
        scene_edit.restyle_molecule(owner, colors='vmd')
        self.assertAlmostEqual(oxygen.active_material.diffuse_color[0], .3)
        bpy.ops.object.select_all(action='DESELECT')
        halves[0].select_set(True)
        bpy.context.view_layer.objects.active = halves[0]
        self.assertEqual(bpy.ops.molecular_blender.style_selected(restore=True), {'FINISHED'})
        palette = get_styler({'colors': 'vmd'})
        for obj in halves:
            self.assertNotIn('mb_accent_color', obj)
            expected = palette.element_color(elements[obj['mb_element']])
            for actual, want in zip(obj.active_material.diffuse_color, expected):
                self.assertAlmostEqual(actual, want)
        self.assertIn('mb_accent_color', h)

    def test_invalid_style_leaves_scene_unchanged(self):
        from molecular_blender.plotter import BlendMoleculeSetup
        from molecular_blender.scene_edit import restyle_molecule
        molecule, _, _ = BlendMoleculeSetup(bpy.context, 'examples/bond_orders.sdf', isovalues='')
        before = molecule.scene_owner['mb_molecule']
        materials = set(bpy.data.materials)
        with self.assertRaises(ValueError):
            restyle_molecule(molecule.scene_owner, surface_opacity=1.1)
        self.assertEqual(molecule.scene_owner['mb_molecule'], before)
        self.assertEqual(set(bpy.data.materials), materials)

    def test_reset_restores_import_style_after_bond_and_surface_edits(self):
        from molecular_blender import scene_edit
        from molecular_blender.bond_edit import set_manual_bond
        from molecular_blender.plotter import BlendMoleculeSetup, draw_surfaces
        from molecular_blender.scene_state import owned_objects, read_record
        scene_edit.register()
        self.addCleanup(scene_edit.unregister)
        model, styler, options = BlendMoleculeSetup(bpy.context, 'examples/water-sto-3g.molden',
            plot_style='sticks', ignore_hydrogen=False, atom_scale=1.25, colors='vmd',
            bond_thickness=.13, volume='density', cumulative=False, isovalues='.05', resolution=.2, remesh=0)
        draw_surfaces(model, styler, bpy.context, options)
        owner = model.scene_owner
        bond_id = read_record(owner)['bonds'][0]['id']
        set_manual_bond(owner, 2, bond_id=bond_id)
        atoms, bonds = owned_objects(owner, 'atom'), owned_objects(owner, 'bond')
        expected = [(obj, tuple(obj.active_material.diffuse_color)) for obj in atoms + bonds]
        sizes = [atom.scale.x for atom in atoms]
        surface, = owned_objects(owner, 'surface')
        surface.hide_render = True
        bpy.ops.object.select_all(action='DESELECT')
        for obj in atoms + bonds:
            obj.select_set(True)
        bpy.context.view_layer.objects.active = atoms[0]
        bpy.ops.molecular_blender.style_selected(color=(.1, .4, .8))
        scene_edit.restyle_molecule(owner, colors='default', atom_scale=1.8, bond_thickness=.3, surface_opacity=.6)
        surface_state = surface.data, surface.active_material, read_record(owner)['surfaces']
        for _ in range(2):
            self.assertEqual(bpy.ops.molecular_blender.reset_atom_bond_style(), {'FINISHED'})
            for obj, rgba in expected:
                self.assertNotIn('mb_accent_color', obj)
                self.assertEqual(tuple(obj.active_material.diffuse_color), rgba)
            for atom, size in zip(atoms, sizes):
                self.assertAlmostEqual(atom.scale.x, size, places=5)
            self.assertTrue(all(abs(obj.data.bevel_depth - .13) < 1e-6 for obj in bonds))
            record = read_record(owner)
            self.assertEqual(record['settings']['colors'], 'vmd')
            self.assertEqual(record['bonds'][0]['order'], 2)
            self.assertTrue(record['bonds'][0]['locked'])
            self.assertEqual((surface.data, surface.active_material, record['surfaces']), surface_state)
            self.assertTrue(surface.hide_render)
            self.assertAlmostEqual(record['settings']['surface_opacity'], .6)

    def test_reset_uses_existing_saved_style_for_legacy_record(self):
        from molecular_blender.plotter import BlendMoleculeSetup
        from molecular_blender.scene_edit import restyle_molecule
        from molecular_blender.scene_state import read_record, write_record
        model, _, _ = BlendMoleculeSetup(bpy.context, 'examples/bond_orders.sdf', colors='vmd', isovalues='')
        owner = model.scene_owner
        record = read_record(owner)
        del record['initial_style']
        write_record(owner, record)
        restyle_molecule(owner, colors='default')
        restyle_molecule(owner, reset_style=True)
        self.assertEqual(read_record(owner)['settings']['colors'], 'vmd')

    def test_operator_undo_redo_in_fresh_blender(self):
        root = Path(__file__).resolve().parents[1]
        script = '''
import bpy, sys
sys.path.insert(0, ROOT)
import molecular_blender as mb
from molecular_blender.scene_state import read_record, owned_objects
mb.register()
owner = mb.BlendMoleculeSetup(bpy.context, SOURCE, isovalues='', ignore_hydrogen=False)[0].scene_owner
name = owner.name
bpy.context.view_layer.objects.active = owner
owner.select_set(True)
bpy.context.preferences.edit.use_global_undo = True
before = tuple(owned_objects(owner, 'atom')[0].scale)
bpy.ops.ed.undo_push(message='Before restyle')
assert bpy.ops.molecular_blender.restyle(atom_scale=1.7, bond_thickness=0.3, colors='vmd') == {'FINISHED'}
after = tuple(owned_objects(owner, 'atom')[0].scale)
assert before != after
bpy.ops.ed.undo_push(message='After restyle')
bpy.ops.ed.undo()
owner = bpy.data.objects[name]
assert read_record(owner)['settings']['colors'] == 'default'
assert tuple(owned_objects(owner, 'atom')[0].scale) == before
bpy.ops.ed.redo()
owner = bpy.data.objects[name]
assert read_record(owner)['settings']['colors'] == 'vmd'
assert tuple(owned_objects(owner, 'atom')[0].scale) == after
'''.replace('ROOT', repr(str(root))).replace('SOURCE', repr(str(root / 'examples/bond_orders.sdf')))
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / 'undo.py'
            path.write_text(script)
            result = subprocess.run([bpy.app.binary_path, '--background', '--factory-startup',
                                     '--python-exit-code', '1', '--python', str(path)],
                                    capture_output=True, text=True, timeout=60)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
