"""Saved molecule identity, duplication and settings contracts."""
import json
from pathlib import Path
import tempfile
import unittest

from molecular_blender.options import validated_options

try:
    import bpy
except ImportError:
    bpy = None


class TestSettings(unittest.TestCase):
    def test_user_settings_are_distinct_from_runtime_values(self):
        source = {'hook_atoms': 'auto', 'isovalues': '0.2,0.5'}
        saved = validated_options(source)
        self.assertEqual(saved['hook_atoms'], 'auto')
        self.assertEqual(source, {'hook_atoms': 'auto', 'isovalues': '0.2,0.5'})
        self.assertEqual(validated_options(json.loads(json.dumps(saved))), saved)

    def test_invalid_settings_fail_before_scene_creation(self):
        for options in ({'resolution': 0}, {'bond_thickness': float('nan')},
                        {'colors': 'missing'}, {'keystride': 1.5},
                        {'isovalues': 'nan'}, {'isovalues': '1.5'},
                        {'hook_atoms': True}, {'cube_dataset': -1}):
            with self.subTest(options=options), self.assertRaises(ValueError):
                validated_options(options)


@unittest.skipUnless(bpy, 'requires Blender')
class TestSceneState(unittest.TestCase):
    def setUp(self):
        self.before = set(bpy.data.objects)
        self.before_collections = set(bpy.data.collections)

    def tearDown(self):
        for obj in set(bpy.data.objects) - self.before:
            bpy.data.objects.remove(obj, do_unlink=True)
        for coll in set(bpy.data.collections) - self.before_collections:
            bpy.data.collections.remove(coll)

    def import_molecule(self, **options):
        from molecular_blender.plotter import BlendMoleculeSetup
        options = dict(dict(isovalues='', ignore_hydrogen=False, plot_style='sticks'), **options)
        return BlendMoleculeSetup(bpy.context, 'examples/water-sto-3g.molden', **options)

    def test_rename_and_duplicate_hierarchy(self):
        from molecular_blender.scene_state import (read_record, resolve_owner, owned_objects,
                                                  normalize_duplicate_ids, ID)
        molecule, _, runtime = self.import_molecule()
        owner = molecule.scene_owner
        original = read_record(owner)
        self.assertIsInstance(runtime['hook_atoms'], bool)
        self.assertEqual(original['settings']['hook_atoms'], 'auto')
        self.assertEqual(len(original['atoms']), 3)
        self.assertEqual(len(original['orbitals']), 7)
        owner.name = 'Renamed root'
        molecule.scene_collection.name = 'Renamed collection'
        for obj in owned_objects(owner):
            obj.name = 'Renamed child'
            self.assertEqual(resolve_owner(obj), owner)
        copied = owner.copy()
        molecule.scene_collection.objects.link(copied)
        for obj in owned_objects(owner):
            child = obj.copy()
            child.parent = copied
            molecule.scene_collection.objects.link(child)
        normalize_duplicate_ids()
        self.assertNotEqual(read_record(copied)['id'], read_record(owner)['id'])
        self.assertEqual({o[ID] for o in owned_objects(owner, 'atom')},
                         {o[ID] for o in owned_objects(copied, 'atom')})
        for child in owned_objects(copied):
            self.assertEqual(resolve_owner(child), copied)
        self.assertEqual(read_record(owner)['atoms'], original['atoms'])

    def test_surface_definitions_and_owner(self):
        from molecular_blender.plotter import draw_surfaces
        from molecular_blender.scene_state import read_record, owned_objects, resolve_owner
        molecule, styler, options = self.import_molecule(isovalues='0.05', cumulative=False, orbital='1', resolution=0.4, remesh=0)
        draw_surfaces(molecule, styler, bpy.context, options)
        surfaces = owned_objects(molecule.scene_owner, 'surface')
        self.assertTrue(surfaces)
        record = read_record(molecule.scene_owner)
        self.assertEqual(len(surfaces), len(record['surfaces']))
        for surface in surfaces:
            self.assertEqual(resolve_owner(surface), molecule.scene_owner)

    def test_wireframe_has_persistent_owner(self):
        from molecular_blender.scene_state import read_record
        molecule, _, _ = self.import_molecule(object_type='wireframe')
        self.assertEqual(len(read_record(molecule.scene_owner)['atoms']), 3)
        self.assertIn(molecule.scene_owner.name, molecule.scene_collection.objects)

    def test_save_reload_datablocks_keeps_identity(self):
        # Load a saved collection into this scene, avoiding destruction of the
        # runner's scene and exercising Blender's actual serialization/remapping.
        from molecular_blender.scene_state import read_record, resolve_owner, normalize_duplicate_ids
        molecule, _, _ = self.import_molecule()
        original = read_record(molecule.scene_owner)
        with tempfile.TemporaryDirectory() as temporary:
            path = str(Path(temporary) / 'molecule.blend')
            bpy.data.libraries.write(path, {molecule.scene_collection})
            with bpy.data.libraries.load(path, link=False) as (source, target):
                target.collections = [molecule.scene_collection.name]
            restored = target.collections[0]
            bpy.context.scene.collection.children.link(restored)
        owner = next(o for o in restored.all_objects if 'mb_molecule' in o)
        self.assertEqual(read_record(owner), original)
        normalize_duplicate_ids()
        self.assertNotEqual(read_record(owner)['id'], read_record(molecule.scene_owner)['id'])
        for obj in restored.all_objects:
            self.assertEqual(resolve_owner(obj), owner)

    def test_unknown_version_is_rejected(self):
        from molecular_blender.scene_state import read_record, RECORD
        molecule, _, _ = self.import_molecule()
        record = read_record(molecule.scene_owner)
        record['version'] = 999
        molecule.scene_owner[RECORD] = json.dumps(record)
        with self.assertRaisesRegex(ValueError, 'version'):
            read_record(molecule.scene_owner)
