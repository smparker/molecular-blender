"""Selection resolves metadata belonging to the correct molecule."""
import unittest
from unittest.mock import patch

try:
    import bpy
except ImportError:
    bpy = None


@unittest.skipUnless(bpy, 'requires Blender')
class TestMoleculePanel(unittest.TestCase):
    def setUp(self):
        self.before = set(bpy.data.objects)
        self.collections = set(bpy.data.collections)

    def tearDown(self):
        for obj in set(bpy.data.objects) - self.before:
            bpy.data.objects.remove(obj, do_unlink=True)
        for coll in set(bpy.data.collections) - self.collections:
            bpy.data.collections.remove(coll)

    def test_any_component_resolves_its_molecule_without_source(self):
        from molecular_blender.plotter import BlendMoleculeSetup, draw_surfaces
        from molecular_blender.scene_panel import molecule_summary
        from molecular_blender.scene_state import owned_objects
        first, styler, options = BlendMoleculeSetup(
            bpy.context, 'examples/water-sto-3g.molden', ignore_hydrogen=False,
            orbital='1', isovalues='0.05', cumulative=False, resolution=0.4, remesh=0)
        draw_surfaces(first, styler, bpy.context, options)
        second, _, _ = BlendMoleculeSetup(bpy.context, 'examples/tetracene_dimer.xyz', isovalues='')
        with patch('molecular_blender.importers.molecule_from_file', side_effect=AssertionError('must not reload')):
            for role in ('atom', 'bond', 'surface'):
                objects = owned_objects(first.scene_owner, role)
                self.assertTrue(objects)
                objects[0].name = 'renamed ' + role
                summary = molecule_summary(objects[0])
                self.assertEqual(summary['owner'], first.scene_owner)
                self.assertEqual(summary['atoms'], 3)
                self.assertEqual(len(summary['orbitals']), 7)
                self.assertEqual(summary['units'], 'Angstrom')
            self.assertEqual(molecule_summary(second.scene_owner)['atoms'], 60)
        self.assertIsNone(molecule_summary(None))
        unrelated = bpy.data.objects.new('unrelated', None)
        self.assertIsNone(molecule_summary(unrelated))

    def test_addon_registers_and_unregisters_panel(self):
        import molecular_blender
        molecular_blender.register()
        try:
            self.assertTrue(hasattr(bpy.types, 'MB_PT_molecule'))
            from molecular_blender.scene_state import normalize_duplicate_ids
            self.assertIn(normalize_duplicate_ids, bpy.app.handlers.load_post)
        finally:
            molecular_blender.unregister()
        self.assertFalse(hasattr(bpy.types, 'MB_PT_molecule'))
