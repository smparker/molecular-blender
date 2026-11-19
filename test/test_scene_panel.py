"""Selection resolves metadata belonging to the correct molecule."""
import unittest
from pathlib import Path
import tempfile
from types import SimpleNamespace
from unittest.mock import Mock, patch

try:
    import bpy
except ImportError:
    bpy = None


def panel_labels(obj):
    from molecular_blender.scene_panel import MB_PT_molecule
    layout = Mock()
    for method in (layout.box, layout.column, layout.row):
        method.return_value = layout
    MB_PT_molecule.draw(SimpleNamespace(layout=layout), SimpleNamespace(object=obj))
    return [call.kwargs['text'] for call in layout.label.call_args_list]


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
                if role != 'atom':
                    self.assertIsNone(summary['selected_atom'])
            self.assertEqual(molecule_summary(second.scene_owner)['atoms'], 60)
        self.assertIsNone(molecule_summary(None))
        unrelated = bpy.data.objects.new('unrelated', None)
        self.assertIsNone(molecule_summary(unrelated))

    def test_selected_atom_shows_its_own_index_and_both_charges_after_rename(self):
        from molecular_blender.plotter import BlendMoleculeSetup
        from molecular_blender.scene_state import owned_objects
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / 'chloride.xyz'
            path.write_text('1\nchloride\nCl 0 0 0 -0.875\n')
            molecule, _, _ = BlendMoleculeSetup(bpy.context, str(path), isovalues='',
                                                perceive_bond_orders=True, molecular_charge=-1)
        atom, = owned_objects(molecule.scene_owner, 'atom')
        atom.name = 'renamed'
        labels = panel_labels(atom)
        for label in ('Selected atom: Cl', 'Atom index: 0', 'Partial charge: -0.875 e', 'Formal charge: -1 e'):
            self.assertIn(label, labels)

    def test_cube_atom_charge_has_its_own_label(self):
        from molecular_blender.plotter import BlendMoleculeSetup
        from molecular_blender.scene_state import owned_objects
        molecule, _, _ = BlendMoleculeSetup(bpy.context, 'examples/water.cube', isovalues='')
        oxygen = next(o for o in owned_objects(molecule.scene_owner, 'atom') if o['mb_element'] == 'o')
        labels = panel_labels(oxygen)
        for label in ('Selected atom: O', 'Atom index: 1', 'Cube atom charge: 0 e'):
            self.assertIn(label, labels)

    def test_browse_saved_orbitals_after_source_removed(self):
        from molecular_blender.plotter import BlendMoleculeSetup
        from molecular_blender import scene_panel
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / 'water.molden'
            source = Path('examples/water-sto-3g.molden').read_text()
            path.write_text(source.replace(' Ene=    0.603242', ' Sym= A1'))
            molecule, _, _ = BlendMoleculeSetup(bpy.context, str(path), isovalues='')
        owner = molecule.scene_owner
        saved = owner['mb_molecule']
        scene_panel.register()
        try:
            self.assertNotIn('Orbital 1', panel_labels(owner))
            owner.mb_show_orbitals = True
            first = panel_labels(owner)
            self.assertIn('Orbital 1', first)
            self.assertIn('Occ: 2', first)
            self.assertNotIn('Orbital 6', first)
            owner.mb_orbital_page = 99
            self.assertEqual(owner.mb_orbital_page, 2)
            last = panel_labels(owner)
            self.assertNotIn('Orbital 5', last)
            for label in ('Orbital 6', 'Orbital 7', 'Energy: —', 'Energy: 0.753849',
                          'Occ: 0', 'Spin: Alpha', 'Sym: A1', 'Sym: —'):
                self.assertIn(label, last)
            self.assertEqual(owner['mb_molecule'], saved)
        finally:
            scene_panel.unregister()

    def test_addon_registers_and_unregisters_panel(self):
        import addon_utils
        import molecular_blender
        self.assertIs(addon_utils.enable('molecular_blender'), molecular_blender)
        try:
            self.assertTrue(hasattr(bpy.types, 'MB_PT_molecule'))
            from molecular_blender.scene_state import normalize_duplicate_ids
            self.assertIn(normalize_duplicate_ids, bpy.app.handlers.load_post)
        finally:
            addon_utils.disable('molecular_blender')
        self.assertFalse(hasattr(bpy.types, 'MB_PT_molecule'))
