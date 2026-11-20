"""The file browser follows the same format and dataset rules as the importer."""
import json
from pathlib import Path
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import Mock, patch

from molecular_blender import importers

try:
    import bpy
except ImportError:
    bpy = None


@unittest.skipUnless(bpy, 'requires Blender')
class TestAdaptiveImport(unittest.TestCase):
    def setUp(self):
        from molecular_blender import MolecularBlender, import_dialog
        self.ui = import_dialog
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.folder = Path(self.tmp.name)
        if not MolecularBlender.is_registered:
            bpy.utils.register_class(MolecularBlender)
            self.addCleanup(bpy.utils.unregister_class, MolecularBlender)
        props = getattr(bpy.ops, 'import').molecular_blender.get_rna_type().properties
        self.op = SimpleNamespace(**{p.identifier: p.default for p in props
                                     if p.identifier != 'dataset_choice' and hasattr(p, 'default')})
        self.op.filepath = str(self.folder / 'stale.xyz')
        self.params = SimpleNamespace(directory=bytes(self.folder), filename='')
        self.context = SimpleNamespace(space_data=SimpleNamespace(type='FILE_BROWSER', params=self.params),
                                       active_file=None)

    def visible(self, filename):
        self.context.active_file = SimpleNamespace(relative_path=filename, name=Path(filename).name)
        box = Mock()
        self.ui.draw_surfaces(SimpleNamespace(box=lambda: box), self.op, self.context)
        self.labels = [call.kwargs['text'] for call in box.label.call_args_list]
        return {call.args[1] for call in box.prop.call_args_list}

    def test_molden_fields_follow_selection(self):
        self.assertIn('orbital_spin', self.visible('ORCA.MOLDEN.INPUT'))
        self.op.orbital_choice = 'custom'
        for expression, spin in [('3, 4', False), ('beta:lumo+1', False), ('2, homo-1', True), ('-1', True)]:
            self.op.orbital = expression
            props = self.visible('wave.MOLDEN')
            self.assertEqual('orbital_spin' in props, spin)
            self.assertIn('orbital', props)
            self.assertFalse(props & {'cube_dataset', 'dataset_choice', 'cube_field_type'})
        for field in ('density', 'alpha', 'beta', 'spin'):
            self.op.volume = field
            props = self.visible('wave.molden')
            self.assertIn('volume', props)
            self.assertIn('resolution', props)
            self.assertFalse(props & {'orbital', 'orbital_choice', 'orbital_spin'})
        self.op.surface_advanced = True
        self.assertIn('molden_convention', self.visible('wave.molden'))
        for name in ('atoms.XYZ', 'mol.SD', 'mol.sdf', 'mol.mol', ''):
            self.assertFalse(self.visible(name))

    def test_orca_title_preview_including_manifest(self):
        self.op.surface_advanced = True
        source = self.folder / 'wave.molden'
        source.write_text('[Molden Format]\n[Title]\n Molden file created by orca_2mkl for BaseName=test\n[Atoms] AU\n')
        (self.folder / 'trajectory.json').write_text(json.dumps({'molecules': [{'filename': source.name}]}))
        for name in (source.name, 'trajectory.json'):
            self.visible(name)
            self.assertIn('Auto uses ORCA', self.labels)
        self.op.molden_convention = 'standard'
        self.visible(source.name)
        self.assertNotIn('Auto uses ORCA', self.labels)

    def test_cube_and_manifest_previews_track_both_files(self):
        cube = self.folder / 'fields.cub'
        manifest = self.folder / 'animation.json'
        manifest.write_text(json.dumps({'molecules': [{'filename': cube.name}]}))
        # A header alone is enough. Neither preview reads the voxel grid.
        cube.write_text('comment\ncomment\n-1 0 0 0\n2 1 0 0\n2 0 1 0\n2 0 0 1\n1 0 0 0 0\n2 7\n23\n')
        with patch.object(importers, 'molecule_from_cube', side_effect=AssertionError('full read')):
            props = self.visible(manifest.name)
            self.assertIn('dataset_choice', props)
            self.assertIn('cube_field_type', props)
            self.assertFalse(props & {'volume', 'orbital', 'resolution', 'molden_convention'})
            items = self.ui.dataset_items(self.op, self.context)
            self.assertEqual([item[1] for item in items], ['Dataset 1 (ID 7)', 'Dataset 2 (ID 23)'])
            self.op.cumulative = False
            self.assertNotIn('cube_field_type', self.visible(cube.name))
            cube.write_text('comment\ncomment\n1 0 0 0\n')
            self.assertNotIn('dataset_choice', self.visible(manifest.name))
            manifest.write_text('{"molecules":[{"filename":"frame.molden"}]}')
            self.assertIn('volume', self.visible(manifest.name))
            manifest.write_text('{')
            self.assertFalse(self.visible(manifest.name))

    def test_combination_exposes_only_result_dataset_even_after_sampling(self):
        manifest = self.folder / 'fields.json'
        manifest.write_text(json.dumps({'molecules': [
            {'filename': 'first.molden'},
            {'units': 'e/bohr^3', 'linear_combination': [{'filename': 'a.cube', 'coefficient': 1}]}]}))
        self.op.frame_start = 1
        self.assertEqual(self.ui.surface_source(str(manifest), 1), ('cube', range(1)))
        props = self.visible(manifest.name)
        self.assertIn('cube_field_type', props)
        self.assertFalse(props & {'dataset_choice', 'volume', 'orbital'})

    def test_execute_uses_confirmed_file_and_ignores_hidden_settings(self):
        from molecular_blender import MolecularBlender
        from molecular_blender.plotter import process_options
        from molecular_blender.options import validated_options
        self.op.cube_dataset = 3
        self.op.volume = 'density'
        self.op.cube_both_signs = True
        self.op.isovalues = '0.25'
        opts = MolecularBlender._build_options(self.op, ('cube', range(1)), from_browser=True)
        self.assertEqual(opts['cube_dataset'], 0)
        opts = validated_options(opts)
        process_options(opts, cube=True)
        self.assertEqual(set(opts['isovalues']), {-0.25, 0.25})
        self.op.volume = 'orbital'
        self.op.cube_both_signs = False
        opts = validated_options(MolecularBlender._build_options(self.op, ('cube', range(5)), from_browser=True))
        self.assertEqual(opts['cube_dataset'], 3)
        process_options(opts, cube=True)
        self.assertEqual(opts['isovalues'], [0.25])
        self.op.isovalues = 'stale invalid surface levels'
        opts = MolecularBlender._build_options(self.op, ('xyz', None), from_browser=True)
        self.assertEqual(validated_options(opts)['isovalues'], '')
        self.params.filename = 'typed.molden'
        self.assertEqual(self.ui.browser_path(self.op, self.context), str(self.folder / 'typed.molden'))
        self.assertEqual(self.ui.browser_path(self.op, None), self.op.filepath)
