"""Cube datasets remain available after import and select by zero-based index."""
import tempfile
from pathlib import Path
import unittest
import numpy as np
from molecular_blender.importers import molecule_from_cube


def cube_text(header='-1 0 0 0', labels='2 5\n9\n'):
    values = [value for z in (0, 1, 0, 1, 0, 1, 0, 1) for value in (0, 2*z)]
    return ('fields\nunits unspecified\n' + header + '\n2 1 0 0\n2 0 1 0\n2 0 0 1\n'
            '1 0 0 0 0\n' + labels + ' '.join(map(str, values)) + '\n')


class TestCubeDatasets(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.directory.cleanup)
        self.path = Path(self.directory.name) / 'fields.cube'
        self.path.write_text(cube_text())

    def test_all_fields_retained_and_default_unchanged(self):
        volume = molecule_from_cube(str(self.path), {})['volume']
        self.assertEqual(volume['labels'], [5, 9])
        self.assertEqual(volume['datasets'].shape, (2, 2, 2, 2))
        self.assertEqual(volume['dataset_index'], 0)
        np.testing.assert_array_equal(volume['data'], 0)
        selected = molecule_from_cube(str(self.path), {'cube_dataset': 1})['volume']
        np.testing.assert_array_equal(selected['data'], selected['datasets'][..., 1])
        self.assertTrue(selected['data'].flags.c_contiguous)
        np.testing.assert_array_equal(selected['data'][:, :, 1], 2)

    def test_optional_header_value_count(self):
        self.path.write_text(cube_text(header='1 0 0 0 2', labels=''))
        volume = molecule_from_cube(str(self.path), {})['volume']
        self.assertEqual(volume['labels'], [0, 1])
        self.assertEqual(volume['datasets'].shape[-1], 2)

    def test_invalid_selection_and_data_lengths(self):
        for index in (-1, 2, True, .5):
            with self.subTest(index=index), self.assertRaises(ValueError):
                molecule_from_cube(str(self.path), {'cube_dataset': index})
        for text in (cube_text().rsplit(' ', 1)[0], cube_text() + '42\n'):
            self.path.write_text(text)
            with self.assertRaisesRegex(ValueError, 'field values'):
                molecule_from_cube(str(self.path), {})

    def test_scene_model_and_surface_selection(self):
        try:
            import bpy
        except ImportError:
            self.skipTest('requires Blender')
        from molecular_blender.molecule import Molecule
        from molecular_blender.plotter import default_options, compute_surface_vertex_sets
        molecule = Molecule.from_dict('cube_fields', molecule_from_cube(str(self.path), {}))
        molecule.volume.select_dataset(1)
        self.assertEqual(molecule.volume.labels, [5, 9])
        options = default_options({'cube_dataset': 1, 'isosurfaces': True,
                                   'isovalues': [1.], 'cumulative': False})
        surfaces, _ = compute_surface_vertex_sets(molecule, options)
        self.assertGreater(len(surfaces[0]['triangles']), 0)
        self.assertEqual(molecule.volume.dataset_index, 1)
        np.testing.assert_array_equal(molecule.volume.datasets[..., 0], 0)


    def test_explicit_containment_semantics_and_voxel_volume(self):
        try:
            import bpy
        except ImportError:
            self.skipTest('requires Blender')
        from molecular_blender.molecule import VolumeData
        from molecular_blender.containing_isovalues import isovalue_containing_proportion
        data = np.arange(1, 9, dtype=np.float32).reshape(2, 2, 2)
        axes = np.array([[1., 1., 0.], [-2., 2., 0.], [0., 0., 3.]])
        volume = VolumeData(np.zeros(3), axes, [2,2,2], data)
        self.assertAlmostEqual(volume.dvol, 12.)
        for kind, square in (('auto', False), ('orbital', True), ('density', False), ('signed', False)):
            expected = isovalue_containing_proportion([.5], data, 12., square=square)
            self.assertEqual(volume.isovalue_containing_proportion([.5], kind), expected)
        self.assertNotEqual(volume.isovalue_containing_proportion([.5], 'orbital'),
                            volume.isovalue_containing_proportion([.5], 'density'))
        volume.field_type = 'orbital'
        self.assertEqual(volume.isovalue_containing_proportion([.5]),
                         volume.isovalue_containing_proportion([.5], 'orbital'))
        with self.assertRaises(ValueError):
            volume.isovalue_containing_proportion([.5], 'typo')
