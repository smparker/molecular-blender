"""Cube datasets remain available after import and select by zero-based index."""
import tempfile
from pathlib import Path
import unittest
import numpy as np
from molecular_blender.importers import molecule_from_cube
from molecular_blender.constants import bohr2ang


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

    def test_rotated_grid_surface_and_atom_coincide_in_both_unit_conventions(self):
        from molecular_blender.isosurfaces import cube_isosurface
        # The field is the first grid index. At 0.5 the four corners form
        # this rectangle. Its first corner is also the atom position.
        corners = np.array([[4.5, -1.5, 3], [2.5, 0.5, 3],
                            [4.5, -1.5, 6], [2.5, 0.5, 6]])
        for count, unit in ((2, bohr2ang), (-2, 1.0)):
            with self.subTest(count=count):
                self.path.write_text(f'rotated\ngrid\n1 4 -2 3\n{count} 1 1 0\n'
                                     f'{count} -2 2 0\n{count} 0 0 3\n'
                                     '1 0 4.5 -1.5 3\n0 0 0 0 1 1 1 1\n')
                result = molecule_from_cube(str(self.path), {})
                grid = result['volume']
                np.testing.assert_allclose(result['atoms'][0]['position'], corners[0] * unit)
                np.testing.assert_allclose(grid['origin'], np.array([4, -2, 3]) * unit)
                surface, = cube_isosurface(grid['data'], grid['origin'], grid['axes'], [0.5])
                points = np.unique(np.reshape(surface['triangles'], (-1, 3)), axis=0)
                np.testing.assert_allclose(points, np.unique(corners * unit, axis=0), atol=1e-6)

    def test_reject_ambiguous_units_and_invalid_axes(self):
        for axes, error in ((['2 1 0 0', '-2 0 1 0', '2 0 0 1'], 'same sign'),
                            (['0 1 0 0', '2 0 1 0', '2 0 0 1'], 'nonzero'),
                            (['2 0 0 0', '2 0 1 0', '2 0 0 1'], 'nonzero'),
                            (['2 1 0 0', '2 .3 1 0', '2 0 0 1'], 'Sheared')):
            with self.subTest(axes=axes), self.assertRaisesRegex(ValueError, error):
                self.path.write_text('axes\ngrid\n1 0 0 0\n' + '\n'.join(axes) + '\n')
                molecule_from_cube(str(self.path), {})

    def test_left_handed_grid_preserves_surface_normals(self):
        from molecular_blender.isosurfaces import cube_isosurface
        data = np.indices((2, 2, 2), dtype=np.float32)[2].copy()
        normals = []
        for xstep in (1., -1.):
            axes = np.diag([xstep, 1., 1.])
            mesh, = cube_isosurface(data, np.zeros(3), axes, [.5])
            a, b, c = np.array(mesh['triangles']).transpose(1, 0, 2)
            normals.append(np.cross(b-a, c-a))
        self.assertTrue(np.all(np.abs(normals[0][:, 2]) > 0))
        np.testing.assert_allclose(normals[0], normals[1])

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
        np.testing.assert_allclose(volume.axis_norms, [np.sqrt(2), np.sqrt(8), 3])
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
