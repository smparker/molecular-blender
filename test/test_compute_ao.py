"""Tests for compute_ao implementations (pure Python and Cython)"""

import os
import unittest
import numpy as np

import molecular_blender.compute_ao_py as ao_py

try:
    import molecular_blender.compute_ao_cy as ao_cy
    HAS_CYTHON_AO = True
except ImportError:
    HAS_CYTHON_AO = False

from molecular_blender.constants import ang2bohr

EXAMPLES = os.path.join(os.path.dirname(__file__), "../examples")


# -- gamma2n (python only) --

class TestShellBounds(unittest.TestCase):
    def test_s_shell_has_nonzero_bounds(self):
        shell = ao_py.Shell(np.zeros(3), 0, [0.2], [1.0])
        threshold = 1e-5
        radius = shell.bounding_box_size(threshold)
        self.assertGreater(radius, 0)
        expected = np.sqrt(np.log(shell.mxnorm / threshold) / 0.2)
        self.assertAlmostEqual(radius, expected)

    @unittest.skipUnless(HAS_CYTHON_AO, 'requires compiled AO backend')
    def test_bounds_match_cython_for_s_through_g(self):
        for angular in range(5):
            with self.subTest(angular=angular):
                center = np.zeros(3, dtype=np.float32)
                exponents = np.array([0.2], dtype=np.float32)
                coefficients = np.array([1.0], dtype=np.float32)
                python = ao_py.Shell(center, angular, exponents, coefficients)
                cython = ao_cy.Shell(center, angular, exponents, coefficients, 0)
                self.assertAlmostEqual(python.bounding_box_size(1e-5),
                                       cython.bounding_box_size(1e-5), places=5)

class TestGamma2n(unittest.TestCase):
    def test_gamma2n_1(self):
        assert abs(ao_py.gamma2n_impl(1) - np.sqrt(np.pi)) < 1e-10

    def test_gamma2n_2(self):
        assert ao_py.gamma2n_impl(2) == 1

    def test_gamma2n_4(self):
        assert abs(ao_py.gamma2n_impl(4) - 1.0) < 1e-10

    def test_gamma2n_negative_raises(self):
        with self.assertRaises(Exception):
            ao_py.gamma2n_impl(0)

    def test_gamma2n_array(self):
        result = ao_py.gamma2n(np.array([1.0, 2.0, 4.0]))
        assert abs(result[0] - np.sqrt(np.pi)) < 1e-10
        assert abs(result[1] - 1.0) < 1e-10


# -- gaussian_overlap (python only) --

class TestGaussianOverlap(unittest.TestCase):
    def test_s_overlap(self):
        center = np.array([0.0, 0.0, 0.0])
        olap = ao_py.gaussian_overlap(center, 1.0, [0, 0, 0], center, 1.0, [0, 0, 0])
        expected = (np.pi / 2.0) ** 1.5
        assert abs(olap - expected) < 1e-5

    def test_p_overlap(self):
        center = np.array([0.0, 0.0, 0.0])
        olap = ao_py.gaussian_overlap(center, 1.0, [1, 0, 0], center, 1.0, [1, 0, 0])
        expected = (np.pi / 2.0) ** 1.5 / 4.0
        assert abs(olap - expected) < 1e-5

    def test_overlap_symmetry(self):
        center = np.array([0.0, 0.0, 0.0])
        olap_xy = ao_py.gaussian_overlap(center, 1.5, [1, 0, 0], center, 2.0, [0, 1, 0])
        olap_yx = ao_py.gaussian_overlap(center, 2.0, [0, 1, 0], center, 1.5, [1, 0, 0])
        assert abs(olap_xy - olap_yx) < 1e-10

    def test_odd_angular_zero(self):
        center = np.array([0.0, 0.0, 0.0])
        olap = ao_py.gaussian_overlap(center, 1.0, [1, 0, 0], center, 1.0, [0, 0, 0])
        assert olap == 0.0


# -- polynomial functions (python only) --

class TestPolynomials(unittest.TestCase):
    def test_s(self):
        result = ao_py.polynomial_s(1.0, 2.0, 3.0)
        assert result[0] == 1.0

    def test_p(self):
        result = ao_py.polynomial_p(1.0, 2.0, 3.0)
        np.testing.assert_array_equal(result, [1.0, 2.0, 3.0])

    def test_d(self):
        result = ao_py.polynomial_d(1.0, 0.0, 0.0)
        assert result[0] == 1.0  # xx
        assert result[1] == 0.0  # yy

    def test_f(self):
        result = ao_py.polynomial_f(1.0, 0.0, 0.0)
        assert result[0] == 1.0  # xxx
        assert result[1] == 0.0  # yyy

    def test_g(self):
        result = ao_py.polynomial_g(1.0, 0.0, 0.0)
        assert result[0] == 1.0  # xxxx
        assert result[1] == 0.0  # yyyy

    def test_plane_polynomial_s(self):
        xx = np.array([[1.0, 2.0], [3.0, 4.0]])
        yy = np.array([[1.0, 2.0], [3.0, 4.0]])
        result = ao_py.plane_polynomial_s(xx, yy, 0.0)
        assert result.shape == (1, 2, 2)
        np.testing.assert_array_equal(result[0], np.ones((2, 2)))

    def test_plane_polynomial_p(self):
        xx = np.array([[1.0, 2.0], [3.0, 4.0]])
        yy = np.array([[5.0, 6.0], [7.0, 8.0]])
        result = ao_py.plane_polynomial_p(xx, yy, 1.0)
        assert result.shape[0] == 3
        np.testing.assert_array_equal(result[0], xx)
        np.testing.assert_array_equal(result[1], yy)


# -- Shell class (parametrized over implementations) --

class _TestShellMixin:
    def _make_s_shell(self):
        center = np.array([0.0, 0.0, 0.0], dtype=np.float32)
        exponents = np.array([1.0], dtype=np.float32)
        coeffs = np.array([1.0], dtype=np.float32)
        return self.ao.Shell(center, 0, exponents, coeffs, start=0)

    def _make_p_shell(self):
        center = np.array([0.0, 0.0, 0.0], dtype=np.float32)
        exponents = np.array([1.0], dtype=np.float32)
        coeffs = np.array([1.0], dtype=np.float32)
        return self.ao.Shell(center, 1, exponents, coeffs, start=0)

    def test_s_shell_size(self):
        sh = self._make_s_shell()
        assert sh.size == 1

    def test_p_shell_size(self):
        sh = self._make_p_shell()
        assert sh.size == 3

    def test_d_shell_size(self):
        center = np.array([0.0, 0.0, 0.0], dtype=np.float32)
        sh = self.ao.Shell(center, 2, np.array([1.0], dtype=np.float32),
                           np.array([1.0], dtype=np.float32), start=0)
        assert sh.size == 6

    def test_bounding_box_size_p(self):
        sh = self._make_p_shell()
        bb = sh.bounding_box_size(1e-6)
        assert bb > 0

    def test_large_coefficient_preserves_s_p_d_tails(self):
        for l in range(3):
            with self.subTest(l=l):
                shell = self.ao.Shell([0., 0., 0.], l, [.2], [1.], 0, denorm=False)
                coeff = np.zeros(shell.size, dtype=np.float32)
                coeff[0] = 64  # first Cartesian term is x^l
                orbital = self.ao.OrbitalCalculater([shell], coeff)
                reference = 64 * 9**l * np.exp(-.2 * 9**2)
                self.assertGreater(reference, 1e-6)
                np.testing.assert_allclose(orbital.value(9., 0., 0.), reference, rtol=3e-6)


class TestShellPy(_TestShellMixin, unittest.TestCase):
    def setUp(self):
        self.ao = ao_py

    def test_s_value_at_origin_py(self):
        """s-shell value at origin should equal the denormalized coefficient."""
        sh = self._make_s_shell()
        val = sh.value(0.0, 0.0, 0.0)
        expected = 1.0 / np.sqrt((np.pi / 2.0) ** 1.5)
        assert abs(val[0] - expected) < 1e-4

    def test_s_value_decays_py(self):
        """s-shell value should decay with distance from center."""
        sh = self._make_s_shell()
        val_origin = sh.value(0.0, 0.0, 0.0)
        val_near = sh.value(0.5, 0.0, 0.0)
        assert abs(val_near[0]) < abs(val_origin[0])

    def test_p_value_at_origin_py(self):
        """p-shell value at its center should be zero (all angular terms vanish)."""
        sh = self._make_p_shell()
        val = sh.value(0.0, 0.0, 0.0)
        for v in val:
            assert abs(v) < 1e-10

    def test_p_value_symmetry_py(self):
        """p_x at (d,0,0) should equal p_y at (0,d,0) and p_z at (0,0,d)."""
        sh = self._make_p_shell()
        d = 0.5
        val_x = sh.value(d, 0.0, 0.0)
        val_y = sh.value(0.0, d, 0.0)
        val_z = sh.value(0.0, 0.0, d)
        assert abs(val_x[0] - val_y[1]) < 1e-6
        assert abs(val_x[0] - val_z[2]) < 1e-6

    def test_s_value_far_away_py(self):
        sh = self._make_s_shell()
        val = sh.value(100.0, 100.0, 100.0)
        assert val is None


@unittest.skipUnless(HAS_CYTHON_AO, "Cython compute_ao not available")
class TestShellCy(_TestShellMixin, unittest.TestCase):
    def setUp(self):
        self.ao = ao_cy

    def test_s_value_at_origin_cy(self):
        """Cython s-shell value at origin should match Python implementation."""
        center = np.array([0.0, 0.0, 0.0], dtype=np.float32)
        sh_cy = ao_cy.Shell(center, 0, np.array([1.0], dtype=np.float32),
                            np.array([1.0], dtype=np.float32), start=0)
        sh_py = ao_py.Shell(center, 0, np.array([1.0], dtype=np.float32),
                            np.array([1.0], dtype=np.float32), start=0)
        coeff = np.array([1.0], dtype=np.float32)
        val_cy = sh_cy.value(0.0, 0.0, 0.0, coeff)
        val_py = sh_py.value(0.0, 0.0, 0.0)
        assert abs(val_cy - val_py[0]) < 1e-4

    def test_s_value_far_away_cy(self):
        """Cython returns 0.0 instead of None for screened-out shells."""
        center = np.array([0.0, 0.0, 0.0], dtype=np.float32)
        sh = ao_cy.Shell(center, 0, np.array([1.0], dtype=np.float32),
                         np.array([1.0], dtype=np.float32), start=0)
        coeff = np.array([1.0], dtype=np.float32)
        val = sh.value(100.0, 100.0, 100.0, coeff)
        assert val is None or val == 0.0


# -- Shell plane_values (python only, different API in cython) --

class TestShellPlaneValues(unittest.TestCase):
    def test_s_plane_matches_pointwise(self):
        """plane_values should be consistent with value() at each grid point."""
        center = np.array([0.0, 0.0, 0.0], dtype=np.float32)
        sh = ao_py.Shell(center, 0, np.array([1.0], dtype=np.float32),
                         np.array([1.0], dtype=np.float32), start=0)
        xs = np.array([0.0, 0.5, 1.0], dtype=np.float32)
        ys = np.array([0.0, 0.5], dtype=np.float32)
        xx, yy = np.meshgrid(xs, ys, indexing='ij')
        z = 0.0
        result = sh.plane_values(xx, yy, z)
        assert result.shape == (1, 3, 2)
        for i, x in enumerate(xs):
            for j, y in enumerate(ys):
                point_val = sh.value(float(x), float(y), z)
                assert point_val is not None
                assert abs(result[0, i, j] - point_val[0]) < 1e-5

    def test_p_plane_matches_pointwise(self):
        """p-shell plane_values should match pointwise value() calls."""
        center = np.array([0.0, 0.0, 0.0], dtype=np.float32)
        sh = ao_py.Shell(center, 1, np.array([1.0], dtype=np.float32),
                         np.array([1.0], dtype=np.float32), start=0)
        xs = np.array([-0.5, 0.0, 0.5], dtype=np.float32)
        ys = np.array([-0.5, 0.5], dtype=np.float32)
        xx, yy = np.meshgrid(xs, ys, indexing='ij')
        z = 0.3
        result = sh.plane_values(xx, yy, z)
        assert result.shape == (3, 3, 2)
        for i, x in enumerate(xs):
            for j, y in enumerate(ys):
                point_val = sh.value(float(x), float(y), z)
                if point_val is not None:
                    for k in range(3):
                        assert abs(result[k, i, j] - point_val[k]) < 1e-5

    def test_far_plane_is_none(self):
        center = np.array([0.0, 0.0, 0.0], dtype=np.float32)
        sh = ao_py.Shell(center, 0, np.array([10.0], dtype=np.float32),
                         np.array([1.0], dtype=np.float32), start=0)
        xx = np.full((2, 2), 100.0, dtype=np.float32)
        yy = np.full((2, 2), 100.0, dtype=np.float32)
        result = sh.plane_values(xx, yy, 100.0)
        assert result is None


# -- OrbitalCalculater (python only; cython doesn't export MOData) --

class TestOrbitalCalculater(unittest.TestCase):
    def setUp(self):
        from molecular_blender.importers import molecule_from_file
        geo = molecule_from_file(os.path.join(EXAMPLES, "water-sto-3g.molden"), {})
        modata = ao_py.MOData.from_dict(geo)
        self.orbital = modata.get_orbital(2)

    def test_value_at_point(self):
        val = self.orbital.value(0.0, 0.0, 0.0)
        assert isinstance(val, float)

    def test_bounding_box(self):
        p0, p1 = self.orbital.bounding_box()
        assert all(p0[i] < p1[i] for i in range(3))

    def test_box_values_shape(self):
        xyz = np.linspace(-2, 2, 5, dtype=np.float32)
        box = self.orbital.box_values(xyz, xyz, xyz)
        assert box.shape == (5, 5, 5)

    def test_plane_values_shape(self):
        xy = np.linspace(-2, 2, 5, dtype=np.float32)
        plane = self.orbital.plane_values(xy, xy, 0.0)
        assert plane.shape == (5, 5)


# -- DensityCalculater (python only) --

class TestDensityCalculater(unittest.TestCase):
    def setUp(self):
        from molecular_blender.importers import molecule_from_file
        geo = molecule_from_file(os.path.join(EXAMPLES, "fake-atom.molden"), {})
        modata = ao_py.MOData.from_dict(geo)
        self.density = modata.get_density()

    def test_value_at_point(self):
        val = self.density.value(0.0, 0.0, 0.0)
        assert val > 0

    def test_bounding_box(self):
        p0, p1 = self.density.bounding_box()
        assert all(p0[i] < p1[i] for i in range(3))

    def test_box_values_shape(self):
        xyz = np.linspace(-5, 5, 8, dtype=np.float32)
        box = self.density.box_values(xyz, xyz, xyz)
        assert box.shape == (8, 8, 8)

    def test_plane_values_shape(self):
        xy = np.linspace(-5, 5, 8, dtype=np.float32)
        plane = self.density.plane_values(xy, xy, 0.0)
        assert plane.shape == (8, 8)


def diffuse_fields(backend):
    alpha = np.array([1., 0., 0., 0.], dtype=np.float32)
    beta = np.array([0., 1., 0., 0.], dtype=np.float32)
    spin = np.outer(alpha, alpha) - .5 * np.outer(beta, beta)
    shells = [backend.Shell(np.array([x, 0, 0], dtype=np.float32), angular,
                            np.array([exp], dtype=np.float32), np.ones(1, dtype=np.float32), start)
              for x, angular, exp, start in ((-1., 0, .055, 0), (1., 1, .08, 1))]
    return [backend.OrbitalCalculater(shells, beta), backend.DensityCalculater(shells, spin, 1.5)]


class TestFractionalDensity(unittest.TestCase):
    def test_public_calculator_aliases(self):
        from molecular_blender import compute_ao
        self.assertIs(compute_ao.DensityCalculator, compute_ao.DensityCalculater)
        self.assertIs(compute_ao.OrbitalCalculator, compute_ao.OrbitalCalculater)

    @unittest.skipUnless(HAS_CYTHON_AO, 'requires compiled AO backend')
    def test_diffuse_open_shell_grid_and_containment_match(self):
        axes = [np.linspace(-10, 10, size, dtype=np.float32) for size in (17, 15, 13)]
        results = []
        for backend in (ao_py, ao_cy):
            fields = diffuse_fields(backend)
            result = []
            for field in fields:
                grid = field.box_values(*axes)
                self.assertEqual(grid.shape, (17, 15, 13))
                self.assertLess(grid.min(), 0)
                self.assertGreater(grid.max(), 0)
                np.testing.assert_allclose(grid[:, :, 5], field.plane_values(*axes[:2], axes[2][5]), atol=1e-6)
                self.assertAlmostEqual(grid[6, 7, 5], field.value(axes[0][6], axes[1][7], axes[2][5]), places=6)
                self.assertAlmostEqual(grid[0, 0, 0], field.value(axes[0][0], axes[1][0], axes[2][0]), places=7)
                levels = field.isovalue_containing_proportion([.5, -.5], resolution=.7,
                                                              box=([-10.]*3, [10.]*3))
                result.append((grid, levels))
            results.append(result)
        for (grid_a, levels_a), (grid_b, levels_b) in zip(*results):
            np.testing.assert_allclose(grid_a, grid_b, atol=2e-6)
            np.testing.assert_allclose(levels_a, levels_b, rtol=2e-4, atol=1e-6)

    @unittest.skipUnless(HAS_CYTHON_AO, 'requires compiled AO backend')
    def test_signed_lobes_match_on_finished_blender_meshes(self):
        try:
            import bpy
            from molecular_blender import marching_cube_cy
        except ImportError:
            self.skipTest('requires Blender and the compiled mesher')
        from mathutils.kdtree import KDTree
        from unittest.mock import patch
        from molecular_blender import isosurfaces, marching_cube_py
        from molecular_blender.plotter import create_geometry, create_surface_mesh_data, field_weld_distance
        resolution = .23
        results = []
        for backend, marcher in ((ao_py, marching_cube_py), (ao_cy, marching_cube_cy)):
            shapes = []
            for field, level, method in zip(diffuse_fields(backend), (.038, .0015), ('tiled', 'adaptive')):
                with patch.object(isosurfaces, 'marching_cube_box', marcher.marching_cube_box):
                    surfaces = isosurfaces.molden_isosurface(field, [level, -level], resolution, method=method)
                for surface in surfaces:
                    mesh = create_surface_mesh_data('Parity', *create_geometry(surface['triangles']),
                                                    field_weld_distance(resolution))
                    try:
                        points = [v.co.copy() for v in mesh.vertices]
                        self.assertTrue(points, (backend, method, surface['isovalue']))
                        edges = np.array([p.vertices[:] for p in mesh.polygons])[:, [0, 1, 1, 2, 2, 0]]
                        edges, counts = np.unique(np.sort(edges.reshape(-1, 2), axis=1), axis=0, return_counts=True)
                        self.assertTrue(np.all(counts == 2), 'open or nonmanifold signed lobe')
                        adjacent = [set() for _ in points]
                        for u, v in edges:
                            adjacent[u].add(v)
                            adjacent[v].add(u)
                        unseen, components = set(range(len(points))), 0
                        while unseen:
                            components += 1
                            pending = [unseen.pop()]
                            while pending:
                                neighbors = adjacent[pending.pop()] & unseen
                                unseen.difference_update(neighbors)
                                pending.extend(neighbors)
                        shapes.append((points, components, len(mesh.polygons)))
                    finally:
                        bpy.data.meshes.remove(mesh)
            results.append(shapes)
        for (a, parts_a, faces_a), (b, parts_b, faces_b) in zip(*results):
            self.assertEqual((parts_a, faces_a), (parts_b, faces_b))
            for source, target in ((a, b), (b, a)):
                tree = KDTree(len(target))
                for i, point in enumerate(target):
                    tree.insert(point, i)
                tree.balance()
                self.assertLess(max(tree.find(p)[2] for p in source), 2e-5)

    @unittest.skipUnless(HAS_CYTHON_AO, 'requires compiled AO backend')
    def test_fractional_electron_count_and_density_match_python(self):
        coefficients = np.array([1.0], dtype=np.float32)
        density = np.array([[0.5]], dtype=np.float32)
        results = []
        xyz = np.linspace(-4, 4, 21, dtype=np.float32)
        for backend in (ao_py, ao_cy):
            shell = backend.Shell(np.zeros(3, dtype=np.float32), 0,
                                  coefficients, coefficients, 0)
            calculator = backend.DensityCalculater([shell], density, 0.5)
            self.assertEqual(calculator.nelec, 0.5)
            results.append(calculator.box_values(xyz, xyz, xyz))
        np.testing.assert_allclose(results[0], results[1], atol=1e-6)
        self.assertAlmostEqual(float(np.sum(results[1]) * (xyz[1]-xyz[0])**3), 0.5, places=4)
