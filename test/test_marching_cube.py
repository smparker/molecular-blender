"""Tests for marching cube implementations (pure Python and Cython)"""

import unittest
import numpy as np

import molecular_blender.marching_cube_py as mc_py

try:
    import molecular_blender.marching_cube_cy as mc_cy
    HAS_CYTHON_MC = True
except ImportError:
    HAS_CYTHON_MC = False

DTYPE = np.float32


class TestSurfaceCoordinates(unittest.TestCase):
    def test_linear_field_crossings_in_grid_coordinates(self):
        for backend in ([mc_py, mc_cy] if HAS_CYTHON_MC else [mc_py]):
            for axis in range(3):
                with self.subTest(backend=backend.__name__, axis=axis):
                    data = np.indices((3, 3, 3), dtype=DTYPE)[axis].copy()
                    vertices = np.asarray(backend.marching_cube_box(data, [0.25])[0])
                    self.assertGreater(vertices.size, 0)
                    np.testing.assert_allclose(vertices[:, :, axis], 0.25, atol=1e-6)
                    self.assertGreaterEqual(vertices.min(), 0)
                    self.assertLessEqual(vertices.max(), 2)


# -- vertexinterp (python only, cython may not expose it) --

class TestVertexInterp(unittest.TestCase):
    def test_at_p1(self):
        result = mc_py.vertexinterp(0.0, [0, 0, 0], [1, 1, 1], 0.0, 1.0)
        assert result == [0, 0, 0]

    def test_at_p2(self):
        result = mc_py.vertexinterp(1.0, [0, 0, 0], [1, 1, 1], 0.0, 1.0)
        assert result == [1, 1, 1]

    def test_equal_values(self):
        result = mc_py.vertexinterp(0.5, [0, 0, 0], [1, 1, 1], 0.5, 0.5)
        assert result == [0, 0, 0]

    def test_midpoint(self):
        x, y, z = mc_py.vertexinterp(0.5, [0, 0, 0], [2, 2, 2], 0.0, 1.0)
        assert abs(x - 1.0) < 1e-5
        assert abs(y - 1.0) < 1e-5
        assert abs(z - 1.0) < 1e-5


# -- active_box (python only) --

class TestActiveBox(unittest.TestCase):
    def test_all_inside(self):
        corners = [0.0] * 8
        assert mc_py.active_box(corners, [0.5]) is False

    def test_all_above_isovalue(self):
        """All corners above isovalue means no surface crossing."""
        corners = [1.0] * 8
        assert mc_py.active_box(corners, [0.5]) is False

    def test_crossing_detected(self):
        corners = [0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0]
        assert mc_py.active_box(corners, [0.5]) is True

    def test_mixed(self):
        corners = [0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0]
        assert mc_py.active_box(corners, [0.5]) is True


# -- polygonise (python only) --

class TestPolygonise(unittest.TestCase):
    def test_no_crossing(self):
        corners = [0.0] * 8
        tris = mc_py.polygonise(corners, 0.5, 0, 0, 0, 1, 1, 1)
        assert tris == []

    def test_single_corner_inside(self):
        corners = [1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        tris = mc_py.polygonise(corners, 0.5, 0, 0, 0, 1, 1, 1)
        assert len(tris) > 0
        for tri in tris:
            assert len(tri) == 3


# -- marching_cube_box (parametrized over implementations) --

class _TestMarchingCubeBoxMixin:
    def test_empty_grid(self):
        data = np.zeros((3, 3, 3), dtype=DTYPE)
        result = self.mc.marching_cube_box(data, [0.5])
        assert len(result) == 1
        assert len(result[0]) == 0

    def test_single_peak(self):
        data = np.zeros((5, 5, 5), dtype=DTYPE)
        data[2, 2, 2] = 1.0
        result = self.mc.marching_cube_box(data, [0.5])
        assert len(result) == 1
        assert len(result[0]) > 0

    def test_multiple_isovalues(self):
        data = np.zeros((5, 5, 5), dtype=DTYPE)
        data[2, 2, 2] = 1.0
        result = self.mc.marching_cube_box(data, [0.3, 0.7])
        assert len(result) == 2
        assert len(result[0]) >= len(result[1])

    def test_negative_isovalue(self):
        data = np.zeros((5, 5, 5), dtype=DTYPE)
        data[2, 2, 2] = -1.0
        result = self.mc.marching_cube_box(data, [-0.5])
        assert len(result[0]) > 0

    def test_full_grid(self):
        data = np.ones((3, 3, 3), dtype=DTYPE)
        result = self.mc.marching_cube_box(data, [0.5])
        assert len(result[0]) == 0


class TestMarchingCubeBoxPy(_TestMarchingCubeBoxMixin, unittest.TestCase):
    def setUp(self):
        self.mc = mc_py


@unittest.skipUnless(HAS_CYTHON_MC, "Cython marching cube not available")
class TestMarchingCubeBoxCy(_TestMarchingCubeBoxMixin, unittest.TestCase):
    def setUp(self):
        self.mc = mc_cy


# -- marching_cube_outline (parametrized over implementations) --

class _TestMarchingCubeOutlineMixin:
    def test_empty_grid(self):
        data = np.zeros((4, 4, 4), dtype=DTYPE)
        xvals = np.linspace(0, 1, 4, dtype=DTYPE)
        yvals = np.linspace(0, 1, 4, dtype=DTYPE)
        zvals = np.linspace(0, 1, 4, dtype=DTYPE)
        result = self.mc.marching_cube_outline(data, xvals, yvals, zvals, [0.5])
        assert len(result) == 0

    def test_with_peak(self):
        data = np.zeros((4, 4, 4), dtype=DTYPE)
        data[2, 2, 2] = 1.0
        xvals = np.linspace(0, 1, 4, dtype=DTYPE)
        yvals = np.linspace(0, 1, 4, dtype=DTYPE)
        zvals = np.linspace(0, 1, 4, dtype=DTYPE)
        result = self.mc.marching_cube_outline(data, xvals, yvals, zvals, [0.5])
        assert len(result) > 0
        for p0, p1 in result:
            assert len(p0) == 3
            assert len(p1) == 3


class TestMarchingCubeOutlinePy(_TestMarchingCubeOutlineMixin, unittest.TestCase):
    def setUp(self):
        self.mc = mc_py


@unittest.skipUnless(HAS_CYTHON_MC, "Cython marching cube not available")
class TestMarchingCubeOutlineCy(_TestMarchingCubeOutlineMixin, unittest.TestCase):
    def setUp(self):
        self.mc = mc_cy


@unittest.skipUnless(HAS_CYTHON_MC, 'requires compiled mesher')
class TestNativeGridBuffers(unittest.TestCase):
    def test_array_and_list_interfaces_match(self):
        grid = np.indices((3, 4, 5), dtype=DTYPE)[0].copy()
        levels = [0.25, 8.]
        arrays = mc_cy.marching_cube_arrays(grid, levels)
        self.assertIsInstance(mc_cy.marching_cube_box(grid, levels)[0], list)
        for array, legacy in zip(arrays, mc_cy.marching_cube_box(grid, levels)):
            self.assertEqual(array.dtype, np.float32)
            self.assertTrue(array.flags.c_contiguous)
            np.testing.assert_array_equal(array, np.asarray(legacy).reshape(-1, 3, 3))
        self.assertEqual(arrays[1].shape, (0, 3, 3))

    def test_outline_checks_axis_sizes_and_empty_levels(self):
        grid = np.zeros((2, 2, 2), dtype=DTYPE)
        axis = np.arange(2, dtype=DTYPE)
        with self.assertRaisesRegex(ValueError, 'dimensions'):
            mc_cy.marching_cube_outline(grid, axis, axis[:1], axis, [0.])
        self.assertEqual(mc_cy.marching_cube_outline(grid, axis, axis, axis, []), [])
