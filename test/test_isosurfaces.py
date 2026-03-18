"""Tests for molecular_blender.isosurfaces"""

import os
import unittest
import numpy as np

from molecular_blender.isosurfaces import (
    transform_triangles, cube_isosurface, isosurface_simple,
    isosurface_adaptive, isosurface, molden_isosurface,
)
from molecular_blender.constants import bohr2ang

EXAMPLES = os.path.join(os.path.dirname(__file__), "../examples")
DTYPE = np.float32


class TestTransformTriangles(unittest.TestCase):
    def test_empty(self):
        assert transform_triangles([], np.zeros(3), np.eye(3)) == []

    def test_identity(self):
        tri = [[0, 0, 0], [1, 0, 0], [0, 1, 0]]
        origin = np.array([0, 0, 0], dtype=DTYPE)
        axes = np.eye(3, dtype=DTYPE)
        result = transform_triangles(tri, origin, axes)
        assert len(result) == 3
        np.testing.assert_allclose(result[0], [0, 0, 0], atol=1e-5)

    def test_translation(self):
        tri = [[0, 0, 0], [1, 0, 0], [0, 1, 0]]
        origin = np.array([10, 20, 30], dtype=DTYPE)
        axes = np.eye(3, dtype=DTYPE)
        result = transform_triangles(tri, origin, axes)
        np.testing.assert_allclose(result[0], [10, 20, 30], atol=1e-5)

    def test_scaling(self):
        tri = [[1, 0, 0]]
        origin = np.zeros(3, dtype=DTYPE)
        axes = np.eye(3, dtype=DTYPE) * 2.0
        result = transform_triangles(tri, origin, axes)
        np.testing.assert_allclose(result[0], [2, 0, 0], atol=1e-5)


class TestCubeIsosurface(unittest.TestCase):
    def test_basic(self):
        data = np.zeros((5, 5, 5), dtype=DTYPE)
        data[2, 2, 2] = 1.0
        origin = np.zeros(3, dtype=DTYPE)
        axes = np.eye(3, dtype=DTYPE)
        result = cube_isosurface(data, origin, axes, [0.5])
        assert len(result) == 1
        assert "triangles" in result[0]
        assert "isovalue" in result[0]
        assert len(result[0]["triangles"]) > 0

    def test_no_crossing(self):
        data = np.zeros((3, 3, 3), dtype=DTYPE)
        origin = np.zeros(3, dtype=DTYPE)
        axes = np.eye(3, dtype=DTYPE)
        result = cube_isosurface(data, origin, axes, [0.5])
        assert len(result[0]["triangles"]) == 0


class TestIsosurfaceSimple(unittest.TestCase):
    def _gaussian_box(self, xvals, yvals, zvals):
        """Simple Gaussian blob for testing."""
        out = np.empty((len(xvals), len(yvals), len(zvals)), dtype=DTYPE)
        for i, x in enumerate(xvals):
            for j, y in enumerate(yvals):
                for k, z in enumerate(zvals):
                    r2 = x*x + y*y + z*z
                    out[i, j, k] = np.exp(-r2)
        return out

    def test_basic(self):
        p0 = np.array([-3.0, -3.0, -3.0])
        p1 = np.array([3.0, 3.0, 3.0])
        npoints = [10, 10, 10]
        axes = np.eye(3, dtype=DTYPE) * bohr2ang
        result = isosurface_simple(
            p0, p1, npoints, [0.5], self._gaussian_box, axes, "test"
        )
        assert len(result) == 1
        assert len(result[0]["triangles"]) > 0


class TestIsosurfaceAdaptive(unittest.TestCase):
    def _gaussian_box(self, xvals, yvals, zvals):
        out = np.empty((len(xvals), len(yvals), len(zvals)), dtype=DTYPE)
        for i, x in enumerate(xvals):
            for j, y in enumerate(yvals):
                for k, z in enumerate(zvals):
                    r2 = x*x + y*y + z*z
                    out[i, j, k] = np.exp(-r2)
        return out

    def test_basic(self):
        p0 = np.array([-3.0, -3.0, -3.0])
        p1 = np.array([3.0, 3.0, 3.0])
        axes = np.eye(3, dtype=DTYPE) * bohr2ang
        result = isosurface_adaptive(
            p0, p1, 1.0, 0.5, 4, [0.5], self._gaussian_box, axes, "test"
        )
        assert len(result) == 1
        assert len(result[0]["triangles"]) > 0

    def test_no_crossing(self):
        def zero_box(xvals, yvals, zvals):
            return np.zeros((len(xvals), len(yvals), len(zvals)), dtype=DTYPE)

        p0 = np.array([-1.0, -1.0, -1.0])
        p1 = np.array([1.0, 1.0, 1.0])
        axes = np.eye(3, dtype=DTYPE) * bohr2ang
        result = isosurface_adaptive(
            p0, p1, 0.5, 0.3, 4, [0.5], zero_box, axes, "test"
        )
        assert len(result[0]["triangles"]) == 0


class TestIsosurfaceDispatch(unittest.TestCase):
    def _gaussian_box(self, xvals, yvals, zvals):
        out = np.empty((len(xvals), len(yvals), len(zvals)), dtype=DTYPE)
        for i, x in enumerate(xvals):
            for j, y in enumerate(yvals):
                for k, z in enumerate(zvals):
                    r2 = x*x + y*y + z*z
                    out[i, j, k] = np.exp(-r2)
        return out

    def test_simple_method(self):
        p0 = np.array([-2.0, -2.0, -2.0])
        p1 = np.array([2.0, 2.0, 2.0])
        axes = np.eye(3, dtype=DTYPE) * bohr2ang
        result = isosurface(
            p0, p1, 0.5, [0.5], self._gaussian_box, axes, "test",
            method="simple"
        )
        assert len(result) == 1
        assert result[0]["isovalue"] == 0.5
        assert len(result[0]["triangles"]) > 0
        for tri in result[0]["triangles"]:
            assert len(tri) == 3

    def test_adaptive_method(self):
        p0 = np.array([-2.0, -2.0, -2.0])
        p1 = np.array([2.0, 2.0, 2.0])
        axes = np.eye(3, dtype=DTYPE) * bohr2ang
        result = isosurface(
            p0, p1, 0.5, [0.5], self._gaussian_box, axes, "test",
            method="adaptive"
        )
        assert len(result) == 1
        assert result[0]["isovalue"] == 0.5
        assert len(result[0]["triangles"]) > 0
        for tri in result[0]["triangles"]:
            assert len(tri) == 3

    def test_simple_and_adaptive_same_topology(self):
        """Both methods should produce a closed surface with same topology for a simple blob."""
        p0 = np.array([-2.0, -2.0, -2.0])
        p1 = np.array([2.0, 2.0, 2.0])
        axes = np.eye(3, dtype=DTYPE) * bohr2ang
        simple = isosurface(
            p0, p1, 0.5, [0.5], self._gaussian_box, axes, "test",
            method="simple"
        )
        adaptive = isosurface(
            p0, p1, 0.5, [0.5], self._gaussian_box, axes, "test",
            method="adaptive"
        )
        assert len(simple[0]["triangles"]) > 0
        assert len(adaptive[0]["triangles"]) > 0


class TestMoldenIsosurface(unittest.TestCase):
    def test_small_orbital(self):
        from molecular_blender.importers import molecule_from_file
        from molecular_blender.orbitals import MOData
        geo = molecule_from_file(os.path.join(EXAMPLES, "fake-atom.molden"), {})
        modata = MOData.from_dict(geo)
        orb = modata.get_orbital(1)
        result = molden_isosurface(orb, [0.05], resolution=0.3, name="test_orb")
        assert len(result) == 1
        assert len(result[0]["triangles"]) > 0
