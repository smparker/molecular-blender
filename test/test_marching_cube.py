"""Tests for marching cube implementations (pure Python and Cython)"""

import pytest
import numpy as np

import molecular_blender.marching_cube_py as mc_py

try:
    import molecular_blender.marching_cube_cy as mc_cy
    HAS_CYTHON_MC = True
except ImportError:
    HAS_CYTHON_MC = False

DTYPE = np.float32


def _implementations():
    """Yield (label, module) for each available implementation."""
    yield pytest.param(mc_py, id="python")
    if HAS_CYTHON_MC:
        yield pytest.param(mc_cy, id="cython")
    else:
        yield pytest.param(
            pytest.importorskip("molecular_blender.marching_cube_cy"),
            id="cython",
            marks=pytest.mark.skip(reason="Cython marching cube not available"),
        )


@pytest.fixture(params=list(_implementations()))
def mc(request):
    return request.param


# -- vertexinterp (python only, cython may not expose it) --

class TestVertexInterp:
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

class TestActiveBox:
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

class TestPolygonise:
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


# -- marching_cube_box --

class TestMarchingCubeBox:
    def test_empty_grid(self, mc):
        data = np.zeros((3, 3, 3), dtype=DTYPE)
        result = mc.marching_cube_box(data, [0.5])
        assert len(result) == 1
        assert len(result[0]) == 0

    def test_single_peak(self, mc):
        data = np.zeros((5, 5, 5), dtype=DTYPE)
        data[2, 2, 2] = 1.0
        result = mc.marching_cube_box(data, [0.5])
        assert len(result) == 1
        assert len(result[0]) > 0

    def test_multiple_isovalues(self, mc):
        data = np.zeros((5, 5, 5), dtype=DTYPE)
        data[2, 2, 2] = 1.0
        result = mc.marching_cube_box(data, [0.3, 0.7])
        assert len(result) == 2
        assert len(result[0]) >= len(result[1])

    def test_negative_isovalue(self, mc):
        data = np.zeros((5, 5, 5), dtype=DTYPE)
        data[2, 2, 2] = -1.0
        result = mc.marching_cube_box(data, [-0.5])
        assert len(result[0]) > 0

    def test_full_grid(self, mc):
        data = np.ones((3, 3, 3), dtype=DTYPE)
        result = mc.marching_cube_box(data, [0.5])
        assert len(result[0]) == 0


# -- marching_cube_outline --

class TestMarchingCubeOutline:
    def test_empty_grid(self, mc):
        data = np.zeros((4, 4, 4), dtype=DTYPE)
        xvals = np.linspace(0, 1, 4, dtype=DTYPE)
        yvals = np.linspace(0, 1, 4, dtype=DTYPE)
        zvals = np.linspace(0, 1, 4, dtype=DTYPE)
        result = mc.marching_cube_outline(data, xvals, yvals, zvals, [0.5])
        assert len(result) == 0

    def test_with_peak(self, mc):
        data = np.zeros((4, 4, 4), dtype=DTYPE)
        data[2, 2, 2] = 1.0
        xvals = np.linspace(0, 1, 4, dtype=DTYPE)
        yvals = np.linspace(0, 1, 4, dtype=DTYPE)
        zvals = np.linspace(0, 1, 4, dtype=DTYPE)
        result = mc.marching_cube_outline(data, xvals, yvals, zvals, [0.5])
        assert len(result) > 0
        for p0, p1 in result:
            assert len(p0) == 3
            assert len(p1) == 3
