"""Tests for compute_ao implementations (pure Python and Cython)"""

import pytest
import numpy as np

import molecular_blender.compute_ao_py as ao_py

try:
    import molecular_blender.compute_ao_cy as ao_cy
    HAS_CYTHON_AO = True
except ImportError:
    HAS_CYTHON_AO = False

from molecular_blender.constants import ang2bohr


def _implementations():
    """Yield (label, module) for each available implementation."""
    yield pytest.param(ao_py, id="python")
    if HAS_CYTHON_AO:
        yield pytest.param(ao_cy, id="cython")
    else:
        yield pytest.param(
            pytest.importorskip("molecular_blender.compute_ao_cy"),
            id="cython",
            marks=pytest.mark.skip(reason="Cython compute_ao not available"),
        )


@pytest.fixture(params=list(_implementations()))
def ao(request):
    return request.param


# -- gamma2n (python only) --

class TestGamma2n:
    def test_gamma2n_1(self):
        assert abs(ao_py.gamma2n_impl(1) - np.sqrt(np.pi)) < 1e-10

    def test_gamma2n_2(self):
        assert ao_py.gamma2n_impl(2) == 1

    def test_gamma2n_4(self):
        assert abs(ao_py.gamma2n_impl(4) - 1.0) < 1e-10

    def test_gamma2n_negative_raises(self):
        with pytest.raises(Exception):
            ao_py.gamma2n_impl(0)

    def test_gamma2n_array(self):
        result = ao_py.gamma2n(np.array([1.0, 2.0, 4.0]))
        assert abs(result[0] - np.sqrt(np.pi)) < 1e-10
        assert abs(result[1] - 1.0) < 1e-10


# -- gaussian_overlap (python only) --

class TestGaussianOverlap:
    def test_s_overlap(self):
        center = np.array([0.0, 0.0, 0.0])
        olap = ao_py.gaussian_overlap(center, 1.0, [0, 0, 0], center, 1.0, [0, 0, 0])
        expected = (np.pi / 2.0) ** 1.5
        assert abs(olap - expected) < 1e-5

    def test_p_overlap(self):
        center = np.array([0.0, 0.0, 0.0])
        olap = ao_py.gaussian_overlap(center, 1.0, [1, 0, 0], center, 1.0, [1, 0, 0])
        assert olap > 0

    def test_odd_angular_zero(self):
        center = np.array([0.0, 0.0, 0.0])
        olap = ao_py.gaussian_overlap(center, 1.0, [1, 0, 0], center, 1.0, [0, 0, 0])
        assert olap == 0.0


# -- polynomial functions (python only) --

class TestPolynomials:
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


# -- Shell class --

class TestShell:
    def _make_s_shell(self, ao):
        center = np.array([0.0, 0.0, 0.0], dtype=np.float32)
        exponents = np.array([1.0], dtype=np.float32)
        coeffs = np.array([1.0], dtype=np.float32)
        return ao.Shell(center, 0, exponents, coeffs, start=0)

    def _make_p_shell(self, ao):
        center = np.array([0.0, 0.0, 0.0], dtype=np.float32)
        exponents = np.array([1.0], dtype=np.float32)
        coeffs = np.array([1.0], dtype=np.float32)
        return ao.Shell(center, 1, exponents, coeffs, start=0)

    def test_s_shell_size(self, ao):
        sh = self._make_s_shell(ao)
        assert sh.size == 1

    def test_p_shell_size(self, ao):
        sh = self._make_p_shell(ao)
        assert sh.size == 3

    def test_d_shell_size(self, ao):
        center = np.array([0.0, 0.0, 0.0], dtype=np.float32)
        sh = ao.Shell(center, 2, np.array([1.0], dtype=np.float32),
                      np.array([1.0], dtype=np.float32), start=0)
        assert sh.size == 6

    def test_s_value_at_origin_py(self):
        sh = self._make_s_shell(ao_py)
        val = sh.value(0.0, 0.0, 0.0)
        assert val is not None
        assert abs(val[0]) > 0

    def test_s_value_far_away_py(self):
        sh = self._make_s_shell(ao_py)
        val = sh.value(100.0, 100.0, 100.0)
        assert val is None

    @pytest.mark.skipif(not HAS_CYTHON_AO, reason="Cython compute_ao not available")
    def test_s_value_at_origin_cy(self):
        sh = self._make_s_shell(ao_cy)
        coeff = np.array([1.0], dtype=np.float32)
        val = sh.value(0.0, 0.0, 0.0, coeff)
        assert val is not None

    @pytest.mark.skipif(not HAS_CYTHON_AO, reason="Cython compute_ao not available")
    def test_s_value_far_away_cy(self):
        """Cython returns 0.0 instead of None for screened-out shells."""
        sh = self._make_s_shell(ao_cy)
        coeff = np.array([1.0], dtype=np.float32)
        val = sh.value(100.0, 100.0, 100.0, coeff)
        assert val is None or val == 0.0

    def test_bounding_box_size_p(self, ao):
        """Use p shell since l=0 shells give bounding_box_size=0 by formula."""
        sh = self._make_p_shell(ao)
        bb = sh.bounding_box_size(1e-6)
        assert bb > 0


# -- Shell plane_values (python only, different API in cython) --

class TestShellPlaneValues:
    def test_s_plane(self):
        center = np.array([0.0, 0.0, 0.0], dtype=np.float32)
        sh = ao_py.Shell(center, 0, np.array([1.0], dtype=np.float32),
                         np.array([1.0], dtype=np.float32), start=0)
        xx = np.array([[0.0, 0.0], [1.0, 1.0]], dtype=np.float32)
        yy = np.array([[0.0, 1.0], [0.0, 1.0]], dtype=np.float32)
        result = sh.plane_values(xx, yy, 0.0)
        assert result is not None
        assert result.shape == (1, 2, 2)

    def test_far_plane_is_none(self):
        center = np.array([0.0, 0.0, 0.0], dtype=np.float32)
        sh = ao_py.Shell(center, 0, np.array([10.0], dtype=np.float32),
                         np.array([1.0], dtype=np.float32), start=0)
        xx = np.full((2, 2), 100.0, dtype=np.float32)
        yy = np.full((2, 2), 100.0, dtype=np.float32)
        result = sh.plane_values(xx, yy, 100.0)
        assert result is None


# -- OrbitalCalculater (python only; cython doesn't export MOData) --

class TestOrbitalCalculater:
    @pytest.fixture
    def water_orbital(self):
        from molecular_blender.importers import molecule_from_file
        geo = molecule_from_file(os.path.join(EXAMPLES, "water-sto-3g.molden"), {})
        modata = ao_py.MOData.from_dict(geo)
        return modata.get_orbital(2)

    def test_value_at_point(self, water_orbital):
        val = water_orbital.value(0.0, 0.0, 0.0)
        assert isinstance(val, float)

    def test_bounding_box(self, water_orbital):
        p0, p1 = water_orbital.bounding_box()
        assert all(p0[i] < p1[i] for i in range(3))

    def test_box_values_shape(self, water_orbital):
        xyz = np.linspace(-2, 2, 5, dtype=np.float32)
        box = water_orbital.box_values(xyz, xyz, xyz)
        assert box.shape == (5, 5, 5)

    def test_plane_values_shape(self, water_orbital):
        xy = np.linspace(-2, 2, 5, dtype=np.float32)
        plane = water_orbital.plane_values(xy, xy, 0.0)
        assert plane.shape == (5, 5)


# -- DensityCalculater (python only) --

class TestDensityCalculater:
    @pytest.fixture
    def fake_density(self):
        from molecular_blender.importers import molecule_from_file
        geo = molecule_from_file(os.path.join(EXAMPLES, "fake-atom.molden"), {})
        modata = ao_py.MOData.from_dict(geo)
        return modata.get_density()

    def test_value_at_point(self, fake_density):
        val = fake_density.value(0.0, 0.0, 0.0)
        assert val > 0

    def test_bounding_box(self, fake_density):
        p0, p1 = fake_density.bounding_box()
        assert all(p0[i] < p1[i] for i in range(3))

    def test_box_values_shape(self, fake_density):
        xyz = np.linspace(-5, 5, 8, dtype=np.float32)
        box = fake_density.box_values(xyz, xyz, xyz)
        assert box.shape == (8, 8, 8)

    def test_plane_values_shape(self, fake_density):
        xy = np.linspace(-5, 5, 8, dtype=np.float32)
        plane = fake_density.plane_values(xy, xy, 0.0)
        assert plane.shape == (8, 8)


import os
EXAMPLES = os.path.join(os.path.dirname(__file__), "examples")
