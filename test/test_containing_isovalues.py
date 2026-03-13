"""Tests for molecular_blender.containing_isovalues"""

import pytest
import numpy as np

from molecular_blender.containing_isovalues import (
    isovalue_containing_proportion,
    isovalue_containing_proportion_binary,
    isovalue_containing_proportion_sorted,
)


def _make_gaussian_grid(n=30, sigma=1.0):
    """Create a 3D Gaussian on a grid, returns (data, dV)."""
    x = np.linspace(-5, 5, n, dtype=np.float32)
    dx = x[1] - x[0]
    dV = dx ** 3
    xx, yy, zz = np.meshgrid(x, x, x, indexing="ij")
    data = np.exp(-(xx**2 + yy**2 + zz**2) / (2 * sigma**2)).astype(np.float32)
    norm = np.sqrt(np.sum(data**2) * dV)
    data /= norm
    return data, dV


class TestBinaryMethod:
    def _verify_containment(self, data, dV, isovalue, target, square=True):
        """Check that the isovalue actually contains ~target proportion."""
        mask = np.abs(data) >= abs(isovalue)
        if square:
            captured = np.sum(data[mask] ** 2) * dV
            total = np.sum(data ** 2) * dV
        else:
            captured = np.sum(np.abs(data[mask])) * dV
            total = np.sum(np.abs(data)) * dV
        ratio = captured / total
        assert abs(ratio - abs(target)) < 0.05, f"Expected ~{target}, got {ratio}"

    def test_50_percent(self):
        data, dV = _make_gaussian_grid()
        isovals = isovalue_containing_proportion_binary([0.5], data, dV)
        assert isovals[0] > 0
        self._verify_containment(data, dV, isovals[0], 0.5)

    def test_90_percent(self):
        data, dV = _make_gaussian_grid()
        iso50 = isovalue_containing_proportion_binary([0.5], data, dV)[0]
        iso90 = isovalue_containing_proportion_binary([0.9], data, dV)[0]
        assert iso90 < iso50
        self._verify_containment(data, dV, iso90, 0.9)

    def test_ordering(self):
        data, dV = _make_gaussian_grid()
        isovals = isovalue_containing_proportion_binary([0.3, 0.6, 0.9], data, dV)
        assert isovals[0] > isovals[1] > isovals[2]

    def test_non_square(self):
        data, dV = _make_gaussian_grid()
        isovals = isovalue_containing_proportion_binary([0.5], data, dV, square=False)
        assert isovals[0] > 0
        self._verify_containment(data, dV, isovals[0], 0.5, square=False)


class TestSortedMethod:
    def test_50_percent(self):
        data, dV = _make_gaussian_grid(n=20)
        isovals = isovalue_containing_proportion_sorted([0.5], data, dV)
        assert isovals[0] > 0
        mask = np.abs(data) >= abs(isovals[0])
        ratio = np.sum(data[mask] ** 2) * dV / (np.sum(data ** 2) * dV)
        assert abs(ratio - 0.5) < 0.1

    def test_ordering(self):
        data, dV = _make_gaussian_grid(n=20)
        isovals = isovalue_containing_proportion_sorted([0.3, 0.9], data, dV)
        assert isovals[0] > isovals[1]


class TestDispatch:
    def test_binary_dispatch_produces_valid_isovalue(self):
        data, dV = _make_gaussian_grid(n=20)
        result = isovalue_containing_proportion([0.5], data, dV, method='binary')
        assert 0 < result[0] <= np.max(np.abs(data))

    def test_sorted_dispatch_produces_valid_isovalue(self):
        data, dV = _make_gaussian_grid(n=20)
        result = isovalue_containing_proportion([0.5], data, dV, method='sorted')
        assert 0 < result[0] <= np.max(np.abs(data))

    def test_invalid_method_falls_back_to_binary(self):
        data, dV = _make_gaussian_grid(n=20)
        result_invalid = isovalue_containing_proportion([0.5], data, dV, method='invalid')
        result_binary = isovalue_containing_proportion([0.5], data, dV, method='binary')
        assert abs(result_invalid[0] - result_binary[0]) < 1e-10

    def test_both_methods_agree(self):
        data, dV = _make_gaussian_grid(n=20)
        binary = isovalue_containing_proportion([0.5], data, dV, method='binary')
        sorted_ = isovalue_containing_proportion([0.5], data, dV, method='sorted')
        assert abs(binary[0] - sorted_[0]) < 0.05
