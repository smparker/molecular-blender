"""Scalar fields sampled in source-space angstroms. Values retain their native units."""
from __future__ import annotations

from dataclasses import dataclass
from itertools import product
from typing import Callable, Optional, TYPE_CHECKING

import numpy as np
if TYPE_CHECKING:
    from numpy.typing import ArrayLike, NDArray

from .constants import ang2bohr, bohr2ang
from .esp import sample_cube


@dataclass(frozen=True)
class FieldInfo:
    source: str
    selection: str
    units: Optional[str]
    signed: Optional[bool]


class ScalarField:
    """Point, batch (N,3) and Cartesian-product grid (nx,ny,nz) sampling.

    Grid axes are monotonic. Bounds are a pair of xyz corners in angstroms. An optional gradient takes
    (N,3) points and returns (N,3) derivatives in field units per angstrom.
    """
    info: FieldInfo
    gradient: Optional[Callable[[ArrayLike], NDArray]] = None

    def point(self, xyz: ArrayLike) -> float:
        return float(self.batch(np.asarray(xyz).reshape(1, 3))[0])

    def batch(self, points: ArrayLike) -> NDArray:
        raise NotImplementedError

    def grid(self, x: ArrayLike, y: ArrayLike, z: ArrayLike) -> NDArray:
        raise NotImplementedError

    def bounds(self, threshold: float = 1e-5) -> tuple[NDArray, NDArray]:
        raise NotImplementedError


def _points(points):
    xyz = np.asarray(points, dtype=float)
    if xyz.ndim != 2 or xyz.shape[1] != 3 or not np.isfinite(xyz).all():
        raise ValueError("Expected finite xyz points with shape (N,3)")
    return xyz


def _axes(x, y, z):
    axes = tuple(np.asarray(axis) for axis in (x, y, z))
    if any(axis.ndim != 1 or not np.isfinite(axis).all() for axis in axes):
        raise ValueError("Expected three finite one-dimensional Cartesian grid axes")
    if any(len(a) > 1 and not (np.all(a[1:] >= a[:-1]) or np.all(a[1:] <= a[:-1])) for a in axes):
        raise ValueError("Cartesian grid axes must be monotonic")
    return axes


class CubeField(ScalarField):
    """A selected dataset on a regular lattice. Cube values have file-defined units."""
    def __init__(self, grid, *, source='', field_type='auto'):
        self.lattice = grid
        index = grid.get('dataset_index', 0)
        label = grid['labels'][index] if 'labels' in grid else index
        kind = grid.get('field_type', 'unknown') if field_type == 'auto' else field_type
        if kind not in ('density', 'orbital', 'signed', 'unknown'):
            raise ValueError("Cube field type must be auto, orbital, density, or signed")
        self.info = FieldInfo(source or grid.get('source', ''), f'dataset {label}', grid.get('units'),
                              {'density': False, 'orbital': True, 'signed': True}.get(kind))

    def batch(self, points):
        return sample_cube(self.lattice, _points(points))

    def grid(self, x, y, z):
        axes = _axes(x, y, z)
        shape = tuple(len(axis) for axis in axes)
        result = np.empty(shape, dtype=np.float32)
        for start in range(0, result.size, 65536):
            chunk = slice(start, min(start + 65536, result.size))
            indices = np.unravel_index(np.arange(chunk.start, chunk.stop), shape)
            result.flat[chunk] = self.batch(np.column_stack([a[i] for a, i in zip(axes, indices)]))
        return result

    def bounds(self, threshold=1e-5):
        grid = self.lattice
        corners = np.array(list(product(*((0, n - 1) for n in grid['nres']))))
        points = grid['origin'] + corners @ grid['axes']
        return points.min(axis=0), points.max(axis=0)


class WavefunctionField(ScalarField):
    """Adapt an AO evaluator in Bohr to the common angstrom coordinate system."""
    def __init__(self, evaluator, *, source='', selection, orbital=False):
        self.evaluator = evaluator
        units = 'bohr^-3/2' if orbital else 'e/bohr^3'
        if selection.startswith(('response ', 'polar')):
            units = None  # Response vector normalization is producer-specific.
        self.info = FieldInfo(source, selection, units, orbital or selection not in ('density', 'alpha', 'beta'))

    def batch(self, points):
        return np.asarray([self.evaluator.value(*p) for p in _points(points) * ang2bohr], dtype=np.float32)

    def grid(self, x, y, z):
        axes = [np.ascontiguousarray(a * ang2bohr, dtype=np.float32) for a in _axes(x, y, z)]
        if any(axis.size == 0 for axis in axes):
            return np.empty(tuple(a.size for a in axes), dtype=np.float32)
        return self.evaluator.box_values(*axes)

    def bounds(self, threshold=1e-5):
        if not np.isfinite(threshold) or threshold <= 0:
            raise ValueError("Bounding threshold must be finite and positive")
        lo, hi = self.evaluator.bounding_box(threshold)
        return np.asarray(lo) * bohr2ang, np.asarray(hi) * bohr2ang
