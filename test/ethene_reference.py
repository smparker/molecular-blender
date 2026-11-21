"""Reference positions for the planar and twisted ethene pi orbital."""
from pathlib import Path

import numpy as np

from molecular_blender.constants import bohr2ang


FIXTURES = Path(__file__).resolve().parent / 'fixtures' / 'pyscf_ethene'
LEVELS = (.05, -.05, .005, -.005)


def outer_crossing_errors(triangles, crossings, level):
    """Compare the outermost mesh hit on each vertical ray to PySCF's continuous root."""
    triangles = np.asarray(triangles, dtype=np.float64)
    a, b, c = triangles.transpose(1, 0, 2)
    ab, ac = b - a, c - a
    det = ab[:, 0] * ac[:, 1] - ab[:, 1] * ac[:, 0]
    errors = []
    for x, y, signed_level, z in crossings:
        if signed_level != level:
            continue
        dx, dy = x * bohr2ang - a[:, 0], y * bohr2ang - a[:, 1]
        valid = det != 0
        u = np.divide(dx * ac[:, 1] - dy * ac[:, 0], det, where=valid, out=np.zeros_like(det))
        v = np.divide(ab[:, 0] * dy - ab[:, 1] * dx, det, where=valid, out=np.zeros_like(det))
        hit = valid & (u >= -1e-9) & (v >= -1e-9) & (u + v <= 1 + 1e-9)
        heights = (a[:, 2] + u * ab[:, 2] + v * ac[:, 2])[hit]
        if not len(heights):
            raise AssertionError(f'No contour crossing at x={x}, y={y}, level={level}')
        outer = max(heights) if z > 0 else min(heights)
        errors.append(abs(outer - z * bohr2ang))
    if not errors:
        raise AssertionError(f'No reference crossings for level {level}')
    return np.array(errors)
