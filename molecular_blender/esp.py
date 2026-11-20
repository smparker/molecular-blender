"""Trilinear sampling of cube fields in source coordinates (angstrom)."""
from itertools import product

import numpy as np


def sample_cube(grid, positions):
    """Trilinear interpolation on the cube's row-vector axes, including boundary faces."""
    points = np.asarray(positions, dtype=np.float64).reshape(-1, 3)
    coordinates = (points - grid['origin']) @ np.linalg.inv(grid['axes'].astype(np.float64))
    return _interpolate(grid, coordinates).astype(np.float32)


def _interpolate(grid, coordinates):
    upper = np.array(grid['nres']) - 1
    if not np.isfinite(coordinates).all() or ((coordinates < -1e-5) | (coordinates > upper + 1e-5)).any():
        raise ValueError('Source cube does not cover all requested points')
    coordinates = np.clip(coordinates, 0, upper)
    low = np.floor(coordinates).astype(np.intp)
    high = np.minimum(low + 1, upper)
    fraction = coordinates - low
    values = np.zeros(len(coordinates))
    for corner in product((0, 1), repeat=3):
        indices = tuple((high if bit else low)[:, axis] for axis, bit in enumerate(corner))
        values += grid['data'][indices] * np.prod(np.where(corner, fraction, 1 - fraction), axis=1)
    return values


def add_resampled_cube(total, target, source, weight):
    """Accumulate on the first term's grid in bounded batches, without extrapolation."""
    inverse = np.linalg.inv(source['axes'].astype(np.float64))
    offset = (target['origin'].astype(np.float64) - source['origin']) @ inverse
    transform = target['axes'] @ inverse
    flat = total.reshape(-1)
    for start in range(0, flat.size, 65536):
        end = min(start + 65536, flat.size)
        indices = np.column_stack(np.unravel_index(np.arange(start, end), total.shape))
        flat[start:end] += weight * _interpolate(source, offset + indices @ transform)
