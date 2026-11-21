# -*- coding: utf-8 -*-
#
#  Molecular Blender
#  Filename: containing_isovalues.py
#  Copyright (C) 2025 Shane Parker
#
#  This file is part of Molecular Blender.
#
#  Molecular Blender is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 3, or (at your option)
#  any later version.
#
#  Molecular Blender is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU Library General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with Molecular Blender; see COPYING.
#  If not, see <http://www.gnu.org/licenses/>.
#

from __future__ import annotations

from typing import TYPE_CHECKING
import numpy as np

if TYPE_CHECKING:
    from numpy.typing import ArrayLike


def _prepare(values, field, dV, square):
    fractions = np.asarray(values, dtype=np.float64)
    amplitudes = np.abs(np.asarray(field, dtype=np.float64)).ravel()
    if fractions.ndim != 1 or not np.all(np.isfinite(fractions)) or np.any(np.abs(fractions) > 1):
        raise ValueError('Signed containment fractions must be finite and within [-1, 1]')
    if not np.isfinite(dV) or dV <= 0:
        raise ValueError('Voxel volume must be finite and positive')
    if not amplitudes.size or not np.all(np.isfinite(amplitudes)):
        raise ValueError('Containment requires a nonempty finite field')
    maximum = amplitudes.max()
    if maximum == 0:
        raise ValueError('Containment is undefined for an all-zero field')
    # Common volume and amplitude factors cancel in a fraction. Scale first to
    # avoid overflow/underflow when squaring large or tiny field values.
    weights = amplitudes / maximum
    if square:
        weights = weights * weights
    weights /= weights.sum(dtype=np.float64)
    return fractions, amplitudes, weights


def isovalue_containing_proportion(values: ArrayLike, volume_data_sh: ArrayLike, dV: float,
                                  square: bool = True, method: str = 'binary') -> list[float]:
    """Return signed thresholds for fractions of a sampled field's integral.

    square=True integrates psi squared; False integrates absolute field values
    (density or signed-density magnitude). Fraction signs affect output signs
    only. Zero containment returns a threshold above the largest magnitude;
    full containment includes every nonzero sample. Discrete grids and tied
    magnitudes may prevent an exact target. Thresholds are not rounded, since
    rounding can change which samples are included. Unknown methods retain the
    historical binary fallback.
    """
    if method == 'sorted':
        return isovalue_containing_proportion_sorted(values, volume_data_sh, dV, square)
    return isovalue_containing_proportion_binary(values, volume_data_sh, dV, square)


def isovalue_containing_proportion_sorted(values, volume_data_sh, dV, square=True):
    fractions, amplitudes, weights = _prepare(values, volume_data_sh, dV, square)
    order = np.argsort(-amplitudes, kind='stable')
    levels = amplitudes[order]
    cumulative = np.cumsum(weights[order], dtype=np.float64)
    # Only whole groups of equal magnitudes can be included by a threshold.
    ends = np.r_[np.flatnonzero(levels[:-1] != levels[1:]), len(levels)-1]
    thresholds = np.r_[np.nextafter(levels[0], np.inf), levels[ends]]
    captured = np.r_[0., cumulative[ends]]
    out = []
    for fraction in fractions:
        if abs(fraction) == 1:
            threshold = amplitudes[amplitudes > 0].min()
        else:
            errors = np.abs(captured - abs(fraction))
            # Include the larger mass when equally close.
            threshold = thresholds[np.flatnonzero(errors == errors.min())[-1]]
        out.append(float(np.copysign(threshold, fraction)))
    return out


def isovalue_containing_proportion_binary(values, volume_data_sh, dV, square=True,
                                         tolerance=1e-3, max_iterations=100):
    if not np.isfinite(tolerance) or tolerance <= 0:
        raise ValueError('Tolerance must be finite and positive')
    if isinstance(max_iterations, bool) or not isinstance(max_iterations, (int, np.integer)) or max_iterations < 1:
        raise ValueError('max_iterations must be a positive integer')
    fractions, amplitudes, weights = _prepare(values, volume_data_sh, dV, square)
    minimum = amplitudes[amplitudes > 0].min()
    maximum = amplitudes.max()
    outside = np.nextafter(maximum, np.inf)
    out = []
    for fraction in fractions:
        target = abs(fraction)
        if target == 0:
            out.append(float(np.copysign(outside, fraction)))
            continue
        if target == 1:
            out.append(float(np.copysign(minimum, fraction)))
            continue
        low, high = minimum, outside
        best, best_mass = (minimum, 1.) if target >= .5 else (outside, 0.)
        best_error = abs(best_mass - target)
        for _ in range(max_iterations):
            threshold = low / 2 + high / 2
            if threshold <= low or threshold >= high:
                break
            mass = weights[amplitudes >= threshold].sum(dtype=np.float64)
            error = abs(mass - target)
            if error < best_error or (error == best_error and mass > best_mass):
                best, best_mass, best_error = threshold, mass, error
            if best_error <= tolerance:
                break
            if mass > target:
                low = threshold
            else:
                high = threshold
        out.append(float(np.copysign(best, fraction)))
    return out
