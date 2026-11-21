# -*- coding: utf-8 -*-
#
#  Molecular Blender
#  Filename: isosurfaces.py
#  Copyright (C) 2017-2025 Shane Parker, Joshua Szekely
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
# pylint: disable=bad-whitespace,multiple-statements,line-too-long,invalid-name

"""
Drivers for isosurface construction from quantum chemistry
information. For the meat of the isosurface construction,
see marching_cube.py
"""

from itertools import product
import numpy as np
import math

from .scalar_fields import CubeField, ScalarField
from .constants import ang2bohr, bohr2ang
from . import marching_cube

# Older binaries and the Python mesher expose the list-returning interface.
marching_cube_box = getattr(marching_cube, 'marching_cube_arrays', marching_cube.marching_cube_box)
marching_cube_outline = marching_cube.marching_cube_outline

DTYPE = np.float32

def transform_triangles(triangles, origin, axes):
    """Transforms and translates set of triangles to match input origin/axes"""
    if not len(triangles):
        return []
    coords = np.array(triangles, dtype=DTYPE)
    out = np.dot(coords, axes) + origin
    if np.linalg.det(axes) < 0:
        # A left-handed grid reflects the mesh. Preserve the field's normal.
        out = out[:, ::-1, :]
    return [ out[i,:] for i in range(out.shape[0]) ]


def cube_isosurface(data, origin, axes, isovalues, name="cube", wm=None):
    """Return set of triangles from cube file"""
    triangle_sets = [{"isovalue": iso, "name" : name} for iso in isovalues]

    tri_list = marching_cube_box(data, isovalues)

    for iiso in range(len(isovalues)):
        triangle_sets[iiso]["triangles"] = transform_triangles(
            tri_list[iiso], origin, axes)

    return triangle_sets


def molden_isosurface(orbital, isovalues, resolution, name="iso", wm=None, method="tiled", rough_resolution=0.6):
    """Return orbital triangles; resolution and rough_resolution are in angstrom."""
    p0, p1 = orbital.bounding_box(min([abs(x) for x in isovalues]) * 0.001)
    axes = np.eye(3, dtype=DTYPE) * bohr2ang

    return isosurface(p0, p1, resolution, isovalues, orbital.box_values, axes, name, wm,
                      method=method, rough_resolution=rough_resolution)


def field_isosurface(field: ScalarField, levels, resolution: float, name="iso", method="tiled"):
    """Contour a ScalarField on its angstrom bounds at the requested angstrom spacing."""
    if isinstance(field, CubeField):
        lattice = field.lattice
        return cube_isosurface(lattice['data'], lattice['origin'], lattice['axes'], levels, name)
    threshold = min((abs(x) for x in levels if x), default=0.01) * 0.001
    lo, hi = field.bounds(threshold)
    return isosurface(lo, hi, resolution, levels, field.grid, np.eye(3, dtype=DTYPE),
                      name, method=method, rough_resolution=0.6, units_per_angstrom=1.)


def isosurface(p0, p1, resolution, isovalues, box_func, axes, name, wm=None, method="adaptive",
        rough_resolution=0.3, max_subdivide=8, *, units_per_angstrom=ang2bohr):
    """Sample Bohr-coordinate bounds with an angstrom target grid spacing.

    axes maps sampled coordinates to the output geometry. The low-level simple
    and adaptive routines continue to use grid counts and Bohr spacing directly.
    """
    if method not in ('simple', 'adaptive', 'tiled'):
        raise ValueError("Isosurface method must be 'simple', 'adaptive' or 'tiled'")
    if not math.isfinite(resolution) or resolution <= 0:
        raise ValueError('Isosurface resolution must be finite and positive')
    if not math.isfinite(rough_resolution) or rough_resolution <= 0:
        raise ValueError('Rough resolution must be finite and positive')

    if method == "adaptive":
        return isosurface_adaptive(p0, p1, rough_resolution*units_per_angstrom, resolution*units_per_angstrom,
                                  max_subdivide, isovalues, box_func, axes, name, wm)
    npoints = [max(1, math.ceil((b-a)/(resolution*units_per_angstrom))) for a, b in zip(p0, p1)]
    if method == 'tiled':
        return isosurface_tiled(p0, p1, npoints, isovalues, box_func, axes, name)
    return isosurface_simple(p0, p1, npoints, isovalues, box_func, axes, name, wm)


def isosurface_tiled(p0, p1, counts, isovalues, box_func, axes, name):
    """Sample every grid point so even lobes inside a single coarse cell survive."""
    grid = [np.linspace(a, b, n + 1, dtype=DTYPE) for a, b, n in zip(p0, p1, counts)]
    step = (np.asarray(p1) - p0) / counts
    surfaces = [{"isovalue": iso, "name": name, "triangles": []} for iso in isovalues]
    # Share one layer of samples between adjacent tiles; every cell belongs to
    # exactly one tile. Shell screening and marching cubes keep their native kernels.
    width = 32
    for start in product(*(range(0, n, width) for n in counts)):
        values = box_func(*(axis[i:i+width+1] for axis, i in zip(grid, start)))
        low, high = values.min(), values.max()
        if not any(low <= level <= high for level in isovalues) or low == high:
            continue
        _mesh_tile(surfaces, values, isovalues, start, p0, step, axes)
    return surfaces


def _mesh_tile(surfaces, values, levels, start, origin, step, axes):
    for surface, triangles in zip(surfaces, marching_cube_box(values, levels)):
        if len(triangles):
            points = origin + (np.asarray(triangles) + start) * step
            surface["triangles"].extend((points @ axes).astype(DTYPE))


def _crossing_cells(values, levels):
    """Indices of grid cells with a sampled crossing at any requested level."""
    corners = [values[tuple(slice(i, None if i else -1) for i in corner)]
               for corner in product((0, 1), repeat=3)]
    low, high = corners[0].copy(), corners[0].copy()
    for corner in corners[1:]:
        np.minimum(low, corner, out=low)
        np.maximum(high, corner, out=high)
    active = np.zeros(low.shape, dtype=bool)
    for level in levels:
        active |= (low <= level) & (level <= high) & (low != high)
    return np.argwhere(active)


def isosurface_adaptive(p0, p1, resolution_start, resolution_end, max_subdivide, isovalues, box_func, axes, name, wm=None):
    """Refine coarse crossings, following each surface across sampled block faces."""
    if max_subdivide < 2:
        raise ValueError('Adaptive refinement requires at least two subdivisions')
    p0 = np.asarray(p0, dtype=float)
    extent = np.asarray(p1, dtype=float) - p0
    counts = np.maximum(1, np.ceil(extent / resolution_start).astype(int))
    target = np.maximum(1, np.ceil(extent / resolution_end).astype(int))
    grid = [np.linspace(a, b, n + 1, dtype=DTYPE) for a, b, n in zip(p0, p1, counts)]
    cells = _crossing_cells(box_func(*grid), isovalues)
    surfaces = [{"isovalue": iso, "name": name, "triangles": []} for iso in isovalues]
    while len(cells):
        divisions = np.clip((target + counts - 1) // counts, 1, max_subdivide)
        fine_counts = counts * divisions
        final = np.all(fine_counts >= target)
        # Every block uses slices of the same lattice, including its shared faces.
        grid = [np.linspace(a, b, n + 1, dtype=DTYPE) for a, b, n in zip(p0, p1, fine_counts)]
        step = extent / fine_counts
        pending = [tuple(cell) for cell in cells]
        visited = set(pending)
        children = []
        while pending:
            cell = pending.pop()
            start = divisions * cell
            samples = [axis[i:i+d+1] for axis, i, d in zip(grid, start, divisions)]
            values = box_func(*samples)
            # A lobe may cross a face twice between its coarse corners. Follow
            # fine crossings into every such neighbor, even beyond the outline.
            for axis in range(3):
                for side in (-1, 1):
                    index = cell[axis] + side
                    if not 0 <= index < counts[axis]:
                        continue
                    neighbor = cell[:axis] + (index,) + cell[axis+1:]
                    if neighbor in visited:
                        continue
                    face = values.take(0 if side < 0 else -1, axis=axis)
                    low, high = face.min(), face.max()
                    if low != high and any(low <= level <= high for level in isovalues):
                        pending.append(neighbor)
                        visited.add(neighbor)
            if final:
                _mesh_tile(surfaces, values, isovalues, start, p0, step, axes)
            else:
                children.extend(start + _crossing_cells(values, isovalues))
        if final:
            break
        counts = fine_counts
        cells = children
    return surfaces


def isosurface_outline(p0, p1, nvoxels, isovalues, box_func, axes, name, wm=None):
    """Returns low resolution outline of an isosurface"""
    xvals, xstep = np.linspace(p0[0], p1[0], num=nvoxels[0]+1, retstep=True, endpoint=True, dtype=DTYPE)
    yvals, ystep = np.linspace(p0[1], p1[1], num=nvoxels[1]+1, retstep=True, endpoint=True, dtype=DTYPE)
    zvals, zstep = np.linspace(p0[2], p1[2], num=nvoxels[2]+1, retstep=True, endpoint=True, dtype=DTYPE)

    box_values = box_func(xvals, yvals, zvals)
    return marching_cube_outline(box_values, xvals, yvals, zvals, isovalues)


def isosurface_simple(p0, p1, nvoxels, isovalues, box_func, axes, name, wm=None):
    """Return set of triangles from function object in single pass algorithm"""
    triangle_sets = [{"isovalue": iso, "name" : name} for iso in isovalues]

    xvals, xstep = np.linspace(p0[0], p1[0], num=nvoxels[0]+1, retstep=True, endpoint=True, dtype=DTYPE)
    yvals, ystep = np.linspace(p0[1], p1[1], num=nvoxels[1]+1, retstep=True, endpoint=True, dtype=DTYPE)
    zvals, zstep = np.linspace(p0[2], p1[2], num=nvoxels[2]+1, retstep=True, endpoint=True, dtype=DTYPE)

    scaled_axes = np.array(axes, dtype=DTYPE) * np.array([xstep, ystep, zstep], dtype=DTYPE)[:, None]

    box_values = box_func(xvals, yvals, zvals)
    tri_list = marching_cube_box(box_values, isovalues)

    origin = np.array(p0, dtype=DTYPE) @ axes
    for iiso in range(len(isovalues)):
        triangle_sets[iiso]["triangles"] = transform_triangles(
            tri_list[iiso], origin, scaled_axes)

    return triangle_sets
