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

import multiprocessing
from itertools import starmap
import numpy as np
import math

from .constants import ang2bohr, bohr2ang
from .marching_cube import marching_cube_box, marching_cube_outline

DTYPE = np.float32

def transform_triangles(triangles, origin, axes):
    """Transforms and translates set of triangles to match input origin/axes"""
    if not triangles:
        return triangles
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


def molden_isosurface(orbital, isovalues, resolution, name="iso", wm=None, method="adaptive", rough_resolution=0.6):
    """Return orbital triangles; resolution and rough_resolution are in angstrom."""
    p0, p1 = orbital.bounding_box(min([abs(x) for x in isovalues]) * 0.001)
    axes = np.eye(3, dtype=DTYPE) * bohr2ang

    return isosurface(p0, p1, resolution, isovalues, orbital.box_values, axes, name, wm,
                      method=method, rough_resolution=rough_resolution)


def isosurface(p0, p1, resolution, isovalues, box_func, axes, name, wm=None, method="adaptive",
        rough_resolution=0.3, max_subdivide=8):
    """Sample Bohr-coordinate bounds with an angstrom target grid spacing.

    axes maps sampled coordinates to the output geometry. The low-level simple
    and adaptive routines continue to use grid counts and Bohr spacing directly.
    """
    if method not in ('simple', 'adaptive'):
        raise ValueError("Isosurface method must be 'simple' or 'adaptive'")
    if not math.isfinite(resolution) or resolution <= 0:
        raise ValueError('Isosurface resolution must be finite and positive')
    if not math.isfinite(rough_resolution) or rough_resolution <= 0:
        raise ValueError('Rough resolution must be finite and positive')

    if method == "adaptive":
        return isosurface_adaptive(p0, p1, rough_resolution*ang2bohr, resolution*ang2bohr,
                                  max_subdivide, isovalues, box_func, axes, name, wm)
    else:
        npoints = [ int(math.ceil((b-a)/(resolution*ang2bohr))) for a, b in zip(p0, p1) ]
        return isosurface_simple(p0, p1, npoints, isovalues, box_func, axes, name, wm)


def isosurface_adaptive(p0, p1, resolution_start, resolution_end, max_subdivide, isovalues, box_func, axes, name, wm=None):
    """Return set of triangles from function object from multi-pass adaptive algorithm"""

    np_start = [ int(math.ceil((b - a)/resolution_start))+1 for a, b in zip(p0, p1) ]
    np_min_end = [ int(math.ceil((b - a)/resolution_end))+1 for a, b in zip(p0, p1) ]

    subdivide = [ min(max_subdivide, max(2, int(math.ceil(want / start)))) for start, want in zip(np_start, np_min_end) ]

    outlines = isosurface_outline(p0, p1, np_start, isovalues, box_func, axes, name)

    np_now = [ a * s for a, s in zip(np_start, subdivide) ]

    while not all([now > want for now, want in zip(np_now, np_min_end)]):
        boxgen = ( (pp0, pp1, subdivide, isovalues, box_func, axes, name, wm) for pp0, pp1 in outlines )
        new_outlines = starmap(isosurface_outline, boxgen)

        outlines = sum(new_outlines, [])

        np_now = [ a * s for a, s in zip(np_now, subdivide) ]

    if not outlines:
        return [{"isovalue": iso, "name" : name, "triangles": []} for iso in isovalues]

    # optimal fine_subdivide: only subdivide as much as needed to reach target resolution
    fine_subdivide = [ min(max_subdivide, max(2, int(math.ceil(want * s / now))))
                       for want, now, s in zip(np_min_end, np_now, subdivide) ]

    # pre-compute grid parameters (same for all boxes in this pass)
    pp0_ref, pp1_ref = outlines[0]
    xstep = DTYPE((pp1_ref[0] - pp0_ref[0]) / fine_subdivide[0])
    ystep = DTYPE((pp1_ref[1] - pp0_ref[1]) / fine_subdivide[1])
    zstep = DTYPE((pp1_ref[2] - pp0_ref[2]) / fine_subdivide[2])

    scaled_axes = np.array(axes, dtype=DTYPE) * np.array([xstep, ystep, zstep], dtype=DTYPE)

    dx = np.arange(fine_subdivide[0] + 1, dtype=DTYPE) * xstep
    dy = np.arange(fine_subdivide[1] + 1, dtype=DTYPE) * ystep
    dz = np.arange(fine_subdivide[2] + 1, dtype=DTYPE) * zstep

    niso = len(isovalues)
    raw_tris = [[] for _ in range(niso)]
    tri_box_info = [[] for _ in range(niso)]

    for pp0, pp1 in outlines:
        xvals = DTYPE(pp0[0]) + dx
        yvals = DTYPE(pp0[1]) + dy
        zvals = DTYPE(pp0[2]) + dz

        box_data = box_func(xvals, yvals, zvals)
        tri_list = marching_cube_box(box_data, isovalues)

        origin_x = DTYPE(pp0[0]) * bohr2ang
        origin_y = DTYPE(pp0[1]) * bohr2ang
        origin_z = DTYPE(pp0[2]) * bohr2ang

        for iiso in range(niso):
            tris = tri_list[iiso]
            if tris:
                start = len(raw_tris[iiso])
                raw_tris[iiso].extend(tris)
                tri_box_info[iiso].append((start, len(tris), origin_x, origin_y, origin_z))

    triangle_sets = [{"isovalue": iso, "name" : name, "triangles": []} for iso in isovalues]

    for iiso in range(niso):
        if raw_tris[iiso]:
            coords = np.array(raw_tris[iiso], dtype=DTYPE)
            transformed = np.dot(coords, scaled_axes)
            for start, count, ox, oy, oz in tri_box_info[iiso]:
                transformed[start:start+count, :, 0] += ox
                transformed[start:start+count, :, 1] += oy
                transformed[start:start+count, :, 2] += oz
            triangle_sets[iiso]["triangles"] = [ transformed[i] for i in range(transformed.shape[0]) ]

    return triangle_sets

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

    scaled_axes = np.array(axes, dtype=DTYPE) * np.array([xstep, ystep, zstep], dtype=DTYPE)

    box_values = box_func(xvals, yvals, zvals)
    tri_list = marching_cube_box(box_values, isovalues)

    origin = np.array(p0, dtype=DTYPE) * bohr2ang
    for iiso in range(len(isovalues)):
        triangle_sets[iiso]["triangles"] = transform_triangles(
            tri_list[iiso], origin, scaled_axes)

    return triangle_sets
