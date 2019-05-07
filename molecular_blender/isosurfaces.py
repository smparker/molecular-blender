# -*- coding: utf-8 -*-
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
    out = np.dot(coords, axes)
    return [out[i, :] + origin for i in range(coords.shape[0])]


def cube_isosurface(data, origin, axes, isovalues, name="cube", wm=None):
    """Return set of triangles from cube file"""
    triangle_sets = [{"isovalue": iso, "name" : name} for iso in isovalues]

    tri_list = marching_cube_box(data, isovalues)

    for iiso in range(len(isovalues)):
        triangle_sets[iiso]["triangles"] = transform_triangles(
            tri_list[iiso], origin, axes)

    return triangle_sets


def molden_isosurface(orbital, isovalues, resolution, name="iso", wm=None, method="adaptive", rough_resolution=0.3):
    """Return set of triangles from Molden file"""
    p0, p1 = orbital.bounding_box(min([abs(x) for x in isovalues]) * 0.001)
    axes = np.eye(3, dtype=DTYPE) * bohr2ang

    return isosurface(p0, p1, resolution, isovalues, orbital.box_values, axes, name, wm, rough_resolution=rough_resolution)


def isosurface(p0, p1, resolution, isovalues, box_func, axes, name, wm=None, method="adaptive",
        rough_resolution=0.3, max_subdivide=8):
    """Return set of triangles from function object"""

    if method == "adaptive":
        return isosurface_adaptive(p0, p1, rough_resolution, resolution, max_subdivide, isovalues, box_func, axes, name, wm)
    else:
        npoints = [ int(math.ceil((b-a)/(resolution*ang2bohr))) for a, b in zip(p0, p1) ]
        return isosurface_simple(p0, p1, npoints, isovalues, box_func, axes, name, wm)


def isosurface_adaptive(p0, p1, resolution_start, resolution_end, max_subdivide, isovalues, box_func, axes, name, wm=None):
    """Return set of triangles from function object from multi-pass adaptive algorithm"""

    # let's start by finding the nearest integer that provides a rough layer of at least 0.25 A
    np_start = [ int(math.ceil((b - a)/resolution_start))+1 for a, b in zip(p0, p1) ]
    np_min_end = [ int(math.ceil((b - a)/resolution_end))+1 for a, b in zip(p0, p1) ]

    subdivide = [ max_subdivide for x in range(3) ]

    outlines = isosurface_outline(p0, p1, np_start, isovalues, box_func, axes, name)

    # np_now is the effective number of points that would be produced by the final pass
    np_now = [ a * s for a, s in zip(np_start, subdivide) ]

    while not all([now > want for now, want in zip(np_now, np_min_end)]):
        new_outlines = []
        boxgen = ( (pp0, pp1, subdivide, isovalues, box_func, axes, name, wm) for pp0, pp1 in outlines )
        new_outlines = starmap(isosurface_outline, boxgen)

        outlines = sum(new_outlines, [])

        np_now = [ a * s for a, s in zip(np_now, subdivide) ]

    # now build actual triangles from the list of active boxes from above
    triangles = None
    vertlist = None
    boxgen = ( (pp0, pp1, subdivide, isovalues, box_func, axes, name) for pp0, pp1 in outlines )
    vertlist = starmap(isosurface_simple, boxgen)

    for verts in vertlist:
        if triangles is None:
            triangles = verts
        else:
            for i in range(len(verts)):
                triangles[i]["triangles"].extend(verts[i]["triangles"])

    return triangles

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

    scaled_axes = np.array(axes, dtype=DTYPE)
    scaled_axes[:,0] *= xstep
    scaled_axes[:,1] *= ystep
    scaled_axes[:,2] *= zstep

    box_values = box_func(xvals, yvals, zvals)
    tri_list = marching_cube_box(box_values, isovalues)

    for iiso in range(len(isovalues)):
        origin = np.array(p0, dtype=DTYPE) * bohr2ang
        triangle_sets[iiso]["triangles"] = transform_triangles(
            tri_list[iiso], origin, scaled_axes)

    return triangle_sets
