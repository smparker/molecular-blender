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
from .marching_cube import marching_cube_box, edgetable

DTYPE = np.float32

def active_box(cornervalues, isolevels):
    """Determines whether any isovalues are crossed based on given corner values"""
    for isolevel in isolevels:
        cubeindex = 0

        if cornervalues[0] < isolevel: cubeindex |= 1
        if cornervalues[1] < isolevel: cubeindex |= 2
        if cornervalues[2] < isolevel: cubeindex |= 4
        if cornervalues[3] < isolevel: cubeindex |= 8
        if cornervalues[4] < isolevel: cubeindex |= 16
        if cornervalues[5] < isolevel: cubeindex |= 32
        if cornervalues[6] < isolevel: cubeindex |= 64
        if cornervalues[7] < isolevel: cubeindex |= 128

        if edgetable[cubeindex] != 0:
            return True

    return False

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


def molden_isosurface(orbital, isovalues, resolution, name="iso", wm=None, method="adaptive", rough_resolution=0.2*ang2bohr):
    """Return set of triangles from Molden file"""
    p0, p1 = orbital.bounding_box(min([abs(x) for x in isovalues]) * 0.001)
    axes = np.eye(3, dtype=DTYPE) * bohr2ang

    return isosurface(p0, p1, resolution, isovalues, orbital.box_values, axes, name, wm, rough_resolution=rough_resolution)


def isosurface(p0, p1, resolution, isovalues, box_func, axes, name, wm=None, method="adaptive",
        rough_resolution=0.2*ang2bohr, max_subdivide=4):
    """Return set of triangles from function object"""

    if method == "adaptive":
        return isosurface_adaptive(p0, p1, rough_resolution, resolution, max_subdivide, isovalues, box_func, axes, name, wm)
    else:
        npoints = [ int(math.ceil((b-a)/(resolution*ang2bohr))) for a, b in zip(p0, p1) ]
        return isosurface_simple(p0, p1, npoints, isovalues, box_func, axes, name, wm)


def isosurface_adaptive(p0, p1, resolution_start, resolution_end, max_subdivide, isovalues, box_func, axes, name, wm=None):
    """Return set of triangles from function object from multi-pass adaptive algorithm"""

    # let's start by finding the nearest integer that provides a rough layer of at least 0.25 A
    np_start = [ int(math.ceil((b - a)/resolution_start)) for a, b in zip(p0, p1) ]
    np_min_end = [ int(math.ceil((b - a)/resolution_end)) for a, b in zip(p0, p1) ]

    subdivide = [ max_subdivide for x in range(3) ]

    outlines = isosurface_outline(p0, p1, np_start, isovalues, box_func, axes, name)

    # np_now is the effective number of points that would be produced by the final pass
    np_now = [ a * s for a, s in zip(np_start, subdivide) ]

    while not all([now > want for now, want in zip(np_now, np_min_end)]):
        new_outlines = []
        boxgen = ( (pp0, pp1, subdivide, isovalues, box_func, axes, name, wm) for pp0, pp1 in outlines )
        #with multiprocessing.Pool(processes=1) as pool:
        new_outlines = starmap(isosurface_outline, boxgen)

        outlines = sum(new_outlines, [])
        #for pp0, pp1 in outlines:
        #    o = isosurface_outline(pp0, pp1, subdivide, isovalues, box_func, axes, name, wm)
        #    new_outlines.extend(o)

        #outlines = new_outlines

        np_now = [ a * s for a, s in zip(np_now, subdivide) ]

    # now build actual triangles from the list of active boxes from above
    triangles = None
    vertlist = None
    boxgen = ( (pp0, pp1, subdivide, isovalues, box_func, axes, name) for pp0, pp1 in outlines )
    #with multiprocessing.Pool() as pool:
    vertlist = starmap(isosurface_simple, boxgen)

    for verts in vertlist:
        if triangles is None:
            triangles = verts
        else:
            for i in range(len(verts)):
                triangles[i]["triangles"].extend(verts[i]["triangles"])

    return triangles

def isosurface_outline(p0, p1, npoints, isovalues, box_func, axes, name, wm=None):
    """Returns low resolution outline of an isosurface"""
    xvals, xstep = np.linspace(p0[0], p1[0], num=npoints[0], retstep=True, endpoint=True, dtype=DTYPE)
    yvals, ystep = np.linspace(p0[1], p1[1], num=npoints[1], retstep=True, endpoint=True, dtype=DTYPE)
    zvals, zstep = np.linspace(p0[2], p1[2], num=npoints[2], retstep=True, endpoint=True, dtype=DTYPE)
    nx, ny, nz = npoints
    r = (xstep, ystep, zstep)

    outlines = [ ]
    z = p0[2]
    box_values = box_func(xvals, yvals, zvals)

    plane_values_1 = box_values[:,:,0]

    cornervalues = [0] * 8

    if wm is not None:
        wm.progress_begin(0, nz)

    for zi in range(1, nz):
        z = zvals[zi-1]
        z2 = zvals[zi]
        plane_values_2 = box_values[:,:,zi]
        for yi in range(ny-1):
            y = yvals[yi]
            y2 = yvals[yi+1]
            for xi in range(nx-1):
                x = xvals[xi]
                x2 = xvals[xi+1]
                cornervalues = [
                    plane_values_1[xi][yi],
                    plane_values_1[xi][yi + 1],
                    plane_values_1[xi + 1][yi + 1],
                    plane_values_1[xi + 1][yi],
                    plane_values_2[xi][yi],
                    plane_values_2[xi][yi + 1],
                    plane_values_2[xi + 1][yi + 1],
                    plane_values_2[xi + 1][yi],
                ]

                if active_box(cornervalues, isovalues):
                    outlines.append( ( (x, y, z), (x2, y2, z2) ) )

        plane_values_1 = plane_values_2

        if wm is not None:
            wm.progress_update(zi)

    if wm is not None:
        wm.progress_end()

    return outlines

def isosurface_simple(p0, p1, npoints, isovalues, box_func, axes, name, wm=None):
    """Return set of triangles from function object in single pass algorithm"""
    triangle_sets = [{"isovalue": iso, "name" : name} for iso in isovalues]

    xvals, xstep = np.linspace(p0[0], p1[0], num=npoints[0], retstep=True, endpoint=True, dtype=DTYPE)
    yvals, ystep = np.linspace(p0[1], p1[1], num=npoints[1], retstep=True, endpoint=True, dtype=DTYPE)
    zvals, zstep = np.linspace(p0[2], p1[2], num=npoints[2], retstep=True, endpoint=True, dtype=DTYPE)

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
