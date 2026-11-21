# Copyright (C) 2014-2025 Shane Parker, Joshua Szekely
# SPDX-License-Identifier: GPL-3.0-or-later

"""Index marching-cubes triangles in source-space angstroms."""
from itertools import product
from math import dist, floor
import numpy as np

from . import _HAS_BPY


def _weld_representatives(vertices, tolerance):
    """Pick the first unmerged neighbor in stable vertex order."""
    representatives = np.arange(len(vertices), dtype=np.int32)
    if _HAS_BPY:
        from mathutils.kdtree import KDTree
        tree = KDTree(len(vertices))
        for index, point in enumerate(vertices):
            tree.insert(point, index)
        tree.balance()
        for index, point in enumerate(vertices):
            if representatives[index] == index:
                for _, other, _ in tree.find_range(point, tolerance):
                    if other > index and representatives[other] == other:
                        representatives[other] = index
        return representatives

    cells = {}
    offsets = tuple(product((-1, 0, 1), repeat=3))
    for index, point in enumerate(vertices):
        x, y, z = (floor(float(value) / tolerance) for value in point)
        matches = (other for dx, dy, dz in offsets
                   for other in cells.get((x+dx, y+dy, z+dz), ())
                   if dist(vertices[other], point) <= tolerance)
        representative = min(matches, default=index)
        representatives[index] = representative
        if representative == index:
            cells.setdefault((x, y, z), []).append(index)
    return representatives


def field_weld_distance(resolution):
    # Intersections are shared to float32 precision. A fixed 0.01 Å weld
    # also merges distinct vertices on fine grids and can tear a closed mesh.
    return resolution * 1e-4


def index_field_triangles(triangles, resolution):
    """Share float32 intersections and close sub-grid seams before allocating a Blender mesh."""
    points = np.ascontiguousarray(triangles, dtype=np.float32).reshape(-1, 3)
    if not len(points):
        return points, np.empty((0, 3), dtype=np.int32)
    # Sort complete XYZ triples as one key. Most shared edges match exactly.
    keys = points.view(np.dtype((np.void, points.strides[0]))).ravel()
    _, first, corners = np.unique(keys, return_index=True, return_inverse=True)
    vertices = points[first]
    faces = corners.reshape(-1, 3)

    # Adjacent marching cells can interpolate the same edge in opposite directions.
    # Merge residual roundoff by distance, including points across rounding-bin boundaries.
    representatives = _weld_representatives(vertices, field_weld_distance(resolution))
    faces = representatives[faces]
    a, b, c = faces.T
    faces = faces[(a != b) & (b != c) & (c != a)]
    # Collapsing a tiny sliver can also make two faces share all three corners.
    face_keys = np.sort(faces, axis=1).view(np.dtype((np.void, 3 * faces.itemsize))).ravel()
    _, first_face = np.unique(face_keys, return_index=True)
    faces = faces[np.sort(first_face)]
    used, indices = np.unique(faces, return_inverse=True)
    return vertices[used], indices.reshape(-1, 3).astype(np.int32)


def field_geometry(definition, resolution):
    if 'faces' in definition:
        return definition['vertices'], definition['faces']
    return index_field_triangles(definition['triangles'], resolution)
