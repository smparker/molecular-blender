"""Pack trajectory triangles into the blend file and decode the displayed sample."""
from bisect import bisect_right
import struct
import zlib

import bpy
from bpy.app.handlers import persistent
import numpy as np


ARCHIVE = 'mb_packed_contours'
TIMES = 'mb_contour_times'
INDEX = 'mb_contour_sample'


def pack_sample(vertices, faces):
    if not len(faces):
        return b''
    header = struct.pack('<II', len(vertices), len(faces))
    return zlib.compress(header + np.asarray(vertices, dtype='<f4').tobytes()
                         + np.asarray(faces, dtype='<i4').tobytes(), level=1)


def show_sample(mesh, frame):
    from .scene_meshes import fill_indexed_triangles
    times = mesh[TIMES]
    index = max(0, bisect_right(times, frame) - 1)
    if mesh.get(INDEX) == index or mesh.is_editmode:
        return
    sample = mesh[ARCHIVE].get(str(index))
    if sample:
        raw = zlib.decompress(sample)
        nv, nf = struct.unpack_from('<II', raw)
        if len(raw) != 8 + 12 * (nv + nf):
            raise ValueError('Invalid packed molecular trajectory mesh')
        vertices = np.frombuffer(raw, dtype='<f4', count=3*nv, offset=8).reshape(nv, 3)
        faces = np.frombuffer(raw, dtype='<i4', count=3*nf, offset=8+12*nv).reshape(nf, 3)
        if faces.min() < 0 or faces.max() >= nv:
            raise ValueError('Invalid face in packed molecular trajectory')
    mesh.clear_geometry()
    if sample:
        fill_indexed_triangles(mesh, vertices, faces)
    else:
        mesh.update()
    mesh[INDEX] = index


@persistent
def update_contours(scene, depsgraph=None):
    if not scene.get(ARCHIVE):
        return
    frame = scene.frame_current + scene.frame_subframe
    for obj in scene.objects:
        if obj.type == 'MESH' and TIMES in obj.data:
            show_sample(obj.data, frame)


def register():
    if update_contours not in bpy.app.handlers.frame_change_pre:
        bpy.app.handlers.frame_change_pre.append(update_contours)


def unregister():
    if update_contours in bpy.app.handlers.frame_change_pre:
        bpy.app.handlers.frame_change_pre.remove(update_contours)
