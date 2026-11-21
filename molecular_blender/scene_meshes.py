# Copyright (C) 2014-2025 Shane Parker, Joshua Szekely
# SPDX-License-Identifier: GPL-3.0-or-later

"""Create and publish Blender meshes for computed fields."""
from itertools import chain

import bmesh
import bpy
import mathutils
import numpy as np

from .geometry_import import track_new
from .scene_materials import make_iso_materials
from .surface_calculation import iter_surface_frames
from .surface_geometry import field_geometry
from .util import stopwatch


def _write_indexed_triangles(mesh, vertices, faces):
    """Fill Blender's vertex, loop and polygon arrays without a second mesh or weld."""
    mesh.vertices.add(len(vertices))
    mesh.vertices.foreach_set('co', vertices.ravel())
    mesh.loops.add(faces.size)
    mesh.loops.foreach_set('vertex_index', faces.ravel())
    mesh.polygons.add(len(faces))
    mesh.polygons.foreach_set('loop_start', np.arange(0, faces.size, 3, dtype=np.int32))
    if hasattr(mesh, 'shade_flat'):
        mesh.shade_flat()
    mesh.update(calc_edges=True)


def create_surface_mesh_data(name, verts, faces, weld_distance):
    me = bpy.data.meshes.new(name)
    try:
        me['mb_generated'] = True
        if weld_distance is None:
            _write_indexed_triangles(me, verts, faces)
            return me
        me.from_pydata(verts, [], faces)
        me.update()
        bm = bmesh.new()
        try:
            bm.from_mesh(me)
            bmesh.ops.remove_doubles(bm, verts=bm.verts, dist=weld_distance)
            bm.to_mesh(me)
        finally:
            bm.free()
    except Exception:
        bpy.data.meshes.remove(me)
        raise
    return me


def create_mesh(name, verts, faces, material, context, remesh=8, weld_distance=0.01):
    mesh = create_surface_mesh_data(name, verts, faces, weld_distance)
    return _surface_object(name, mesh, material, remesh)


def _surface_object(name, me, material, remesh):
    ob = track_new(bpy.data.objects.new(name, me))
    if ob.data.materials:
        ob.data.materials[0] = material
    else:
        ob.data.materials.append(material)

    if remesh > 0:
        mod = ob.modifiers.new('Remesh', 'REMESH')
        mod.mode = 'SMOOTH'
        mod.octree_depth = remesh
        mod.scale = 0.99
        mod.use_smooth_shade = True
        mod.use_remove_disconnected = False

    return ob


def create_geometry(verts):
    """Build geometry out of vertices"""
    faces = []
    faceoffset = 0
    for ver in verts:
        if len(ver) == 4:
            faces.append((faceoffset + 0, faceoffset + 1,
                          faceoffset + 2, faceoffset + 3))
            faceoffset += 4
        elif len(ver) == 3:
            faces.append((faceoffset + 0, faceoffset + 1, faceoffset + 2))
            faceoffset += 3
    return list(chain.from_iterable(verts)), faces


class SurfaceMeshStage:
    """Build final Blender mesh data as frames arrive. Publish the series together."""
    def __init__(self, molecule, styler, options):
        self.molecule = molecule
        self.styler = styler
        self.options = options
        self.definitions = []
        self.frames = {}

    def add(self, surfaces, frames):
        self.frames.update(frames)
        for surface in surfaces:
            name = f"{self.molecule.name}_{surface['name']}_{surface['isovalue']}"
            vertices, faces = field_geometry(surface, surface.get('grid_spacing', self.options['resolution']))
            mesh = create_surface_mesh_data(name, vertices, faces, weld_distance=None)
            self.definitions.append(dict(name=surface['name'], series=surface.get('series', surface['name']),
                                         isovalue=surface['isovalue'], mesh_data=mesh))
            if 'ambiguous_phase_frame' in surface:
                self.definitions[-1]['ambiguous_phase_frame'] = surface['ambiguous_phase_frame']

    def finish(self, context):
        from .geometry_import import rollback_new_ids
        from .scene_state import RECORD
        owner = getattr(self.molecule, 'scene_owner', None)
        old_record = owner.get(RECORD) if owner else None
        try:
            with rollback_new_ids():
                create_surface_meshes(self.molecule, self.styler, context, self.options,
                                      self.definitions, self.frames)
        except Exception:
            if old_record is not None:
                owner[RECORD] = old_record
            raise
        self.definitions.clear()
        self.frames.clear()

    def discard(self):
        unused = [s['mesh_data'] for s in self.definitions
                  if s['mesh_data'].as_pointer() and not s['mesh_data'].users]
        if unused:
            bpy.data.batch_remove(unused)
        self.definitions.clear()
        self.frames.clear()


def create_surface_meshes(molecule, styler, context, options, vertex_sets, vertex_trajectory_map,
                          *, select_surfaces=True):
    """Create Blender mesh objects from precomputed isosurface data. Must run on main thread."""
    if not vertex_sets:
        return

    molecule_collection = getattr(molecule, "scene_collection", None) or bpy.data.collections[molecule.name]
    make_iso_materials(molecule, styler, vertex_sets, options)
    mol_obj = getattr(molecule, "scene_owner", None) or bpy.data.objects[molecule.name]

    # Vertices are in source coordinates. The imported empty starts at the
    # center of mass, while a wireframe root starts at the source origin.
    parent_inverse = mathutils.Matrix.Identity(4)
    if options["object_type"] != "wireframe":
        parent_inverse = mathutils.Matrix.Translation(-molecule.center_of_mass())
    groups = {}
    for surface in vertex_sets:
        groups.setdefault(surface['name'], []).append(surface)
    for func, fsets in groups.items():
        fmeshes = []
        for v in fsets:
            name = "{0}_{1}_{2}".format(molecule.name, v["name"], v["isovalue"])
            mesh = v.get('mesh_data')
            if mesh is None:
                verts, faces = field_geometry(v, v.get('grid_spacing', options['resolution']))
                mesh = create_surface_mesh_data(name, verts, faces, weld_distance=None)
            vmesh = _surface_object(name, mesh, v["material"], options["remesh"])
            v["mesh"] = vmesh
            fmeshes.append(vmesh)

        function_collection = track_new(bpy.data.collections.new(f"{molecule.name}-isosurface-{func}"))
        function_collection['mb_surface_collection'] = True
        molecule_collection.children.link(function_collection)
        for m in fmeshes:
            function_collection.objects.link(m)
            m.parent = mol_obj
            m.matrix_parent_inverse = parent_inverse
            if select_surfaces:
                m.select_set(True)

    if select_surfaces:
        mol_obj.select_set(True)
        context.view_layer.objects.active = mol_obj

    if hasattr(molecule, 'scene_owner'):
        from .scene_state import persist_surfaces
        persist_surfaces(mol_obj, vertex_sets, options)
        if options.get('surface_opacity', 1.0) != 1.0:
            from .scene_edit import _set_opacity
            for definition in vertex_sets:
                obj = definition['mesh']
                material = track_new(obj.active_material.copy())
                material['mb_style_material'] = True
                obj.material_slots[0].link = 'OBJECT'
                obj.material_slots[0].material = material
                _set_opacity(material, options['surface_opacity'])

    if vertex_trajectory_map:
        previous = set()
        for index, frame in enumerate(sorted(vertex_trajectory_map)):
            visible = set(vertex_trajectory_map[frame])
            changed = groups.keys() if index == 0 else previous ^ visible
            for name in changed:
                for surface in groups.get(name, ()):
                    obj = surface['mesh']
                    obj.hide_viewport = obj.hide_render = name not in visible
                    obj.keyframe_insert(data_path='hide_viewport', frame=frame)
                    obj.keyframe_insert(data_path='hide_render', frame=frame)
                    obj.animation_data.action['mb_generated'] = True
            previous = visible
        scene = context.scene
        scene.frame_set(scene.frame_current, subframe=scene.frame_subframe)

    if hasattr(molecule, 'scene_owner'):
        from .scene_surfaces import SOLO_GROUP, apply_solo, surface_groups
        if SOLO_GROUP in mol_obj:
            from .scene_state import read_record
            apply_solo(mol_obj, surface_groups(mol_obj, read_record(mol_obj)))


@stopwatch("draw surfaces")
def draw_surfaces(molecule, styler, context, options):
    """Draws isosurfaces"""
    stage = SurfaceMeshStage(molecule, styler, options)
    try:
        for batch in iter_surface_frames(molecule, options, indexed=True):
            stage.add(*batch)
            del batch
        stage.finish(context)
    finally:
        stage.discard()
