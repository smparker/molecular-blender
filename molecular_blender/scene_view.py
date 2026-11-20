"""Frame displayed molecular geometry in the invoking 3D viewport."""
import bpy
import blf
import numpy as np
from bpy_extras.view3d_utils import location_3d_to_region_2d
from mathutils import Matrix, Vector

from .scene_state import ID, ROLE, owned_objects, read_record, resolve_owner


def atom_index_labels(context, owner):
    legacy = None
    for obj in owned_objects(owner, 'atom'):
        if not obj.visible_get(view_layer=context.view_layer, viewport=context.space_data):
            continue
        index = obj.get('mb_atom_index')
        if index is None:
            if legacy is None:
                legacy = {atom['id']: atom['index'] for atom in read_record(owner)['atoms']}
            index = legacy.get(obj.get(ID))
        if index is not None:
            yield obj.matrix_world.translation, str(index)


def draw_atom_indices():
    context = bpy.context
    owner = resolve_owner(context.object)
    if owner is None or not owner.mb_show_atom_indices or not context.space_data.overlay.show_overlays:
        return
    size = round(13 * context.preferences.system.ui_scale)
    if bpy.app.version < (4, 0, 0):
        blf.size(0, size, 72)
    else:
        blf.size(0, size)
    blf.color(0, 1, 1, 1, 1)
    blf.enable(0, blf.SHADOW)
    blf.shadow(0, 3, 0, 0, 0, 1)
    try:
        for position, label in atom_index_labels(context, owner):
            xy = location_3d_to_region_2d(context.region, context.region_data, position)
            if xy is not None and 0 <= xy.x < context.region.width and 0 <= xy.y < context.region.height:
                blf.position(0, xy.x + 4, xy.y + 4, 0)
                blf.draw(0, label)
    finally:
        blf.disable(0, blf.SHADOW)


def principal_view_rotation(positions, current):
    """Put the longest geometric axis across the screen and the shortest toward the viewer."""
    if len(positions) < 2:
        return current
    centered = np.array(positions, dtype=np.float64)
    centered -= centered.mean(axis=0)
    spread, axes = np.linalg.eigh(centered.T @ centered)
    tolerance = .01 * spread[-1]
    if spread[-1] - spread[0] <= tolerance:
        return current

    # Resolve signs and near-degenerate eigenspaces using the invoking viewport.
    right, up, front = np.asarray(current.to_matrix()).T
    def closest_in(subspace, candidates):
        for candidate in candidates:
            vector = subspace @ (subspace.T @ candidate)
            length = np.linalg.norm(vector)
            if length > 1e-6:
                return vector / length

    # A nearly spherical cloud can have two small adjacent gaps. Choose the
    # better-defined axis and treat the perpendicular plane as degenerate.
    short_gap, long_gap = np.diff(spread)
    if short_gap < long_gap and short_gap <= tolerance:
        normal_space = axes[:, :2]
        long_space = axes[:, 2:]
    else:
        normal_space = axes[:, :1]
        long_space = axes[:, 1:] if long_gap <= tolerance else axes[:, 2:]
    x = closest_in(long_space, (right, up, front))
    z = closest_in(normal_space, (front, up, right))
    return Matrix(np.column_stack((x, np.cross(z, x), z))).to_quaternion()


class _FrameMolecule:
    face_on = False
    principal_axes = False

    @classmethod
    def poll(cls, context):
        return context.area and context.area.type == 'VIEW_3D' and resolve_owner(context.object) is not None

    def execute(self, context):
        owner = resolve_owner(context.object)
        space = context.space_data
        region = context.region
        if region.type != 'WINDOW':
            region = next(r for r in context.area.regions if r.type == 'WINDOW')
        with context.temp_override(region=region):
            view = context.region_data
            rotation = view.view_rotation.copy()
            if view.view_perspective == 'CAMERA':
                rotation = view.view_matrix.to_quaternion().inverted()
            if self.face_on:
                atoms = context.selected_objects
                if len(atoms) != 3 or any(obj.get(ROLE) != 'atom' or resolve_owner(obj) != owner for obj in atoms):
                    self.report({'INFO'}, 'Select three atoms from the same molecule')
                    return {'CANCELLED'}
                a, b, c = (obj.matrix_world.translation for obj in atoms)
                ab, ac = b - a, c - a
                normal = ab.normalized().cross(ac.normalized())
                if normal.length < 1e-5:
                    self.report({'INFO'}, 'The three selected atoms must define a plane')
                    return {'CANCELLED'}
                normal.normalize()
                # Use the closest side and keep the current screen orientation.
                toward_viewer = rotation @ Vector((0, 0, 1))
                if normal.dot(toward_viewer) < 0:
                    normal.negate()
                rotation = toward_viewer.rotation_difference(normal) @ rotation
            if self.principal_axes:
                if owner.type == 'CURVE':
                    # Wireframe imports share junctions across several curve splines.
                    curve = owner.evaluated_get(context.evaluated_depsgraph_get())
                    positions = list({tuple(curve.matrix_world @ point.co.to_3d())
                                      for spline in curve.data.splines
                                      for point in (*spline.points, *spline.bezier_points)})
                else:
                    positions = [o.matrix_world.translation for o in owned_objects(owner, 'atom')
                                 if o.visible_get(view_layer=context.view_layer, viewport=space)]
                if not positions:
                    self.report({'INFO'}, 'Show an atom to align the molecule')
                    return {'CANCELLED'}
                rotation = principal_view_rotation(positions, rotation)
            to_view = rotation.inverted().to_matrix().to_4x4()
            depsgraph = context.evaluated_depsgraph_get()
            low = [float('inf')] * 3
            high = [-float('inf')] * 3
            for obj in [owner, *owned_objects(owner)]:
                if obj.type not in {'MESH', 'CURVE', 'SURFACE', 'META', 'FONT'} or not obj.visible_get(
                        view_layer=context.view_layer, viewport=space):
                    continue
                evaluated = obj.evaluated_get(depsgraph)
                bounds = evaluated.bound_box
                if all(tuple(corner) == (-1., -1., -1.) for corner in bounds):
                    continue
                transform = to_view @ evaluated.matrix_world
                for corner in bounds:
                    point = transform @ Vector(corner)
                    for axis in range(3):
                        low[axis] = min(low[axis], point[axis])
                        high[axis] = max(high[axis], point[axis])
            if low[0] > high[0]:
                self.report({'INFO'}, 'This molecule has no visible geometry in this viewport')
                return {'CANCELLED'}

            # Leave camera view before changing any view parameters, including
            # when Lock Camera to View is on. Use the current lens and aspect.
            if view.view_perspective == 'CAMERA':
                view.view_perspective = 'PERSP'
            view.view_rotation = rotation
            view.update()
            half = (Vector(high) - Vector(low)) * 0.5
            projection = view.window_matrix
            extent = max(half.x * abs(projection[0][0]), half.y * abs(projection[1][1]))
            if view.is_perspective:
                distance = half.z + max(1.15 * extent, 2 * space.clip_start)
            else:
                distance = max(1.15 * extent * view.view_distance, 2 * space.clip_start)
            view.view_location = rotation @ ((Vector(high) + Vector(low)) * 0.5)
            view.view_distance = distance
            view.update()
        context.area.tag_redraw()
        return {'FINISHED'}


class MB_OT_frame_molecule(_FrameMolecule, bpy.types.Operator):
    bl_idname = 'molecular_blender.frame_molecule'
    bl_label = 'Frame Molecule'
    bl_description = 'Fit the visible molecule in this viewport, preserving the selection and render camera'


class MB_OT_view_face_on(_FrameMolecule, bpy.types.Operator):
    bl_idname = 'molecular_blender.view_face_on'
    bl_label = 'View Face-On'
    bl_description = 'Look perpendicular to the plane of three selected atoms and frame their molecule'
    face_on = True


class MB_OT_view_principal_axes(_FrameMolecule, bpy.types.Operator):
    bl_idname = 'molecular_blender.view_principal_axes'
    bl_label = 'View Principal Axes'
    bl_description = 'Frame the molecule face-on to its best-fit plane, with its longest dimension horizontal'
    principal_axes = True
