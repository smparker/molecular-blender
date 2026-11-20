"""Frame displayed molecular geometry in the invoking 3D viewport."""
import bpy
from mathutils import Vector

from .scene_state import owned_objects, resolve_owner


class MB_OT_frame_molecule(bpy.types.Operator):
    bl_idname = 'molecular_blender.frame_molecule'
    bl_label = 'Frame Molecule'
    bl_description = 'Fit the visible molecule in this viewport, preserving the selection and render camera'

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
