"""Release generated Blender IDs when the last molecule using them is removed."""
import bpy

from .scene_state import MOLECULE_ID, RECORD, ROLE, normalize_duplicate_ids, owned_objects, resolve_owner


RESOURCE_OWNER = 'mb_resource_owner'
GENERATED = 'mb_generated'
PROFILE = 'molecular_blender_bond_profile'


def mark_import_resources(owner, blocks):
    for block in blocks:
        if block.as_pointer():
            block[RESOURCE_OWNER] = owner[MOLECULE_ID]
            block[GENERATED] = True
            animation = getattr(block, 'animation_data', None)
            if animation and animation.action:
                animation.action[RESOURCE_OWNER] = owner[MOLECULE_ID]
                animation.action[GENERATED] = True


def discard_profile(profile):
    """Bevel references can survive even when Blender reports only a fake user."""
    if (profile is None or profile.library or not profile.get(PROFILE)
            or profile.users != int(profile.use_fake_user) or profile.users_collection
            or any(curve.bevel_object == profile for curve in bpy.data.curves)):
        return
    data = profile.data
    bpy.data.objects.remove(profile)
    if data and not data.users:
        bpy.data.curves.remove(data)


def remove_molecule(owner):
    """Remove tagged geometry and unused generated dependencies in dependency order."""
    if owner is None or RECORD not in owner or owner.library:
        raise ValueError('Select a local molecule to remove')
    normalize_duplicate_ids()
    identifier = owner[MOLECULE_ID]
    objects = {obj for obj in owned_objects(owner)
               if obj.get(GENERATED) or obj.get(ROLE) in ('atom', 'bond', 'surface')}
    objects.add(owner)
    if any(obj.library for obj in objects):
        raise ValueError('Make linked molecule objects local before removing them')

    # Include dependencies shared with a previous import and generated blocks
    # that never acquired a user (e.g. a bond material in a stick model).
    geometry = {obj.data for obj in objects if obj.data and obj.data.get(GENERATED)}
    for kind in (bpy.data.meshes, bpy.data.curves, bpy.data.metaballs):
        geometry.update(data for data in kind if data.get(RESOURCE_OWNER) == identifier and not data.users)
    geometry.update(obj['mb_esp_original_mesh'] for obj in objects if 'mb_esp_original_mesh' in obj)
    materials = {slot.material for obj in objects for slot in obj.material_slots if slot.material}
    materials.update(obj['mb_esp_original_material'] for obj in objects if 'mb_esp_original_material' in obj)
    materials.update(mat for data in geometry for mat in data.materials if mat)
    materials.update(mat for mat in bpy.data.materials if mat.get(RESOURCE_OWNER) == identifier)
    profiles = {data.bevel_object for data in geometry if getattr(data, 'bevel_object', None)}
    profiles.update(obj for obj in bpy.data.objects if obj.get(PROFILE) and obj.get(RESOURCE_OWNER) == identifier)
    actions = {block.animation_data.action for block in objects | geometry
               if block.animation_data and block.animation_data.action}
    actions.update(action for action in bpy.data.actions if action.get(RESOURCE_OWNER) == identifier)
    collections = {coll for coll in bpy.data.collections if coll.get(RESOURCE_OWNER) == identifier}
    for obj in objects:
        for coll in obj.users_collection:
            if coll.get(GENERATED) or coll.get('mb_surface_collection') or 'mb_surface_visibility_home' in coll:
                collections.add(coll)
            if 'mb_surface_visibility_home' in coll:
                collections.add(coll['mb_surface_visibility_home'])

    # Keep user decorations and nested molecules where the user placed them.
    for obj in objects:
        for child in obj.children:
            if child not in objects:
                world = child.matrix_world.copy()
                child.parent = owner.parent
                child.matrix_world = world
    bpy.data.batch_remove(objects)
    unused = {block for block in geometry if not block.users and not block.library}
    if unused:
        bpy.data.batch_remove(unused)
    for profile in profiles:
        discard_profile(profile)
    for mat in materials:
        if (not mat.users and not mat.library and
                (mat.get(GENERATED) or mat.get('mb_style_material') or mat.get('mb_surface_material'))):
            bpy.data.materials.remove(mat)
    for action in actions:
        if action.get(GENERATED) and not action.users and not action.library:
            bpy.data.actions.remove(action)
    for coll in sorted(collections, key=lambda c: len(c.children_recursive)):
        if not coll.objects and not coll.children and coll.users <= 1 and not coll.library and not coll.use_fake_user:
            bpy.data.collections.remove(coll)


class MB_OT_remove_molecule(bpy.types.Operator):
    bl_idname = 'molecular_blender.remove_molecule'
    bl_label = 'Remove Molecule'
    bl_description = 'Remove this molecule and its unused generated data'
    bl_options = {'UNDO'}

    @classmethod
    def poll(cls, context):
        return context.mode == 'OBJECT' and resolve_owner(context.object) is not None

    def invoke(self, context, event):
        return context.window_manager.invoke_confirm(self, event)

    def execute(self, context):
        try:
            remove_molecule(resolve_owner(context.object))
        except ValueError as exc:
            self.report({'ERROR'}, str(exc))
            return {'CANCELLED'}
        return {'FINISHED'}
