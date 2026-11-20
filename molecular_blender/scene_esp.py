"""Map an external ESP cube onto the displayed mesh of a density surface."""
import math
import os

import bpy
import numpy as np
from bpy.props import EnumProperty, FloatProperty, IntProperty, StringProperty

from .esp import sample_cube
from .importers import file_format, molecule_from_cube
from .scene_edit import _set_opacity
from .scene_state import read_record, resolve_owner
from .scene_surfaces import surface_groups, surface_opacity


ATTRIBUTE = 'mb_esp'
BASE_MESH = 'mb_esp_original_mesh'
BASE_MATERIAL = 'mb_esp_original_material'
BASE_SLOT = 'mb_esp_original_slot'
REMESH = 'mb_esp_original_remesh'


def density_group(group):
    density = (group['field'] in ('density', 'alpha', 'beta') or
               group['field'] == 'cube' and group['cube_kind'] in ('density', 'unknown'))
    return density and not any(
        obj.animation_data and obj.animation_data.action for obj in group['objects'])


def _esp_material(path, dataset, units, low, high, opacity):
    mat = bpy.data.materials.new('ESP · ' + os.path.basename(path))
    try:
        mat['mb_surface_material'] = True
        mat['mb_esp_path'] = path
        mat['mb_esp_dataset'] = dataset
        mat['mb_esp_units'] = units
        mat['mb_esp_limits'] = (low, high)
        mat.use_nodes = True
        tree = mat.node_tree
        shader = tree.nodes.get('Principled BSDF')
        attr = tree.nodes.new('ShaderNodeAttribute')
        attr.attribute_name = ATTRIBUTE
        attr.location = (-720, 260)
        limits = tree.nodes.new('ShaderNodeMapRange')
        limits.label = 'ESP in ' + units
        limits.location = (-500, 260)
        limits.inputs['From Min'].default_value = low
        limits.inputs['From Max'].default_value = high
        limits.clamp = True
        ramp = tree.nodes.new('ShaderNodeValToRGB')
        ramp.label = 'Negative · zero · positive'
        ramp.location = (-260, 260)
        ramp.color_ramp.elements[0].color = (0.8, .025, .015, 1)
        ramp.color_ramp.elements[1].color = (.025, .12, .85, 1)
        ramp.color_ramp.elements.new(-low / (high - low)).color = (1., 1., 1., 1.)
        tree.links.new(attr.outputs['Fac'], limits.inputs['Value'])
        tree.links.new(limits.outputs['Result'], ramp.inputs['Fac'])
        tree.links.new(ramp.outputs['Color'], shader.inputs['Base Color'])
        if opacity != 1.:
            _set_opacity(mat, opacity)
        return mat
    except Exception:
        bpy.data.materials.remove(mat)
        raise


def _display_mesh(context, obj):
    """Evaluate only the imported surface's remesh, even if the object is hidden."""
    if any(mod.type != 'REMESH' for mod in obj.modifiers) or len(obj.modifiers) > 1:
        raise ValueError('ESP mapping supports the imported surface and its Remesh modifier')
    if BASE_MESH in obj or not obj.modifiers:
        return obj.data.copy()
    mod = obj.modifiers[0]
    if mod.show_viewport != mod.show_render:
        raise ValueError('Set Remesh to the same visibility for viewport and render before mapping')
    if not mod.show_viewport:
        return obj.data.copy()
    # Local mesh coordinates are source angstroms. Parent and scene placement are presentation.
    proxy = obj.copy()
    try:
        proxy.parent = None
        proxy.hide_viewport = False
        context.scene.collection.objects.link(proxy)
        depsgraph = context.evaluated_depsgraph_get()
        return bpy.data.meshes.new_from_object(proxy.evaluated_get(depsgraph), depsgraph=depsgraph)
    finally:
        bpy.data.objects.remove(proxy, do_unlink=True)


def _discard_unused(meshes, materials):
    for mesh in meshes:
        if not mesh.users:
            bpy.data.meshes.remove(mesh)
    for mat in materials:
        if mat and not mat.users and (mat.get('mb_style_material') or mat.get('mb_surface_material')):
            bpy.data.materials.remove(mat)


def map_esp(context, record, group, path, dataset, units, low, high):
    if not density_group(group):
        raise ValueError('Select a single-frame density surface')
    if not (math.isfinite(low) and math.isfinite(high) and low < 0 < high):
        raise ValueError('ESP limits must be finite and span zero')
    if file_format(path) != 'cube':
        raise ValueError('Choose a cube containing electrostatic potential')
    cube = molecule_from_cube(path, {'cube_dataset': dataset})
    if ([atom['symbol'] for atom in cube['atoms']] != [a['element'] for a in record['atoms']]
            or not np.allclose([a['position'] for a in cube['atoms']], [a['position'] for a in record['atoms']],
                               atol=1e-4, rtol=0)):
        raise ValueError('ESP cube atom order and source coordinates must match the molecule (0.0001 Å tolerance)')
    staged, meshes, materials = [], [], {}
    try:
        for obj in group['objects']:
            if obj.type != 'MESH' or len(obj.material_slots) != 1 or obj.active_material is None:
                raise ValueError('ESP mapping requires a surface mesh with one material')
            mesh = _display_mesh(context, obj)
            meshes.append(mesh)
            points = np.empty(len(mesh.vertices) * 3, dtype=np.float32)
            mesh.vertices.foreach_get('co', points)
            attr = mesh.attributes.get(ATTRIBUTE) or mesh.attributes.new(ATTRIBUTE, 'FLOAT', 'POINT')
            attr.data.foreach_set('value', sample_cube(cube['volume'], points))
            opacity = surface_opacity(obj)
            if opacity not in materials:
                materials[opacity] = _esp_material(path, dataset, units, low, high, opacity)
            staged.append((obj, mesh, materials[opacity]))
        if not any(len(mesh.vertices) for mesh in meshes):
            raise ValueError('The selected density surface has no vertices')
    except Exception:
        _discard_unused(meshes, materials.values())
        raise
    old_meshes, old_materials = set(), set()
    for obj, mesh, mat in staged:
        slot = obj.material_slots[0]
        if BASE_MESH not in obj:
            obj[BASE_MESH] = obj.data
            obj[BASE_MATERIAL] = slot.material
            obj[BASE_SLOT] = slot.link
            if obj.modifiers:
                mod = obj.modifiers[0]
                obj[REMESH] = mod.show_viewport
                mod.show_viewport = mod.show_render = False
        else:
            old_meshes.add(obj.data)
            old_materials.add(slot.material)
            old_materials.update(obj.data.materials)
        obj.data = mesh
        # Evaluated meshes can lose empty slots. Keep the original slot on the backup mesh.
        mesh.materials.clear()
        mesh.materials.append(mat)
        obj.material_slots[0].link = 'DATA'
    _discard_unused(old_meshes, old_materials)


def clear_esp(group):
    meshes, materials = set(), set()
    for obj in group['objects']:
        if BASE_MESH not in obj:
            continue
        meshes.add(obj.data)
        materials.add(obj.active_material)
        materials.update(obj.data.materials)
        obj.data = obj[BASE_MESH]
        obj.material_slots[0].link = obj[BASE_SLOT]
        if obj[BASE_SLOT] == 'OBJECT':
            obj.material_slots[0].material = obj[BASE_MATERIAL]
        if REMESH in obj:
            if obj.modifiers:
                obj.modifiers[0].show_viewport = obj.modifiers[0].show_render = bool(obj[REMESH])
            del obj[REMESH]
        for key in (BASE_MESH, BASE_MATERIAL, BASE_SLOT):
            del obj[key]
    _discard_unused(meshes, materials)


class MB_OT_map_esp(bpy.types.Operator):
    bl_idname = 'molecular_blender.map_esp'
    bl_label = 'Map ESP Cube onto Density'
    bl_description = 'Color the density surface from a potential cube in the same molecular coordinate system'
    bl_options = {'UNDO'}
    group: StringProperty(options={'HIDDEN', 'SKIP_SAVE'})
    filepath: StringProperty(subtype='FILE_PATH', options={'SKIP_SAVE'})
    filter_glob: StringProperty(default='*.cube;*.cub', options={'HIDDEN'})
    dataset: IntProperty(name='Cube dataset index', default=0, min=0)
    units: EnumProperty(name='Potential values in', items=(
        ('Hartree/e', 'Hartree/e', 'Atomic units of electrostatic potential'),
        ('V', 'Volts', 'The supplied cube already contains volts')))
    minimum: FloatProperty(name='Negative limit', default=-.05, precision=4)
    maximum: FloatProperty(name='Positive limit', default=.05, precision=4)

    def invoke(self, context, event):
        self._owner = resolve_owner(context.object)
        if self._owner is None:
            return {'CANCELLED'}
        group = surface_groups(self._owner, read_record(self._owner)).get(self.group)
        if not group:
            return {'CANCELLED'}
        mat = group['objects'][0].active_material
        if mat and 'mb_esp_path' in mat:
            self.filepath, self.dataset, self.units = mat['mb_esp_path'], mat['mb_esp_dataset'], mat['mb_esp_units']
            self.minimum, self.maximum = mat['mb_esp_limits']
        context.window_manager.fileselect_add(self)
        return {'RUNNING_MODAL'}

    def draw(self, context):
        for key in ('dataset', 'units', 'minimum', 'maximum'):
            self.layout.prop(self, key)

    def execute(self, context):
        try:
            owner = self._owner if hasattr(self, '_owner') else resolve_owner(context.object)
            if owner is None:
                raise ValueError('Select a density surface or its molecule')
            record = read_record(owner)
            group = surface_groups(owner, record)[self.group]
            path = os.path.abspath(bpy.path.abspath(self.filepath))
            map_esp(context, record, group, path, self.dataset, self.units, self.minimum, self.maximum)
        except (ValueError, OSError, EOFError, KeyError, RuntimeError, ReferenceError) as exc:
            self.report({'ERROR'}, str(exc) or 'Cannot map ESP onto the surface')
            return {'CANCELLED'}
        return {'FINISHED'}


class MB_OT_clear_esp(bpy.types.Operator):
    bl_idname = 'molecular_blender.clear_esp'
    bl_label = 'Remove ESP Map'
    bl_description = 'Restore the original surface material, geometry and remeshing'
    bl_options = {'UNDO'}
    group: StringProperty(options={'HIDDEN'})

    def execute(self, context):
        owner = resolve_owner(context.object)
        group = surface_groups(owner, read_record(owner)).get(self.group) if owner else None
        if not group:
            return {'CANCELLED'}
        clear_esp(group)
        return {'FINISHED'}
