"""Map a supplied scalar cube onto the displayed mesh of a density surface."""
import math
import os

import bpy
import numpy as np
from bpy.props import EnumProperty, FloatProperty, IntProperty, StringProperty

from .esp import sample_cube
from .import_dialog import browser_path
from .importers import file_format, molecule_from_cube, molecule_from_json
from .scene_edit import _set_opacity, _surface_output
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


def scale_nodes(material):
    # Follow the active shader so renamed nodes and unrelated ramps are safe.
    if material is None or not material.use_nodes:
        raise ValueError("The field material needs shader nodes")
    shader = _surface_output(material).inputs['Surface'].links[0].from_node
    if shader.type == 'MIX_SHADER' and shader.inputs[2].is_linked:
        shader = shader.inputs[2].links[0].from_node  # The mapped surface's opacity wrapper.
    if shader.type != 'BSDF_PRINCIPLED' or not shader.inputs['Base Color'].is_linked:
        raise ValueError('The mapped surface must use a Principled shader with a connected Base Color')
    ramp = shader.inputs['Base Color'].links[0].from_node
    if ramp.type == 'VALTORGB' and ramp.inputs['Fac'].is_linked:
        scale = ramp.inputs['Fac'].links[0].from_node
        if scale.type == 'MAP_RANGE' and scale.data_type == 'FLOAT' and scale.interpolation_type == 'LINEAR':
            bounds = [scale.inputs[key] for key in ('From Min', 'From Max', 'To Min', 'To Max')]
            low, high, start, end = (socket.default_value for socket in bounds)
            if not any(socket.is_linked for socket in bounds) and low < high and (start, end) == (0., 1.):
                return ramp, scale
    raise ValueError('Field scale needs a Color Ramp connected to Base Color with a linear Map Range from 0 to 1')


def field_map_settings(mat):
    if mat is None:
        return None
    kind = mat.get('mb_field_kind', 'esp')
    prefix = 'mb_esp_' if kind == 'esp' else 'mb_field_'
    if prefix + 'limits' not in mat:
        return None
    try:
        _, scale = scale_nodes(mat)
        limits = tuple(scale.inputs[key].default_value for key in ('From Min', 'From Max'))
    except ValueError:
        limits = mat[prefix + 'limits']
    return dict(kind=kind, label=mat.get('mb_field_label', 'ESP'),
                path=mat[prefix + 'path'], dataset=mat[prefix + 'dataset'],
                units=mat[prefix + 'units'], limits=limits)


def validate_limits(kind, low, high):
    if not (math.isfinite(low) and math.isfinite(high) and low < high):
        raise ValueError('Color limits must be finite and increase from minimum to maximum')
    if kind == 'unsigned' and low < 0:
        raise ValueError('The minimum for an unsigned field must be nonnegative')
    if kind != 'unsigned' and not low < 0 < high:
        raise ValueError('Signed color limits must span zero')


def _field_material(path, dataset, units, low, high, opacity, kind, label):
    mat = bpy.data.materials.new(label + ' · ' + os.path.basename(path))
    try:
        mat['mb_surface_material'] = True
        mat['mb_field_kind'] = kind
        mat['mb_field_label'] = label
        prefix = 'mb_esp_' if kind == 'esp' else 'mb_field_'
        mat[prefix + 'path'] = path
        mat[prefix + 'dataset'] = dataset
        mat[prefix + 'units'] = units
        mat[prefix + 'limits'] = (low, high)
        mat.use_nodes = True
        tree = mat.node_tree
        shader = tree.nodes.get('Principled BSDF')
        attr = tree.nodes.new('ShaderNodeAttribute')
        attr.attribute_name = ATTRIBUTE if kind == 'esp' else 'mb_field'
        attr.location = (-720, 260)
        limits = tree.nodes.new('ShaderNodeMapRange')
        limits.label = label + ' (' + units + ')'
        limits.location = (-500, 260)
        limits.inputs['From Min'].default_value = low
        limits.inputs['From Max'].default_value = high
        limits.clamp = True
        ramp = tree.nodes.new('ShaderNodeValToRGB')
        ramp.location = (-260, 260)
        ramp.color_ramp.elements[1].color = (.025, .12, .85, 1)
        if kind == 'unsigned':
            ramp.label = 'Low to high'
            ramp.color_ramp.elements[0].color = (1., 1., 1., 1.)
        else:
            ramp.label = 'Negative · zero · positive'
            ramp.color_ramp.elements[0].color = (0.8, .025, .015, 1)
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


def _display_mesh(context, obj, replacement=None):
    """Evaluate only the imported surface's remesh, even if the object is hidden."""
    if any(mod.type != 'REMESH' for mod in obj.modifiers) or len(obj.modifiers) > 1:
        raise ValueError('Field mapping supports the imported surface and its Remesh modifier')
    data = replacement if replacement is not None else obj.data
    if not obj.modifiers or (replacement is None and BASE_MESH in obj):
        return data.copy()
    mod = obj.modifiers[0]
    if mod.show_viewport != mod.show_render:
        raise ValueError('Set Remesh to the same visibility for viewport and render before mapping')
    enabled = obj.get(REMESH, mod.show_viewport) if replacement is not None else mod.show_viewport
    if not enabled:
        return data.copy()
    # Local mesh coordinates are source angstroms. Parent and scene placement are presentation.
    proxy = obj.copy()
    try:
        proxy.parent = None
        proxy.hide_viewport = False
        proxy.data = data
        proxy.modifiers[0].show_viewport = True
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


def _mapping_grid(record, path, dataset):
    if file_format(path) == 'json':
        cube = molecule_from_json(path, {}, single_combination=True)
    else:
        cube = molecule_from_cube(path, {'cube_dataset': dataset})
    if ([atom['symbol'] for atom in cube['atoms']] != [a['element'] for a in record['atoms']]
            or not np.allclose([a['position'] for a in cube['atoms']], [a['position'] for a in record['atoms']],
                               atol=1e-4, rtol=0)):
        raise ValueError('Cube atom order and source coordinates must match the molecule (0.0001 Å tolerance)')
    return cube['volume']


def _sample_mesh(mesh, grid, name):
    points = np.empty(len(mesh.vertices) * 3, dtype=np.float32)
    mesh.vertices.foreach_get('co', points)
    attr = mesh.attributes.get(name) or mesh.attributes.new(name, 'FLOAT', 'POINT')
    attr.data.foreach_set('value', sample_cube(grid, points))


def remap_contour(context, record, template, base, grids):
    """Resample the original map on a new displayed mesh, using the edited material."""
    material = template.active_material
    settings = field_map_settings(material)
    if not settings:
        raise ValueError('Cannot find the source field for the mapped surface')
    key = (settings['path'], settings['dataset'])
    if key not in grids:
        grids[key] = _mapping_grid(record, *key)
    mesh = _display_mesh(context, template, replacement=base)
    try:
        _sample_mesh(mesh, grids[key], ATTRIBUTE if settings['kind'] == 'esp' else 'mb_field')
        mesh.materials.clear()
        mesh.materials.append(material)
        return mesh
    except Exception:
        bpy.data.meshes.remove(mesh)
        raise


def map_esp(context, record, group, path, dataset, units, low, high, *, kind='esp', label=None):
    if not density_group(group):
        raise ValueError('Select a single-frame density surface')
    source_format = file_format(path)
    if source_format == 'json':
        kind, dataset = 'signed', 0
    elif source_format != 'cube':
        raise ValueError('Choose a field cube or a JSON linear combination')
    if label is None:
        label = 'Linear combination' if source_format == 'json' else 'ESP'
    if kind not in ('esp', 'signed', 'unsigned') or not label.strip() or (source_format == 'cube' and not units.strip()):
        raise ValueError('Choose a field type and enter its name and scalar units')
    validate_limits(kind, low, high)
    grid = _mapping_grid(record, path, dataset)
    if source_format == 'json':
        units = grid['units']
    staged, meshes, materials = [], [], {}
    try:
        for obj in group['objects']:
            if obj.type != 'MESH' or len(obj.material_slots) != 1 or obj.active_material is None:
                raise ValueError('Field mapping requires a surface mesh with one material')
            mesh = _display_mesh(context, obj)
            meshes.append(mesh)
            name = ATTRIBUTE if kind == 'esp' else 'mb_field'
            if BASE_MESH in obj:
                previous = ATTRIBUTE if obj.active_material.get('mb_field_kind', 'esp') == 'esp' else 'mb_field'
                if previous != name and mesh.attributes.get(previous):
                    mesh.attributes.remove(mesh.attributes[previous])
            _sample_mesh(mesh, grid, name)
            opacity = surface_opacity(obj)
            if opacity not in materials:
                materials[opacity] = _field_material(path, dataset, units, low, high, opacity, kind, label)
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
    bl_label = 'Map Field onto Density'
    bl_description = 'Color a density surface with a cube or a JSON linear combination of cubes'
    bl_options = {'UNDO'}
    group: StringProperty(options={'HIDDEN', 'SKIP_SAVE'})
    filepath: StringProperty(subtype='FILE_PATH', options={'SKIP_SAVE'})
    filter_glob: StringProperty(default='*.cube;*.cub;*.json', options={'HIDDEN'})
    dataset: IntProperty(name='Cube dataset index', default=0, min=0)

    def _update_kind(self, context):
        if self.field_kind == 'unsigned' and self.minimum < 0:
            self.minimum = 0.
        elif self.field_kind != 'unsigned' and self.minimum >= 0:
            self.minimum = -max(self.maximum, .05)

    field_kind: EnumProperty(name='Color by', items=(
        ('esp', 'Electrostatic Potential', 'Map a supplied ESP cube'),
        ('signed', 'Signed Field', 'Red for negative, white at zero, blue for positive'),
        ('unsigned', 'Unsigned Field', 'White at the lower limit and blue at the upper limit')),
        update=_update_kind)
    field_label: StringProperty(name='Field name', default='Scalar field')
    field_units: StringProperty(name='Scalar units', description='Units of the supplied data, or 1 for dimensionless')
    units: EnumProperty(name='Potential values in', items=(
        ('Hartree/e', 'Hartree/e', 'Atomic units of electrostatic potential'),
        ('V', 'Volts', 'The supplied cube already contains volts')))
    minimum: FloatProperty(name='Minimum', default=-.05, precision=4)
    maximum: FloatProperty(name='Maximum', default=.05, precision=4)

    def invoke(self, context, event):
        self._owner = resolve_owner(context.object)
        if self._owner is None:
            return {'CANCELLED'}
        group = surface_groups(self._owner, read_record(self._owner)).get(self.group)
        if not group:
            return {'CANCELLED'}
        saved = field_map_settings(group['objects'][0].active_material)
        if saved:
            self.filepath, self.dataset, self.field_kind = saved['path'], saved['dataset'], saved['kind']
            if self.field_kind == 'esp':
                self.units = saved['units']
            else:
                self.field_label, self.field_units = saved['label'], saved['units']
            self.minimum, self.maximum = saved['limits']
        context.window_manager.fileselect_add(self)
        return {'RUNNING_MODAL'}

    def draw(self, context):
        if file_format(browser_path(self, context)) == 'json':
            self.layout.label(text='Signed field · units and datasets from JSON')
            fields = ('field_label',)
        else:
            self.layout.prop(self, 'field_kind')
            self.layout.prop(self, 'dataset')
            fields = ('units',) if self.field_kind == 'esp' else ('field_label', 'field_units')
        for key in (*fields, 'minimum', 'maximum'):
            self.layout.prop(self, key)

    def execute(self, context):
        try:
            owner = self._owner if hasattr(self, '_owner') else resolve_owner(context.object)
            if owner is None:
                raise ValueError('Select a density surface or its molecule')
            record = read_record(owner)
            group = surface_groups(owner, record)[self.group]
            path = os.path.abspath(bpy.path.abspath(self.filepath))
            esp = file_format(path) != 'json' and self.field_kind == 'esp'
            units = self.units if esp else self.field_units.strip()
            label = 'ESP' if esp else self.field_label.strip()
            map_esp(context, record, group, path, self.dataset, units, self.minimum, self.maximum,
                    kind=self.field_kind, label=label)
        except (ValueError, OSError, EOFError, KeyError, RuntimeError, ReferenceError) as exc:
            self.report({'ERROR'}, str(exc) or 'Cannot map the field onto the surface')
            return {'CANCELLED'}
        return {'FINISHED'}


class MB_OT_clear_esp(bpy.types.Operator):
    bl_idname = 'molecular_blender.clear_esp'
    bl_label = 'Remove Field Map'
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
