"""Browse visible fields and add surfaces to an existing molecule."""
import os

import bpy
from bpy.props import StringProperty

from . import import_dialog
from .importers import file_format
from .molecule import Atom, Molecule
from .options import validated_options
from .plotter import process_options
from .scene_edit import _set_color, _set_opacity, _surface_output
from .scene_state import ID, owned_objects, read_record, resolve_owner, write_record
from .surface_jobs import SavedMoldenSource
from .stylers import get_styler
from .surface_operator import DENSITIES, MoldenSurfaceProperties, SurfaceJobOperator


DENSITY_LABELS = {k: name for k, name, _ in DENSITIES}
VISIBILITY_HOME = 'mb_surface_visibility_home'


def surface_groups(owner, record):
    objects = {}
    for obj in owned_objects(owner, 'surface'):
        objects.setdefault(obj.get(ID), []).append(obj)
    groups = {}
    fields = {entry['field'] for entry in record['surfaces']}
    settings = record['settings']
    for entry in record['surfaces']:
        members = objects.get(entry['id'], ())
        field = entry['field']
        # Early scene records used the frame-specific mesh name as the field.
        if 'group' not in entry and settings['plot_type'] == 'animate':
            base, _, suffix = field.rpartition('-')
            if base in fields and len(suffix) >= 4 and suffix.strip().isdigit():
                field = base
        key = entry.get('group', field)
        if key not in groups:
            if field == 'cube':
                label = 'Cube'
                datasets = (record['fields'] or {}).get('datasets', [])
                if len(datasets) > 1:
                    label += ' · dataset {}'.format(datasets[settings['cube_dataset']])
            else:
                label = DENSITY_LABELS.get(field, 'Orbital ' + field.upper())
                spin = entry.get('orbital_spin', settings['orbital_spin'])
                if spin != 'all' and import_dialog.uses_frontier_spin(field):
                    label += ' (' + spin + ')'
            groups[key] = dict(label=label, levels=set(), objects=[], phases={}, ids=set())
        group = groups[key]
        group['ids'].add(entry['id'])
        if members:
            group['levels'].add(entry['isovalue'])
            group['objects'].extend(members)
            group['phases'].setdefault(entry['isovalue'] >= 0, []).extend(members)
    return {key: group for key, group in groups.items() if group['objects']}


def surface_color(obj):
    material = obj.active_material
    return tuple(material.diffuse_color[:3]) if material else (1., 1., 1.)


def surface_opacity(obj):
    material = obj.active_material
    if material and material.use_nodes:
        mix = material.node_tree.nodes.get('MB Surface Opacity')
        if mix:
            return mix.inputs[0].default_value
    return 1.


def _edit_surface_materials(obj, *, color=None, opacity=None):
    owner = resolve_owner(obj)
    if owner is None:
        return
    for group in surface_groups(owner, read_record(owner)).values():
        if obj.get(ID) in group['ids']:
            break
    else:
        return
    objects = group['objects']
    if color is not None:
        objects = next(phase for phase in group['phases'].values() if obj in phase)
    slots_by_material = {}
    for member in objects:
        for slot in member.material_slots:
            if slot.material:
                slots_by_material.setdefault(slot.material, []).append(slot)
    if opacity is not None:
        for material in slots_by_material:
            if material.use_nodes:
                _surface_output(material)
    for material, slots in slots_by_material.items():
        # A data-linked slot can be shared by many objects with a single user
        # count. Switch to object-linked copies before the first edit.
        private = (material.library is None and material.get('mb_style_material') and material.users == len(slots)
                   and all(slot.link == 'OBJECT' for slot in slots))
        target = material if private else material.copy()
        if color is not None:
            _set_color(target, (*color, target.diffuse_color[3]))
        else:
            _set_opacity(target, opacity)
        if not private:
            target['mb_style_material'] = True
            for slot in slots:
                slot.link = 'OBJECT'
                slot.material = target


def set_surface_color(obj, color):
    _edit_surface_materials(obj, color=color)


def set_surface_opacity(obj, opacity):
    _edit_surface_materials(obj, opacity=opacity)


def surface_visible(obj):
    if obj.animation_data:
        return any(VISIBILITY_HOME not in coll for coll in obj.users_collection)
    return not (obj.hide_viewport and obj.hide_render)


def set_surface_visibility(objects, visible):
    # Animated meshes already keyframe both visibility flags. Gate those meshes
    # in a hidden child collection so playback still chooses the correct frame.
    hidden_collections = {}
    emptied = set()
    for obj in objects:
        if not obj.animation_data:
            obj.hide_viewport = obj.hide_render = not visible
            continue
        for coll in obj.users_collection:
            home = coll.get(VISIBILITY_HOME)
            if visible and home is not None:
                if obj.name not in home.objects:
                    home.objects.link(obj)
                coll.objects.unlink(obj)
                emptied.add(coll)
            elif not visible and home is None:
                if coll not in hidden_collections:
                    hidden = bpy.data.collections.new(coll.name + ' (hidden surfaces)')
                    hidden[VISIBILITY_HOME] = coll
                    hidden.hide_viewport = hidden.hide_render = True
                    coll.children.link(hidden)
                    hidden_collections[coll] = hidden
                hidden_collections[coll].objects.link(obj)
                coll.objects.unlink(obj)
    for coll in emptied:
        if not coll.objects and not coll.children:
            bpy.data.collections.remove(coll)


class MB_OT_surface_visibility(bpy.types.Operator):
    bl_idname = 'molecular_blender.surface_visibility'
    bl_label = 'Toggle surface visibility'
    bl_description = 'Show or hide all signs and levels of this field in the viewport and render'
    bl_options = {'UNDO'}

    group: StringProperty(options={'HIDDEN'})

    def execute(self, context):
        owner = resolve_owner(context.object)
        if owner is None:
            return {'CANCELLED'}
        group = surface_groups(owner, read_record(owner)).get(self.group)
        if not group:
            return {'CANCELLED'}
        objects = group['objects']
        visible = not any(surface_visible(obj) for obj in objects)
        if not visible and context.object in objects:
            owner.select_set(True)
            context.view_layer.objects.active = owner
        set_surface_visibility(objects, visible)
        return {'FINISHED'}


class MB_OT_remove_surface(bpy.types.Operator):
    bl_idname = 'molecular_blender.remove_surface'
    bl_label = 'Remove Surface'
    bl_description = 'Remove all lobes, levels, and frames of this calculation'
    bl_options = {'UNDO'}

    group: StringProperty(options={'HIDDEN'})

    def execute(self, context):
        owner = resolve_owner(context.object)
        if owner is None:
            return {'CANCELLED'}
        record = read_record(owner)
        group = surface_groups(owner, record).get(self.group)
        if group is None:
            return {'CANCELLED'}
        objects = group['objects']
        if context.object in objects:
            owner.select_set(True)
            context.view_layer.objects.active = owner
        meshes = {obj.data for obj in objects if obj.type == 'MESH'}
        materials = {slot.material for obj in objects for slot in obj.material_slots if slot.material}
        materials.update(mat for mesh in meshes for mat in mesh.materials if mat)
        actions = {obj.animation_data.action for obj in objects if obj.animation_data and obj.animation_data.action}
        collections = {coll for obj in objects for coll in obj.users_collection}
        hidden = {coll for coll in collections if VISIBILITY_HOME in coll}
        collections.update(coll[VISIBILITY_HOME] for coll in hidden)
        for obj in objects:
            bpy.data.objects.remove(obj, do_unlink=True)
        for mesh in meshes:
            if not mesh.users:
                bpy.data.meshes.remove(mesh)
        for material in materials:
            if not material.users and (material.get('mb_style_material') or material.get('mb_surface_material')):
                bpy.data.materials.remove(material)
        for action in actions:
            if not action.users:
                bpy.data.actions.remove(action)
        homes = collections - hidden
        for coll in [*hidden, *(c for c in homes if c.get('mb_surface_collection'))]:
            if coll.users == 1 and not coll.objects and not coll.children:
                bpy.data.collections.remove(coll)
        record['surfaces'] = [entry for entry in record['surfaces'] if entry['id'] not in group['ids']]
        write_record(owner, record)
        return {'FINISHED'}


def surface_model(owner, record):
    path = record['source']
    if file_format(path) != 'molden' or not record['orbitals']:
        raise ValueError('Adding a surface requires a single-frame Molden molecule with orbitals')
    if not os.path.isfile(path):
        raise ValueError('Molden source unavailable: ' + path)
    atoms = [(a['index'], a['element'], tuple(a['position'])) for a in record['atoms']]
    # Scene placement only needs the original center of mass.
    model = Molecule(owner.name, [Atom(symbol, pos, index) for index, symbol, pos in atoms])
    model.surface_input = SavedMoldenSource(path, atoms, record['orbitals'])
    model.scene_owner = owner
    model.scene_collection = owner.users_collection[0]
    return model


class MB_OT_add_surface(MoldenSurfaceProperties, SurfaceJobOperator, bpy.types.Operator):
    bl_idname = 'molecular_blender.add_surface'
    bl_label = 'Add Surface'
    bl_description = 'Calculate an orbital or density surface on the selected Molden molecule'
    bl_options = {'UNDO'}

    _settings_names = ('volume', 'orbital', 'orbital_spin', 'cumulative', 'isovalues',
                       'resolution', 'remesh', 'orca_fix')

    def invoke(self, context, event):
        self._target = resolve_owner(context.object)
        if self._target is None:
            return {'CANCELLED'}
        settings = read_record(self._target)['settings']
        for name in self._settings_names:
            if not self.properties.is_property_set(name) and (name != 'isovalues' or settings[name]):
                setattr(self, name, settings[name])
        return context.window_manager.invoke_props_dialog(self, width=350)

    def draw(self, context):
        import_dialog.draw_surface_controls(self.layout, self, 'molden')

    def execute(self, context):
        try:
            owner = self._target if hasattr(self, '_target') else resolve_owner(context.object)
            if owner is None:
                raise ValueError('Select a Molden molecule')
            record = read_record(owner)
            settings = validated_options(dict(record['settings'],
                **{name: getattr(self, name) for name in self._settings_names}))
            if not settings['isovalues'].strip():
                raise ValueError('Enter at least one surface level')
            model = surface_model(owner, record)
            styler = get_styler(settings)
            options = process_options(settings)
        except Exception as exc:
            self.report({'ERROR'}, str(exc))
            return {'CANCELLED'}
        return self.start_surface_job(context, model, styler, options, select_surfaces=False)
