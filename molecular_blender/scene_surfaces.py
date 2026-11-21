"""Browse visible fields and add surfaces to an existing molecule."""
import os

import bpy
from bpy.props import StringProperty

from . import import_dialog
from .importers import FORMAT_SUFFIXES, file_format
from .molecule import Atom, Molecule
from .options import molden_convention, validated_options
from .plotter import create_surface_mesh_data, field_geometry, process_options
from .scene_edit import _set_color, _set_opacity, _surface_output
from .scene_state import ID, new_id, surface_sampling, owned_objects, read_record, resolve_owner, write_record
from .surface_jobs import SavedSurfaceSource
from .stylers import get_styler
from .surface_operator import DENSITIES, MoldenSurfaceProperties, SurfaceJobOperator


DENSITY_LABELS = {k: name for k, name, _ in DENSITIES}
VISIBILITY_HOME = 'mb_surface_visibility_home'
SOLO_GROUP = 'mb_solo_surface'
_PRE_SOLO = 'mb_visibility_before_solo'


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
            sampling = entry.get('sampling') or {}
            if field == 'cube':
                label = 'Linear combination' if (record['fields'] or {}).get('linear_combination') else 'Cube'
                datasets = (record['fields'] or {}).get('datasets', [])
                if len(datasets) > 1:
                    index = sampling.get('cube_dataset', settings['cube_dataset'])
                    label += ' · dataset {}'.format(datasets[index])
            elif field.startswith('response '):
                label = 'Response density · eigenpair ' + field.split()[1]
            else:
                label = DENSITY_LABELS.get(field, 'Orbital ' + field.upper())
                spin = entry.get('orbital_spin', settings['orbital_spin'])
                if spin != 'all' and import_dialog.uses_frontier_spin(field):
                    label += ' (' + spin + ')'
            groups[key] = dict(label=label, field=field, spin=entry.get('orbital_spin', settings['orbital_spin']),
                               sampling=entry.get('sampling'), levels=set(), objects=[], phases={}, ids=set())
        group = groups[key]
        if 'ambiguous_phase_frame' in entry:
            group.setdefault('ambiguous_phase_frames', set()).add(entry['ambiguous_phase_frame'])
        if field == 'cube':
            kind = (group['sampling'] or {}).get('cube_field_type', settings['cube_field_type'])
            group['cube_kind'] = (record['fields'] or {}).get('field_type', 'unknown') if kind == 'auto' else kind
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


def visibility_state(obj):
    if obj.animation_data:
        return 0 if surface_visible(obj) else 3
    return int(obj.hide_render) + 2 * int(obj.hide_viewport)


def restore_visibility(obj, state):
    if obj.animation_data:
        set_surface_visibility([obj], state == 0)
    else:
        obj.hide_render = bool(state & 1)
        obj.hide_viewport = bool(state & 2)


def apply_solo(owner, groups):
    if owner[SOLO_GROUP] not in groups:
        end_solo(owner)
        return
    for key, group in groups.items():
        for obj in group['objects']:
            if _PRE_SOLO not in obj:
                obj[_PRE_SOLO] = visibility_state(obj)
        set_surface_visibility(group['objects'], key == owner[SOLO_GROUP])


def end_solo(owner):
    if SOLO_GROUP not in owner:
        return
    for obj in owned_objects(owner, 'surface'):
        if _PRE_SOLO in obj:
            restore_visibility(obj, obj[_PRE_SOLO])
            del obj[_PRE_SOLO]
    del owner[SOLO_GROUP]


class MB_OT_solo_surface(bpy.types.Operator):
    bl_idname = 'molecular_blender.solo_surface'
    bl_label = 'Solo Surface'
    bl_description = 'Show only this surface on this molecule. Click again to restore previous visibility'
    bl_options = {'UNDO'}

    group: StringProperty(options={'HIDDEN', 'SKIP_SAVE'})

    def execute(self, context):
        owner = resolve_owner(context.object)
        if owner is None:
            return {'CANCELLED'}
        groups = surface_groups(owner, read_record(owner))
        if self.group not in groups:
            return {'CANCELLED'}
        if context.object in owned_objects(owner, 'surface'):
            owner.select_set(True)
            context.view_layer.objects.active = owner
        if owner.get(SOLO_GROUP) == self.group:
            end_solo(owner)
        else:
            owner[SOLO_GROUP] = self.group
            apply_solo(owner, groups)
        return {'FINISHED'}


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
        end_solo(owner)
        if not visible and context.object in objects:
            owner.select_set(True)
            context.view_layer.objects.active = owner
        set_surface_visibility(objects, visible)
        return {'FINISHED'}


def remove_surface_objects(objects):
    meshes = {obj.data for obj in objects if obj.type == 'MESH'}
    materials = {slot.material for obj in objects for slot in obj.material_slots if slot.material}
    materials.update(mat for mesh in meshes for mat in mesh.materials if mat)
    meshes.update(obj['mb_esp_original_mesh'] for obj in objects if 'mb_esp_original_mesh' in obj)
    materials.update(obj['mb_esp_original_material'] for obj in objects if 'mb_esp_original_material' in obj)
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
        if owner.get(SOLO_GROUP) == self.group:
            end_solo(owner)
        objects = group['objects']
        if context.object in objects:
            owner.select_set(True)
            context.view_layer.objects.active = owner
        remove_surface_objects(objects)
        record['surfaces'] = [entry for entry in record['surfaces'] if entry['id'] not in group['ids']]
        write_record(owner, record)
        return {'FINISHED'}


def surface_model(owner, record, *, cube=False):
    path = record['source']
    if cube:
        if not can_recompute_cube(record):
            raise ValueError('Choose a single-frame cube surface')
    elif not record['orbitals'] or (file_format(path) != 'molden' and not import_dialog.saved_response_pairs(record)):
        raise ValueError('Adding a surface requires a single-frame Molden molecule with orbitals')
    if not os.path.isfile(path):
        raise ValueError('Surface source unavailable: ' + path)
    source = SavedSurfaceSource.from_record(path, record)
    # Scene placement only needs the original center of mass.
    model = Molecule(owner.name, [Atom(symbol, pos, index) for index, symbol, pos in source.atoms])
    model.surface_input = source
    model.scene_owner = owner
    model.scene_collection = owner.users_collection[0]
    return model


class MB_OT_locate_source(bpy.types.Operator):
    bl_idname = 'molecular_blender.locate_source'
    bl_label = 'Locate Molden Source'
    bl_description = 'Check the selected Molden file against the saved atoms and orbital metadata, then relink it'
    bl_options = {'UNDO'}

    filepath: StringProperty(subtype='FILE_PATH', options={'SKIP_SAVE'})
    filter_glob: StringProperty(default=';'.join('*' + s for s in FORMAT_SUFFIXES['molden']), options={'HIDDEN'})

    def invoke(self, context, event):
        self._owner = resolve_owner(context.object)
        if self._owner is None:
            return {'CANCELLED'}
        self.filepath = read_record(self._owner)['source']
        context.window_manager.fileselect_add(self)
        return {'RUNNING_MODAL'}

    def execute(self, context):
        try:
            owner = self._owner if hasattr(self, '_owner') else resolve_owner(context.object)
            if owner is None:
                raise ValueError('Select an imported Molden molecule')
            record = read_record(owner)
            path = os.path.abspath(bpy.path.abspath(self.filepath))
            if file_format(record['source']) != 'molden' or file_format(path) != 'molden':
                raise ValueError('Choose a Molden source file for a Molden molecule')
            SavedSurfaceSource.from_record(path, record).load(record['settings'])
            record['source'] = path
            write_record(owner, record)
        except Exception as exc:
            self.report({'ERROR'}, str(exc) or 'Could not read the Molden source')
            return {'CANCELLED'}
        self.report({'INFO'}, 'Linked ' + os.path.basename(path))
        return {'FINISHED'}


class MB_OT_add_surface(MoldenSurfaceProperties, SurfaceJobOperator, bpy.types.Operator):
    bl_idname = 'molecular_blender.add_surface'
    bl_label = 'Add Surface'
    bl_description = 'Calculate another surface from the selected molecule’s source file'
    bl_options = {'UNDO'}

    _molden_names = ('volume', 'orbital', 'orbital_spin', 'cumulative', 'isovalues',
                     'resolution', 'remesh', 'molden_convention')
    _cube_names = ('cube_dataset', 'cube_field_type', 'cube_both_signs', 'cumulative', 'isovalues', 'remesh')

    def invoke(self, context, event):
        self._target = resolve_owner(context.object)
        if self._target is None:
            return {'CANCELLED'}
        record = read_record(self._target)
        settings = record['settings']
        cube = can_recompute_cube(record)
        pairs = import_dialog.saved_response_pairs(record)
        if pairs and not self.response_eigenpair:
            self.response_eigenpair = pairs[0][0]
        names = self._cube_names if cube else self._molden_names
        for name in names:
            if self.properties.is_property_set(name):
                continue
            if name == 'cube_both_signs':
                self.cube_both_signs = settings.get(name, settings['volume'] in ('orbital', 'spin', 'response'))
            elif (name != 'isovalues' or settings[name]) and (name != 'volume' or not pairs):
                setattr(self, name, settings[name])
        return context.window_manager.invoke_props_dialog(self, width=350)

    def draw(self, context):
        record = read_record(self._target)
        if can_recompute_cube(record):
            label = record['fields']['datasets'][self.cube_dataset]
            self.layout.label(text=f'Cube dataset {label}')
            import_dialog.draw_sampling_controls(self.layout, self, 'cube')
        elif import_dialog.saved_response_pairs(record):
            self.layout.label(text=f'Response density · eigenpair {self.response_eigenpair}')
            import_dialog.draw_sampling_controls(self.layout, self, 'molden')
        else:
            path = record['source'] if self.surface_advanced else None
            import_dialog.draw_surface_controls(self.layout, self, 'molden', path=path)

    def execute(self, context):
        try:
            owner = self._target if hasattr(self, '_target') else resolve_owner(context.object)
            if owner is None:
                raise ValueError('Select an imported molecule')
            record = read_record(owner)
            cube = can_recompute_cube(record)
            pairs = import_dialog.saved_response_pairs(record)
            names = self._cube_names if cube else self._molden_names
            if cube or pairs:
                overrides = {k: getattr(self, k) for k in names if self.properties.is_property_set(k)}
                if not record['settings']['isovalues']:
                    overrides.setdefault('isovalues', self.isovalues)
                if cube:
                    index = overrides.get('cube_dataset', record['settings']['cube_dataset'])
                    if index >= len(record['fields']['datasets']):
                        raise ValueError('Cube dataset index is out of range')
                else:
                    pair = self.response_eigenpair or pairs[0][0]
                    if pair not in {p for p, _ in pairs}:
                        raise ValueError(f'DIPL eigenpair {pair} is unavailable')
                    overrides.update(volume='response', response_eigenpair=pair)
            else:
                overrides = {k: getattr(self, k) for k in names}
                if not self.properties.is_property_set('molden_convention') and not self.properties.is_property_set('orca_fix'):
                    overrides['molden_convention'] = record['settings']['molden_convention']
            return start_field_job(self, context, owner, record, overrides, cube=cube)
        except (KeyError, ValueError) as exc:
            self.report({'ERROR'}, str(exc))
            return {'CANCELLED'}


def can_recompute_cube(record):
    kind = file_format(record['source'])
    fields = record.get('fields')
    if fields is None or kind not in ('cube', 'json'):
        return False
    if 'frames' in fields:
        return fields['frames'] == 1
    return kind == 'cube' or record['settings']['plot_type'] != 'animate'


def start_field_job(op, context, owner, record, overrides, *, cube=False):
    settings = validated_options(dict(record['settings'], **overrides))
    if not settings['isovalues'].strip():
        raise ValueError('Enter at least one surface level')
    return op.start_surface_job(context, surface_model(owner, record, cube=cube), get_styler(settings),
                                process_options(settings, cube=cube), select_surfaces=False)


def recompute_settings(record, group):
    field = group['field']
    if field == 'cube':
        # Legacy cubes stored their requested levels in the molecule settings.
        settings = record['settings']
        if group['sampling']:
            return dict(group['sampling'])
        return dict(cumulative=settings['cumulative'], isovalues=settings['isovalues'], resolution=settings['resolution'],
                    cube_both_signs=settings.get('cube_both_signs', settings['volume'] in ('orbital', 'spin')))
    # Older files saved the absolute mesh levels, even for fractional requests.
    signed = field not in DENSITY_LABELS or field == 'spin'
    sampling = group['sampling'] or dict(cumulative=False, isovalues=', '.join(
        str(x) for x in sorted({abs(v) if signed else v for v in group['levels']})),
        resolution=record['settings']['resolution'], molden_convention=molden_convention(record['settings']))
    if field.startswith('response '):
        return dict(sampling, volume='response', response_eigenpair=int(field.split()[1]))
    return dict(sampling, molden_convention=molden_convention(sampling), volume=field if field in DENSITY_LABELS else 'orbital',
                orbital=field if field not in DENSITY_LABELS else 'homo', orbital_spin=group['spin'])


def replace_surface_geometry(context, owner, record, key, group, definitions, options, styler):
    """Prepare the complete replacement before swapping any existing mesh."""
    from .scene_esp import BASE_MESH, remap_contour
    entries = {entry['id']: entry for entry in record['surfaces'] if entry['id'] in group['ids']}
    old_objects = group['objects']
    remaining = set(old_objects)
    sampling = surface_sampling(options, cube=group['field'] == 'cube')
    staged, clones, new_meshes, new_materials = [], [], [], []
    replacements = []
    grids, backups = {}, []
    try:
        # Match unchanged levels first so named objects and external references survive.
        def match_order(item):
            level = item['isovalue']
            return (level not in group['levels'], -abs(level))
        for item in sorted(definitions, key=match_order):
            level = item['isovalue']
            phase = group['phases'].get(level >= 0, [])

            def distance(obj):
                return abs(entries[obj[ID]]['isovalue'] - level)

            obj = min((o for o in phase if o in remaining), key=distance, default=None)
            template = obj or min(phase, key=distance, default=old_objects[0])
            mesh = create_surface_mesh_data(template.data.name, *field_geometry(item, options['resolution']),
                                            weld_distance=None)
            new_meshes.append(mesh)
            old_base = template.get(BASE_MESH)
            for mat in (old_base if old_base is not None else template.data).materials:
                mesh.materials.append(mat)
            if obj is None:
                obj = template.copy()
                clones.append(obj)
                obj[ID] = new_id()
                obj.name = '{}_{}_{}'.format(owner.name, group['field'], level)
                for coll in template.users_collection:
                    coll.objects.link(obj)
            else:
                remaining.remove(obj)
            if old_base is not None:
                base = mesh
                mesh = remap_contour(context, record, template, base, grids)
                new_meshes.append(mesh)
                backups.append((obj, old_base, base))
            elif not phase:
                material = styler.outer_isosurface_material('orbital_mat_' + ('plus' if level >= 0 else 'negative'))
                new_materials.append(material)
                material['mb_surface_material'] = True
                opacity = surface_opacity(template)
                if opacity != 1.:
                    _set_opacity(material, opacity)
                if not mesh.materials:
                    mesh.materials.append(material)
                obj.data = mesh  # only a newly created sign can enter here
                obj.material_slots[0].link = 'OBJECT'
                obj.material_slots[0].material = material
            staged.append((obj, obj.data, mesh))
            replacements.append(dict(id=obj[ID], field=group['field'], group=key, orbital_spin=group['spin'],
                                     isovalue=float(level), sampling=sampling))
        for obj, old_mesh, mesh in staged:
            obj.data = mesh
        for obj, old, base in backups:
            obj[BASE_MESH] = base
        surfaces = record['surfaces']
        first = next(i for i, entry in enumerate(surfaces) if entry['id'] in entries)
        record['surfaces'] = (surfaces[:first] + replacements
                              + [e for e in surfaces[first:] if e['id'] not in entries])
        write_record(owner, record)
    except Exception:
        for obj, old, base in backups:
            obj[BASE_MESH] = old
        for obj, old_mesh, mesh in staged:
            if obj not in clones:
                obj.data = old_mesh
        for obj in clones:
            bpy.data.objects.remove(obj, do_unlink=True)
        for mesh in new_meshes:
            bpy.data.meshes.remove(mesh)
        for material in new_materials:
            bpy.data.materials.remove(material)
        raise
    for old_mesh in {old for _, old, _ in staged + backups}:
        if old_mesh.users == 0:
            bpy.data.meshes.remove(old_mesh)
    if context.view_layer.objects.active in remaining:
        target = staged[0][0]
        target.select_set(True)
        context.view_layer.objects.active = target
    remove_surface_objects(remaining)


class MB_OT_recompute_surface(MoldenSurfaceProperties, SurfaceJobOperator, bpy.types.Operator):
    bl_idname = 'molecular_blender.recompute_surface'
    bl_label = 'Recompute Surface'
    bl_description = 'Recalculate this field at new levels or spacing and keep its styling'
    bl_options = {'REGISTER', 'UNDO'}

    group: StringProperty(options={'HIDDEN', 'SKIP_SAVE'})

    def _group_settings(self, owner):
        if owner is None:
            raise ValueError('Select a molecule with a saved surface')
        record = read_record(owner)
        group = surface_groups(owner, record).get(self.group)
        if group is None:
            raise ValueError('The surface calculation no longer exists')
        return record, group, recompute_settings(record, group)

    def invoke(self, context, event):
        self._target = resolve_owner(context.object)
        try:
            _, group, settings = self._group_settings(self._target)
        except ValueError as exc:
            self.report({'ERROR'}, str(exc))
            return {'CANCELLED'}
        self._label = group['label']
        self._cube = group['field'] == 'cube'
        for name in ('cumulative', 'isovalues', 'resolution'):
            setattr(self, name, settings[name])
        if self._cube:
            self.cube_both_signs = settings['cube_both_signs']
        return context.window_manager.invoke_props_dialog(self)

    def draw(self, context):
        self.layout.label(text=getattr(self, '_label', self.bl_label))
        self.layout.prop(self, 'level_mode')
        self.layout.prop(self, 'isovalues', text='Fractions' if self.cumulative else 'Isovalues')
        if getattr(self, '_cube', False):
            self.layout.prop(self, 'cube_both_signs')
        else:
            self.layout.prop(self, 'resolution', text='Grid spacing (Å)')

    def execute(self, context):
        owner = self._target if hasattr(self, '_target') else resolve_owner(context.object)
        try:
            record, group, settings = self._group_settings(owner)
            for name in ('cumulative', 'isovalues', 'resolution'):
                if self.properties.is_property_set(name):
                    settings[name] = getattr(self, name)
            if group['field'] == 'cube' and self.properties.is_property_set('cube_both_signs'):
                settings['cube_both_signs'] = self.cube_both_signs
            self._original_meshes = {o.data.as_pointer() for o in group['objects']}
            return start_field_job(self, context, owner, record, settings, cube=group['field'] == 'cube')
        except (ValueError, KeyError) as exc:
            self.report({'ERROR'}, str(exc))
            return {'CANCELLED'}

    def apply_surface_result(self, context, vertex_sets, trajectory_map):
        record, group, _ = self._group_settings(self._owner)
        if {o.data.as_pointer() for o in group['objects']} != self._original_meshes:
            raise ValueError('The surface changed while it was being calculated')
        replace_surface_geometry(context, self._owner, record, self.group, group,
                                 vertex_sets, self._options, self._styler)
