"""Edit a mapped surface's shader scale without reopening the scalar field."""
import math

import bpy
from bpy.props import CollectionProperty, EnumProperty, FloatProperty, StringProperty

from .field_palettes import PALETTES, KEEP_PALETTE, apply_palette
from .scene_esp import BASE_MESH, field_map_settings, scale_nodes, validate_limits
from .scene_state import RECORD, read_record, resolve_owner
from .scene_surfaces import surface_groups


def scale_values(material):
    node, scale = scale_nodes(material)
    ramp = node.color_ramp
    low, high = (scale.inputs[key].default_value for key in ('From Min', 'From Max'))
    stops = tuple((stop.position, tuple(stop.color)) for stop in ramp.elements)
    return (low, high, ramp.color_mode, ramp.hue_interpolation, ramp.interpolation, stops)


def _set_limits(ramp, scale, low, high, signed):
    old_low, old_high = (scale.inputs[key].default_value for key in ('From Min', 'From Max'))
    if signed and old_low < 0 < old_high and (old_low, old_high) != (low, high):
        # Anchor a signed ramp's zero stop. Preserve the relative spacing of
        # custom stops on each side. A freely positioned ramp stays as drawn.
        zero, new_zero = -old_low / (old_high-old_low), -low / (high-low)
        stops = [(s, s.position) for s in ramp.elements]
        if any(math.isclose(p, zero, abs_tol=1e-7) for _, p in stops):
            # Move in order away from zero to avoid passing another stop.
            for stop, p in (reversed(stops) if new_zero > zero else stops):
                stop.position = (p * new_zero / zero if p <= zero else
                                 new_zero + (p-zero) * (1-new_zero) / (1-zero))
    scale.inputs['From Min'].default_value = low
    scale.inputs['From Max'].default_value = high


def _write_scale(material, values, low, high, signed):
    src_low, src_high, mode, hue, interp, stops = values
    node, scale = scale_nodes(material)
    ramp = node.color_ramp
    ramp.color_mode, ramp.hue_interpolation, ramp.interpolation = mode, hue, interp
    while len(ramp.elements) > 1:
        ramp.elements.remove(ramp.elements[-1])
    for index, (position, color) in enumerate(stops):
        stop = ramp.elements.new(position) if index else ramp.elements[0]
        stop.position, stop.color = position, color
    scale.inputs['From Min'].default_value = src_low
    scale.inputs['From Max'].default_value = src_high
    _set_limits(ramp, scale, low, high, signed)


def copy_ramp(material, source, low, high):
    node, _ = scale_nodes(material)
    source_node, _ = scale_nodes(source)
    node['mb_anchor_zero'] = source_node.get('mb_anchor_zero', True)
    node.label = source_node.label
    signed = source.get('mb_field_kind', 'esp') != 'unsigned'
    _write_scale(material, scale_values(source), low, high, signed and node['mb_anchor_zero'])


def apply_scale(group, source, low, high, palette='current'):
    """Copy only scale settings, keeping each destination shader and its mesh."""
    settings = field_map_settings(source)
    if not settings:
        raise ValueError('Choose a mapped surface as the scale source')
    signed = settings['kind'] != 'unsigned'
    validate_limits(settings['kind'], low, high)
    slots = {}
    for obj in group['objects']:
        material = obj.active_material
        info = field_map_settings(material)
        if BASE_MESH not in obj or not info:
            raise ValueError('Choose a mapped surface')
        if settings['units'] != info['units'] or signed != (info['kind'] != 'unsigned'):
            raise ValueError('Field scales must have matching units and both be signed or both unsigned')
        scale_nodes(material)
        slots.setdefault(material, []).append(obj.material_slots[0])

    replacements = {}
    try:
        for original in slots:
            material = original.copy()
            replacements[original] = material
            material.use_fake_user = False
            material['mb_style_material'] = True
            copy_ramp(material, source, low, high)
            if palette != 'current':
                apply_palette(scale_nodes(material)[0], palette, settings['kind'], low, high)
            prefix = 'mb_esp_' if material.get('mb_field_kind', 'esp') == 'esp' else 'mb_field_'
            material[prefix + 'limits'] = (low, high)
        for original, members in slots.items():
            for slot in members:
                slot.link = 'OBJECT'
                slot.material = replacements[original]
    finally:
        for mat in replacements.values():
            if not mat.users:
                bpy.data.materials.remove(mat)
    for original in slots:
        if not original.users and original.get('mb_style_material'):
            bpy.data.materials.remove(original)


class MB_ScaleSource(bpy.types.PropertyGroup):
    object_name: StringProperty()


class MB_OT_edit_field_scale(bpy.types.Operator):
    bl_idname = 'molecular_blender.edit_field_scale'
    bl_label = 'Edit Color Scale'
    bl_description = 'Change the range and ramp, or copy them from another mapped surface'
    bl_options = {'UNDO'}

    group: StringProperty(options={'HIDDEN', 'SKIP_SAVE'})
    sources: CollectionProperty(type=MB_ScaleSource, options={'HIDDEN', 'SKIP_SAVE'})

    palette: EnumProperty(name='Color map', items=(KEEP_PALETTE, *PALETTES), options={'SKIP_SAVE'})
    copy_from: StringProperty(name='Copy scale from', options={'SKIP_SAVE'})
    minimum: FloatProperty(name='Minimum', precision=6, options={'SKIP_SAVE'})
    maximum: FloatProperty(name='Maximum', precision=6, options={'SKIP_SAVE'})

    def check(self, context):
        if self.copy_from != self._loaded_source:
            if self.copy_from:
                choice = self.sources.get(self.copy_from)
                obj = context.scene.objects.get(choice.object_name) if choice else None
                if obj is None:
                    return False
            else:
                obj = self._group(context)['objects'][0]
            self._start_preview(obj.active_material)
            self._loaded_source = self.copy_from
            return True
        kind = self._preview.get('mb_field_kind', 'esp')
        try:
            validate_limits(kind, self.minimum, self.maximum)
        except ValueError:
            return False
        node, scale = scale_nodes(self._preview)
        if self.palette != self._loaded_palette:
            self._loaded_palette = self.palette
            if self.palette != 'current':
                _set_limits(node.color_ramp, scale, self.minimum, self.maximum, False)
                apply_palette(node, self.palette, kind, self.minimum, self.maximum)
            else:
                choice = self.sources.get(self.copy_from)
                obj = context.scene.objects.get(choice.object_name) if choice else self._group(context)['objects'][0]
                if obj is None:
                    raise ValueError('The selected scale source is unavailable')
                copy_ramp(self._preview, obj.active_material, self.minimum, self.maximum)
            return True
        if (self.minimum, self.maximum) == tuple(scale.inputs[k].default_value for k in ('From Min', 'From Max')):
            return False
        _set_limits(node.color_ramp, scale, self.minimum, self.maximum, kind != 'unsigned' and node.get('mb_anchor_zero', True))
        return True

    def _start_preview(self, material):
        self.palette = self._loaded_palette = 'current'
        _, scale = scale_nodes(material)
        self.minimum = scale.inputs['From Min'].default_value
        self.maximum = scale.inputs['From Max'].default_value
        if getattr(self, '_preview', None) is None:
            self._preview = material.copy()
            self._preview.use_fake_user = False
        else:
            # Keep the node displayed by the open ramp widget alive.
            copy_ramp(self._preview, material, self.minimum, self.maximum)

    def _clear_preview(self):
        preview = getattr(self, '_preview', None)
        if preview is not None:
            bpy.data.materials.remove(preview)
            self._preview = None

    def _group(self, context):
        owner = getattr(self, '_owner', None) or resolve_owner(context.object)
        if owner:
            group = surface_groups(owner, read_record(owner)).get(self.group)
            if group and BASE_MESH in group['objects'][0] and field_map_settings(group['objects'][0].active_material):
                return group
        raise ValueError('Select a mapped surface')

    def invoke(self, context, event):
        try:
            self._owner = resolve_owner(context.object)
            group = self._group(context)
            target = group['objects'][0].active_material
            info = field_map_settings(target)
            self.sources.clear()
            for owner in context.scene.objects:
                if RECORD not in owner:
                    continue
                for key, other in surface_groups(owner, read_record(owner)).items():
                    if owner == self._owner and key == self.group:
                        continue
                    obj = other['objects'][0]
                    saved = field_map_settings(obj.active_material)
                    if (saved and saved['units'] == info['units'] and
                            (saved['kind'] == 'unsigned') == (info['kind'] == 'unsigned')):
                        try:
                            scale_nodes(obj.active_material)
                        except ValueError:
                            continue
                        item = self.sources.add()
                        item.name = f"{owner.name} · {other['label']} · {saved['label']} ({obj.name})"
                        item.object_name = obj.name
            self._loaded_source = self.copy_from = ''
            self._start_preview(target)
            return context.window_manager.invoke_props_dialog(self, width=440)
        except (ValueError, RuntimeError, ReferenceError) as exc:
            self._clear_preview()
            self.report({'ERROR'}, str(exc))
            return {'CANCELLED'}

    def draw(self, context):
        if self.sources:
            self.layout.prop_search(self, 'copy_from', self, 'sources', icon='EYEDROPPER')
        self.layout.prop(self, 'palette')
        row = self.layout.row(align=True)
        row.prop(self, 'minimum')
        row.prop(self, 'maximum')
        node, _ = scale_nodes(self._preview)
        self.layout.template_color_ramp(node, 'color_ramp', expand=True)

    def cancel(self, context):
        self._clear_preview()

    def execute(self, context):
        try:
            group = self._group(context)
            source = getattr(self, '_preview', None)
            if source is not None and self.copy_from != self._loaded_source:
                if not self.check(context):
                    raise ValueError('The selected scale source is unavailable')
                source = self._preview
            if source is None:
                # Direct API calls can supply the source object's scene name.
                if self.copy_from:
                    obj = context.scene.objects.get(self.copy_from)
                    if obj is None or BASE_MESH not in obj:
                        raise ValueError('Choose a mapped surface as the scale source')
                else:
                    obj = group['objects'][0]
                source = obj.active_material
            _, scale = scale_nodes(source)
            low = self.minimum if self.properties.is_property_set('minimum') else scale.inputs['From Min'].default_value
            high = self.maximum if self.properties.is_property_set('maximum') else scale.inputs['From Max'].default_value
            apply_scale(group, source, low, high, self.palette if getattr(self, '_preview', None) is None else 'current')
        except (ValueError, KeyError, ReferenceError, RuntimeError) as exc:
            self.report({'ERROR'}, str(exc))
            return {'CANCELLED'}
        finally:
            self._clear_preview()
        return {'FINISHED'}
