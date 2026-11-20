"""Save Molden normal coordinates and preview them on Blender's timeline."""
from functools import lru_cache
import math

import bpy
from bpy.props import EnumProperty, FloatProperty, IntProperty

from .scene_state import ID, owned_objects, read_record, resolve_owner

DATA = 'mb_vibrations'
BASE = 'mb_vibration_base'
ACTIVE = 'mb_preview_mode'


def store_modes(owner, data):
    def flatten(vectors):
        return [v for xyz in vectors for v in xyz]

    owner[DATA] = dict(frequencies=data['frequencies'], equilibrium=flatten(data['equilibrium']),
                       vectors={str(i): flatten(mode) for i, mode in enumerate(data['displacements'], 1)})
    owner['mb_vibration_mode'] = next((i for i, f in enumerate(data['frequencies']) if abs(f) > 1.e-8), 0)


@lru_cache(maxsize=None)
def _frequency_items(frequencies):
    # Blender requires Python to retain the strings used by dynamic enum items.
    return [(str(i+1), f'{i+1}   ·   {freq:g} cm⁻¹', '', i) for i, freq in enumerate(frequencies)]


def _modes(owner, context):
    return _frequency_items(tuple(owner[DATA]['frequencies'])) if DATA in owner else []


def _has_delta_animation(obj):
    anim = obj.animation_data
    if anim is None:
        return False
    if any(c.data_path == 'delta_location' for c in anim.drivers):
        return True
    actions = [anim.action]
    def strip_actions(strips):
        for strip in strips:
            if strip.action:
                yield strip.action
            yield from strip_actions(strip.strips)
    for track in anim.nla_tracks:
        actions.extend(strip_actions(track.strips))
    return any(action and any(c.data_path == 'delta_location' for c in action.fcurves)
               for action in actions)


def clear_mode(owner):
    for obj in owned_objects(owner, 'atom'):
        if BASE in obj:
            for axis in range(3):
                obj.driver_remove('delta_location', axis)
            obj.delta_location = obj[BASE]
            del obj[BASE]
    if ACTIVE in owner:
        del owner[ACTIVE]


def apply_mode(owner):
    if DATA not in owner:
        raise ValueError('Import a single Molden geometry containing normal coordinates')
    data = owner[DATA]
    mode = owner.mb_vibration_mode
    vectors = data['vectors'][mode]
    peak = max(math.hypot(*vectors[i:i+3]) for i in range(0, len(vectors), 3))
    if peak == 0:
        raise ValueError('This mode has zero displacement for every atom')
    atoms = owned_objects(owner, 'atom')
    if not atoms:
        raise ValueError('This preview needs displayed atom objects')
    if any(BASE not in obj and _has_delta_animation(obj) for obj in atoms):
        raise ValueError('An atom already has animated Delta Location channels')
    indices = {atom['id']: (i, atom['position']) for i, atom in enumerate(read_record(owner)['atoms'])}
    scale = owner.mb_mode_amplitude / peak
    speed = math.tau / owner.mb_mode_period
    phase = owner.mb_mode_phase
    equilibrium = data['equilibrium']
    try:
        for obj in atoms:
            index, original = indices[obj[ID]]
            if BASE not in obj:
                obj[BASE] = list(obj.delta_location)
            for axis in range(3):
                j = 3 * index + axis
                center = obj[BASE][axis] + equilibrium[j] - original[axis]
                amplitude = scale * vectors[j]
                driver = obj.driver_add('delta_location', axis).driver
                driver.type = 'SCRIPTED'
                driver.expression = f'{center:.12g}+({amplitude:.12g})*sin((frame-1)*{speed:.12g}+({phase:.12g}))'
    except Exception:
        clear_mode(owner)
        raise
    owner[ACTIVE] = int(mode)


class MB_OT_preview_mode(bpy.types.Operator):
    bl_idname = 'molecular_blender.preview_mode'
    bl_label = 'Preview Mode'
    bl_description = 'Apply the selected normal mode to the timeline as a repeating oscillation'
    bl_options = {'UNDO'}

    def execute(self, context):
        owner = resolve_owner(context.object)
        if owner is None:
            return {'CANCELLED'}
        try:
            apply_mode(owner)
        except (ValueError, RuntimeError, KeyError) as exc:
            self.report({'ERROR'}, str(exc))
            return {'CANCELLED'}
        context.scene.frame_set(context.scene.frame_current, subframe=context.scene.frame_subframe)
        return {'FINISHED'}


class MB_OT_clear_mode(bpy.types.Operator):
    bl_idname = 'molecular_blender.clear_mode'
    bl_label = 'Restore Geometry'
    bl_description = 'Remove the mode preview and restore the atom offsets'
    bl_options = {'UNDO'}

    def execute(self, context):
        owner = resolve_owner(context.object)
        if owner is None or ACTIVE not in owner:
            return {'CANCELLED'}
        clear_mode(owner)
        context.view_layer.update()
        return {'FINISHED'}


def draw_modes(layout, owner, context):
    box = layout.box()
    box.label(text='Vibrational modes')
    box.prop(owner, 'mb_vibration_mode')
    box.prop(owner, 'mb_mode_amplitude')
    box.prop(owner, 'mb_mode_period')
    box.prop(owner, 'mb_mode_phase')
    box.operator('molecular_blender.preview_mode', text='Update Preview' if ACTIVE in owner else 'Preview Mode')
    if ACTIVE in owner:
        row = box.row(align=True)
        playing = context.screen.is_animation_playing
        row.operator('screen.animation_play', text='Pause' if playing else 'Play',
                     icon='PAUSE' if playing else 'PLAY')
        row.operator('molecular_blender.clear_mode')


def register_properties():
    bpy.types.Object.mb_vibration_mode = EnumProperty(name='Mode', items=_modes)
    bpy.types.Object.mb_mode_amplitude = FloatProperty(name='Peak displacement (Å)', min=0., default=.3,
        description='Maximum displacement of any atom from mode equilibrium')
    bpy.types.Object.mb_mode_period = IntProperty(name='Frames per cycle', min=2, default=24)
    bpy.types.Object.mb_mode_phase = FloatProperty(name='Starting phase', subtype='ANGLE', default=0.)


def unregister_properties():
    for key in ('mb_mode_phase', 'mb_mode_period', 'mb_mode_amplitude', 'mb_vibration_mode'):
        delattr(bpy.types.Object, key)
    _frequency_items.cache_clear()
