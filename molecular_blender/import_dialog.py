"""File-browser previews and surface controls for the molecular import operator."""
from functools import lru_cache
import json
import os
import re
import stat

import bpy

from .importers import cube_dataset_labels, file_format, preview_molden_convention
from .options import MOLDEN_CONVENTION_LABELS
from .turbomole_reader import response_headers

MAX_DATASET_MENU = 256


def browser_path(operator, context):
    space = getattr(context, 'space_data', None)
    if getattr(space, 'type', None) != 'FILE_BROWSER' or not space.params:
        return operator.filepath
    params = space.params
    entry = getattr(context, 'active_file', None)
    name = params.filename
    if entry:
        name = getattr(entry, 'relative_path', None) or entry.name
    path = os.path.join(os.fsdecode(params.directory), name)
    return path if name and not os.path.isdir(path) else ''


def _manifest_source(path, first_frame):
    with open(path) as stream:
        entries = json.load(stream)['molecules']
    first = entries[0 if len(entries) == 1 else first_frame]
    if 'linear_combination' in first:
        return 'cube', None, len(entries), None
    child = os.path.join(os.path.dirname(path), first['filename'])
    polar = first.get('polar')
    return file_format(child), child, len(entries), os.path.join(os.path.dirname(path), polar) if polar else None


@lru_cache(maxsize=16)
def _read_preview(reader, path, mtime, size, args):
    # Cache failures too. The importer reports the full error on confirmation.
    try:
        return reader(path, *args)
    except (OSError, ValueError, EOFError, IndexError, KeyError, TypeError, OverflowError):
        return None


def _preview(reader, path, *args):
    try:
        info = os.stat(path)
    except (OSError, ValueError):
        return None
    if stat.S_ISREG(info.st_mode):
        return _read_preview(reader, path, info.st_mtime_ns, info.st_size, args)
    return None


def surface_source(path, first_frame=0):
    """The first frame supplies a JSON trajectory's field, just as in the reader."""
    kind = file_format(path)
    if kind == 'json':
        kind, path, _, _ = _preview(_manifest_source, path, first_frame) or (None, None, 0, None)
        if kind == 'cube' and path is None:
            return 'cube', range(1)
    labels = _preview(cube_dataset_labels, path) if kind == 'cube' else None
    return kind, labels


def uses_frontier_spin(expression):
    # Absolute file indices and channel-qualified names ignore the global channel.
    return any(re.fullmatch(r'(?:homo|lumo)(?:\s*[-+]\s*\d+)?|[+-]?0+|-\d+', token.strip(), re.I)
               for token in expression.split(','))


def get_orbital_choice(op):
    if op.get('editing_custom_orbital', False):
        return 2
    return {'homo': 0, 'lumo': 1}.get(op.orbital.strip().lower(), 2)


def set_orbital_choice(op, value):
    op['editing_custom_orbital'] = value == 2
    if value < 2:
        if op.orbital.strip().lower() not in ('homo', 'lumo'):
            op['last_custom_orbital'] = op.orbital
        op.orbital = ('homo', 'lumo')[value]
    else:
        op.orbital = op.get('last_custom_orbital', op.orbital)


@lru_cache(maxsize=None)
def dataset_item(index, label):
    # Blender keeps pointers to the strings returned by dynamic enum callbacks.
    name = 'Dataset {}'.format(index + 1)
    if label is not None:
        name += ' (ID {})'.format(label)
    return str(index), name, '', index


def dataset_items(op, context):
    _, labels = surface_source(browser_path(op, context), op.frame_start)
    if labels is None or len(labels) > MAX_DATASET_MENU:
        return [dataset_item(0, None)]
    has_ids = not isinstance(labels, range)
    return [dataset_item(i, label if has_ids else None) for i, label in enumerate(labels)]


def dataset_index(index, labels):
    if labels is not None and len(labels) <= MAX_DATASET_MENU and index >= len(labels):
        return 0
    return index


def get_dataset_choice(op):
    _, labels = surface_source(browser_path(op, bpy.context), op.frame_start)
    return dataset_index(op.cube_dataset, labels)


def draw_surfaces(layout, op, context):
    path = browser_path(op, context)
    kind, labels = surface_source(path, op.frame_start)
    if kind not in ('molden', 'cube'):
        return
    box = layout.box()
    box.label(text='Surfaces')
    trajectory = False
    responses = False
    if kind == 'molden' and file_format(path) == 'json':
        _, path, count, polar = _preview(_manifest_source, path, op.frame_start)
        responses = bool(polar)
        last = min(count - 1, op.frame_end) if op.frame_end >= 0 else count - 1
        trajectory = last - op.frame_start >= op.frame_step
    draw_surface_controls(box, op, kind, labels, path, trajectory=trajectory, responses=responses)


def draw_surface_controls(box, op, kind, labels=None, path=None, *, trajectory=False, responses=False):
    if kind == 'molden':
        box.prop(op, 'field_with_response' if responses else 'volume')
        if responses and op.field_with_response == 'response':
            box.prop(op, 'response_choice')
        elif op.volume == 'orbital':
            box.prop(op, 'orbital_choice')
            if op.orbital_choice == 'custom':
                box.prop(op, 'orbital', text='Selection')
            if uses_frontier_spin(op.orbital):
                box.prop(op, 'orbital_spin')
            if trajectory:
                box.prop(op, 'align_orbital_phase')
    elif labels is not None and len(labels) > 1:
        if len(labels) <= MAX_DATASET_MENU:
            box.prop(op, 'dataset_choice')
        else:
            box.prop(op, 'cube_dataset', text='Dataset index (0–{})'.format(len(labels) - 1))
            if not isinstance(labels, range) and op.cube_dataset < len(labels):
                box.label(text='Dataset ID {}'.format(labels[op.cube_dataset]))

    draw_sampling_controls(box, op, kind, path)


def draw_sampling_controls(box, op, kind, path=None):
    box.prop(op, 'level_mode')
    box.prop(op, 'isovalues', text='Fractions' if op.cumulative else 'Isovalues')
    if kind == 'cube' and op.cumulative:
        box.prop(op, 'cube_field_type')
    if kind == 'molden':
        box.prop(op, 'resolution', text='Grid spacing (Å)')
    box.prop(op, 'surface_advanced', icon='TRIA_DOWN' if op.surface_advanced else 'TRIA_RIGHT',
             emboss=False)
    if op.surface_advanced:
        if kind == 'molden':
            box.prop(op, 'molden_convention')
            if op.molden_convention == 'auto' and path:
                convention = _preview(preview_molden_convention, path)
                label = MOLDEN_CONVENTION_LABELS.get(convention)
                box.label(text='Auto uses ' + label if label else 'Auto uses standard Molden unless the title identifies ORCA')
        else:
            box.prop(op, 'cube_both_signs')
        box.prop(op, 'remesh')


def response_path(path, start):
    if file_format(path) == 'json':
        manifest = _preview(_manifest_source, path, start)
        if manifest and manifest[0] == 'molden':
            return manifest[3]


@lru_cache(maxsize=None)
def pair_item(index, pair, value):
    return str(pair), f'Eigenpair {pair}  ·  {value:.7g}', 'DIPL eigenpair and eigenvalue', index


def response_items(op, context):
    path = response_path(browser_path(op, context), op.frame_start)
    pairs = _preview(response_headers, path) if path else ()
    return [pair_item(i, pair, value) for i, (pair, value) in enumerate(pairs or ())] or [
        ('0', 'First eigenpair', 'Read the response file at import', 0)]


def saved_response_pairs(record):
    if 'response_pairs' in record:
        return record['response_pairs']
    # Older blends did not store the DIPL headers. Preview just the imported frame.
    if record['orbitals'] and file_format(record['source']) == 'json':
        settings = record['settings']
        start, end, step = (settings[k] for k in ('frame_start', 'frame_end', 'frame_step'))
        manifest = _preview(_manifest_source, record['source'], start)
        if manifest and manifest[0] == 'molden' and manifest[3]:
            last = manifest[2] - 1 if end < 0 else min(end, manifest[2] - 1)
            if last - start < step:
                return _preview(response_headers, manifest[3]) or ()
    return ()


def get_response_choice(op):
    items = response_items(op, bpy.context)
    return next((item[3] for item in items if item[0] == str(op.response_eigenpair)), 0)


def set_response_choice(op, index):
    items = response_items(op, bpy.context)
    if 0 <= index < len(items):
        op.response_eigenpair = int(items[index][0])
