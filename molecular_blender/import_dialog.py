"""File-browser previews and surface controls for the molecular import operator."""
from functools import lru_cache
import json
import os
import re
import stat

import bpy

from .importers import cube_dataset_labels, file_format

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
        return 'cube', None
    child = os.path.join(os.path.dirname(path), first['filename'])
    return file_format(child), child


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
        kind, path = _preview(_manifest_source, path, first_frame) or (None, None)
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
def _dataset_item(index, label):
    # Blender keeps pointers to the strings returned by dynamic enum callbacks.
    name = 'Dataset {}'.format(index + 1)
    if label is not None:
        name += ' (ID {})'.format(label)
    return str(index), name, '', index


def dataset_items(op, context):
    _, labels = surface_source(browser_path(op, context), op.frame_start)
    if labels is None or len(labels) > MAX_DATASET_MENU:
        return [_dataset_item(0, None)]
    has_ids = not isinstance(labels, range)
    return [_dataset_item(i, label if has_ids else None) for i, label in enumerate(labels)]


def dataset_index(index, labels):
    if labels is not None and len(labels) <= MAX_DATASET_MENU and index >= len(labels):
        return 0
    return index


def get_dataset_choice(op):
    _, labels = surface_source(browser_path(op, bpy.context), op.frame_start)
    return dataset_index(op.cube_dataset, labels)


def draw_surfaces(layout, op, context):
    kind, labels = surface_source(browser_path(op, context), op.frame_start)
    if kind not in ('molden', 'cube'):
        return
    box = layout.box()
    box.label(text='Surfaces')
    draw_surface_controls(box, op, kind, labels)


def draw_surface_controls(box, op, kind, labels=None):
    if kind == 'molden':
        box.prop(op, 'volume')
        if op.volume == 'orbital':
            box.prop(op, 'orbital_choice')
            if op.orbital_choice == 'custom':
                box.prop(op, 'orbital', text='Selection')
            if uses_frontier_spin(op.orbital):
                box.prop(op, 'orbital_spin')
    elif labels is not None and len(labels) > 1:
        if len(labels) <= MAX_DATASET_MENU:
            box.prop(op, 'dataset_choice')
        else:
            box.prop(op, 'cube_dataset', text='Dataset index (0–{})'.format(len(labels) - 1))
            if not isinstance(labels, range) and op.cube_dataset < len(labels):
                box.label(text='Dataset ID {}'.format(labels[op.cube_dataset]))

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
            box.prop(op, 'orca_fix')
        else:
            box.prop(op, 'cube_both_signs')
        box.prop(op, 'remesh')
