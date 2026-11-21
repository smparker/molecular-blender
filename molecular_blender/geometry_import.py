"""Track Blender data allocated while building one molecular geometry."""
from contextlib import contextmanager
from contextvars import ContextVar

import bpy


_created = ContextVar('molecular_blender_geometry_ids', default=None)


def track_new(block):
    """Register a newly allocated ID, including objects not yet linked to a collection."""
    blocks = _created.get()
    if blocks is not None:
        blocks.add(block)
    return block


def _discard_ids(blocks):
    # Conversion and profile replacement can already have freed an ID.
    live = [block for block in blocks if block.as_pointer()]
    # Atom transforms and animated bond curves both allocate actions implicitly.
    actions = {obj.animation_data.action for obj in live if isinstance(obj, (bpy.types.Object, bpy.types.Curve))
               and obj.animation_data and obj.animation_data.action}
    if live:
        bpy.data.batch_remove(live)
    for action in actions:
        if not action.users:
            bpy.data.actions.remove(action)


@contextmanager
def rollback_new_ids():
    """Remove only this operation's new Blender data if publication fails."""
    blocks = set()
    token = _created.set(blocks)
    try:
        yield
    except Exception:
        _discard_ids(blocks)
        raise
    finally:
        _created.reset(token)


@contextmanager
def geometry_import(context):
    layer = context.view_layer
    active = layer.objects.active
    mode = active.mode if active else 'OBJECT'
    selection = [obj for obj in layer.objects if obj.select_get(view_layer=layer)]
    blocks = set()
    token = _created.set(blocks)
    try:
        if mode != 'OBJECT':
            bpy.ops.object.mode_set(mode='OBJECT')
        yield
    except Exception:
        if layer.objects.active and layer.objects.active.mode != 'OBJECT':
            bpy.ops.object.mode_set(mode='OBJECT')
        _discard_ids(blocks)
        layer.update()
        for obj in layer.objects:
            if obj.select_get(view_layer=layer):
                obj.select_set(False, view_layer=layer)
        for obj in selection:
            obj.select_set(True, view_layer=layer)
        layer.objects.active = active
        if mode != 'OBJECT':
            bpy.ops.object.mode_set(mode=mode)
        raise
    finally:
        _created.reset(token)
