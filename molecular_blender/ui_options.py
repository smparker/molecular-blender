"""Translate a serializable option declaration to a Blender RNA property."""
from bpy.props import BoolProperty, EnumProperty, FloatProperty, IntProperty, StringProperty

from .options import OPTIONS


def option_property(key, *, api_default=False, **overrides):
    spec = OPTIONS[key]
    factory = {bool: BoolProperty, float: FloatProperty, int: IntProperty, str: StringProperty}
    kwargs = dict(spec.ui, **overrides)
    if api_default or 'default' not in kwargs:
        kwargs['default'] = spec.default
    prop = EnumProperty if spec.choices else factory[type(kwargs['default'])]
    return prop(**kwargs)
