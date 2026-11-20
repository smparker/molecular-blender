"""Export the color scale driving a mapped surface as an editable SVG."""
from pathlib import Path
import xml.etree.ElementTree as ET

import bpy
from bpy.props import StringProperty
from bpy_extras.io_utils import ExportHelper

from .scene_edit import _surface_output
from .scene_esp import field_map_settings
from .scene_state import read_record, resolve_owner
from .scene_surfaces import surface_groups


def _ramp_and_limits(material):
    # Follow links so renaming nodes or adding an unrelated ramp does not change the legend.
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
                return ramp.color_ramp, low, high
    raise ValueError('Legend needs a Color Ramp connected to Base Color with a linear Map Range from 0 to 1')


def _svg_color(rgba):
    def srgb(value):
        value = min(1., max(0., value))
        return round(255 * (12.92 * value if value < .0031308 else 1.055 * value**(1/2.4) - .055))
    return '#' + ''.join(f'{srgb(v):02x}' for v in rgba[:3])


def write_field_legend(material, path):
    settings = field_map_settings(material)
    if not settings or not material.use_nodes:
        raise ValueError('Choose a surface with a field map')
    ramp, low, high = _ramp_and_limits(material)
    title = f"{settings['label']} ({settings['units']})"
    width = max(360, 40 + 13 * len(title))
    left, right = 24, width - 24
    svg = ET.Element('svg', xmlns='http://www.w3.org/2000/svg', width=str(width), height='100',
                     viewBox=f'0 0 {width} 100')
    ET.SubElement(svg, 'title').text = title
    gradient = ET.SubElement(ET.SubElement(svg, 'defs'), 'linearGradient', id='field')
    positions = {0., 1., *(float(element.position) for element in ramp.elements)}
    if ramp.interpolation != 'CONSTANT':
        positions.update(i / 256 for i in range(257))
    for p in sorted(positions):
        offset = format(p, '.9g')
        if ramp.interpolation == 'CONSTANT' and p > 0:
            # Two stops at the same point preserve the discontinuity.
            ET.SubElement(gradient, 'stop', offset=offset, attrib={'stop-color': previous})
        color = _svg_color(ramp.evaluate(p))
        ET.SubElement(gradient, 'stop', offset=offset, attrib={'stop-color': color})
        previous = color
    ET.SubElement(svg, 'rect', x=str(left), y='34', width=str(right-left), height='18', fill='url(#field)')
    labels = ET.SubElement(svg, 'g', fill='#222', attrib={'font-family': 'sans-serif', 'font-size': '14'})
    ET.SubElement(labels, 'text', x=str(left), y='23').text = title

    def tick(x, value, anchor, baseline):
        ET.SubElement(labels, 'path', d=f'M{x:g} 52 V{baseline-14}', stroke='#222')
        ET.SubElement(labels, 'text', x=f'{x:g}', y=str(baseline), attrib={'text-anchor': anchor}).text = value

    tick(left, f'{low:.6g}', 'start', 74)
    tick(right, f'{high:.6g}', 'end', 74)
    if settings['kind'] != 'unsigned' and low < 0 < high:
        x = left + (right-left) * -low / (high-low)
        # Separate the zero label vertically when it falls near an endpoint.
        tick(x, '0', 'middle', 94 if min(x-left, right-x) < 65 else 74)
    ET.ElementTree(svg).write(path, encoding='utf-8', xml_declaration=True)


class MB_OT_export_field_legend(bpy.types.Operator, ExportHelper):
    bl_idname = 'molecular_blender.export_field_legend'
    bl_label = 'Export Field Legend'
    filename_ext = '.svg'
    filter_glob: StringProperty(default='*.svg', options={'HIDDEN'})
    group: StringProperty(options={'HIDDEN', 'SKIP_SAVE'})

    def invoke(self, context, event):
        self._owner = resolve_owner(context.object)
        self.filepath = bpy.path.clean_name(self._owner.name) + '-legend.svg'
        return super().invoke(context, event)

    def execute(self, context):
        try:
            owner = getattr(self, '_owner', None) or resolve_owner(context.object)
            if owner is None:
                raise ValueError('Select a molecule with a mapped surface')
            objects = surface_groups(owner, read_record(owner))[self.group]['objects']
            obj = context.object if context.object in objects else objects[0]
            path = Path(bpy.path.abspath(self.filepath))
            write_field_legend(obj.active_material, path)
        except (ValueError, KeyError, OSError, ReferenceError) as exc:
            self.report({'ERROR'}, str(exc))
            return {'CANCELLED'}
        self.report({'INFO'}, f'Exported {path.name}')
        return {'FINISHED'}
