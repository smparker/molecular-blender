"""Molecule-scoped appearance edits; scientific fields are never recomputed."""
import bpy
from bpy.props import FloatProperty, EnumProperty

from .options import validated_options
from .periodictable import elements
from .scene_state import (read_record, write_record, resolve_owner, owned_objects,
                          normalize_duplicate_ids, ID)
from .stylers import get_styler
from .util import darken_color


def _atom_radius(element, settings):
    if settings['plot_style'] == 'sticks':
        return settings['bond_thickness']
    return elements[element].vdw * (1.0 if settings['plot_style'] == 'vdw' else 0.25)


def _private_curve(obj):
    if obj.data.users > 1:
        obj.data = obj.data.copy()


def update_junctions(owner, record):
    """Scale objects, preserving manual nonuniform scaling and linked atom meshes."""
    from .plotter import _multiple_bond_profile
    settings = record['settings']
    extent = {}
    for obj in owned_objects(owner, 'bond'):
        radius = obj.data.bevel_depth
        order = obj.get('mb_order', 1)
        bond_extent = _multiple_bond_profile(radius, order)[2] if order > 1 else radius
        for endpoint in obj.get('mb_endpoints', []):
            extent[endpoint] = max(extent.get(endpoint, 0), bond_extent)
    atoms = owned_objects(owner, 'atom')
    bpy.context.view_layer.update()
    old_matrices = {obj: obj.matrix_world.copy() for obj in atoms}
    for obj in atoms:
        base_radius = obj.get('mb_atom_radius', _atom_radius(obj['mb_element'], settings))
        if base_radius <= 0:
            raise ValueError('Cannot resize zero-radius atom geometry; reimport with a positive atom radius')
        obj['mb_atom_radius'] = base_radius
        radius = _atom_radius(obj['mb_element'], settings)
        if settings['plot_style'] == 'sticks':
            radius = max(radius, extent.get(obj[ID], 0))
        target = settings.get('atom_scale', 1.0) * radius / base_radius
        previous = obj.get('mb_style_scale', obj.get('mb_base_scale', [1])[0])
        if previous <= 0:
            raise ValueError('Cannot resize an atom with zero saved scale')
        obj.scale *= target / previous
        obj['mb_style_scale'] = target
    bpy.context.view_layer.update()
    # Hooks must follow atom movement, not changes to the display sphere's size.
    # Preserve their current deformation while compensating for style scaling.
    for bond in owned_objects(owner, 'bond'):
        for hook in bond.modifiers:
            if hook.type == 'HOOK' and hook.object in old_matrices:
                atom = hook.object
                hook.matrix_inverse = (atom.matrix_world.inverted() @ old_matrices[atom]
                                       @ hook.matrix_inverse)


def _set_color(material, color):
    material.diffuse_color = color
    if material.use_nodes:
        for node in material.node_tree.nodes:
            if node.type == 'BSDF_PRINCIPLED':
                node.inputs['Base Color'].default_value = color
                for name in ('Emission Color', 'Emission'):
                    if name in node.inputs:
                        node.inputs[name].default_value = darken_color(color, 0.3)
            elif node.type == 'BSDF_DIFFUSE':
                node.inputs['Color'].default_value = color


def _set_opacity(material, opacity):
    """Wrap the existing shader, retaining its original surface appearance."""
    material.use_nodes = True
    tree = material.node_tree
    mix = tree.nodes.get('MB Surface Opacity')
    if mix is None:
        output = next((n for n in tree.nodes if n.type == 'OUTPUT_MATERIAL' and n.is_active_output), None)
        if output is None or not output.inputs['Surface'].links:
            raise ValueError('Surface material has no connected surface shader')
        original = output.inputs['Surface'].links[0].from_socket
        transparent = tree.nodes.new('ShaderNodeBsdfTransparent')
        mix = tree.nodes.new('ShaderNodeMixShader')
        mix.name = 'MB Surface Opacity'
        mix.label = 'Molecule surface opacity'
        tree.links.new(transparent.outputs[0], mix.inputs[1])
        tree.links.new(original, mix.inputs[2])
        tree.links.new(mix.outputs[0], output.inputs['Surface'])
    mix.inputs[0].default_value = opacity
    if hasattr(material, 'surface_render_method'):
        material.surface_render_method = 'DITHERED'
    elif hasattr(material, 'blend_method'):
        material.blend_method = 'BLEND'


def restyle_molecule(owner, *, atom_scale=None, bond_thickness=None,
                     colors=None, surface_opacity=None):
    normalize_duplicate_ids()
    record = read_record(owner)
    old_settings = record['settings']
    atom_scale = old_settings['atom_scale'] if atom_scale is None else atom_scale
    bond_thickness = old_settings['bond_thickness'] if bond_thickness is None else bond_thickness
    colors = old_settings['colors'] if colors is None else colors
    surface_opacity = old_settings['surface_opacity'] if surface_opacity is None else surface_opacity
    settings = validated_options(dict(old_settings, atom_scale=atom_scale,
                                      bond_thickness=bond_thickness, colors=colors,
                                      surface_opacity=surface_opacity))
    if bond_thickness <= 0:
        raise ValueError('Bond thickness must be positive when restyling')
    atoms = owned_objects(owner, 'atom')
    bonds = owned_objects(owner, 'bond')
    surfaces = owned_objects(owner, 'surface')
    # Fail before making changes if a saved/externally edited object cannot be handled.
    for atom in atoms:
        if atom.get('mb_atom_radius', _atom_radius(atom['mb_element'], old_settings)) <= 0:
            raise ValueError('Cannot resize zero-radius atom geometry; reimport with a positive atom radius')
    for obj in bonds:
        if obj.type != 'CURVE' or len(obj.data.splines) != 1 or len(obj.data.splines[0].bezier_points) != 2:
            raise ValueError('Bond editing requires straight two-point curves')
        points = obj.data.splines[0].bezier_points
        if (points[1].co - points[0].co).length < 1e-8:
            raise ValueError('Cannot restyle a zero-length bond')
    if settings['object_type'] == 'wireframe':
        raise ValueError('Restyling individual atoms and bonds requires mesh, NURBS, or metaball imports')

    # Prepare private materials before touching scene links. Object material slots
    # keep linked atom geometry independent across copied/imported molecules.
    copies = {}
    assignments = []
    styler = get_styler(settings)
    try:
        for obj in atoms + bonds + surfaces:
            is_surface = obj in surfaces
            symbol = obj.get('mb_element')
            for slot in obj.material_slots:
                material = slot.material
                if material is None or (not is_surface and symbol is None):
                    continue
                key = (material.as_pointer(), symbol, is_surface)
                if key not in copies:
                    copy = material.copy()
                    copies[key] = copy
                    copy['mb_style_material'] = True
                    if is_surface:
                        _set_opacity(copy, surface_opacity)
                    else:
                        _set_color(copy, styler.element_color(elements[symbol]))
                assignments.append((slot, material, copies[key]))
    except Exception:
        for material in copies.values():
            bpy.data.materials.remove(material)
        raise

    from .plotter import set_bond_order
    bond_records = {b['id']: b for b in record['bonds']}
    for obj in bonds:
        _private_curve(obj)
        # Preserve element-dependent relative radii, even after repeated edits.
        initial = obj.get('mb_radius', obj.data.bevel_depth)
        nominal = obj.get('mb_nominal_radius', old_settings['bond_thickness'])
        obj['mb_nominal_radius'] = nominal
        obj['mb_radius'] = initial
        obj.data.bevel_depth = bond_thickness * initial / nominal if nominal else bond_thickness
        tilts = [p.tilt for p in obj.data.splines[0].bezier_points]
        set_bond_order(obj, obj.get('mb_order', 1), bond_records[obj[ID]]['direction'])
        # Preserve manual curve tilt when there was no recorded direction.
        if bond_records[obj[ID]]['direction'] is None:
            for point, tilt in zip(obj.data.splines[0].bezier_points, tilts):
                point.tilt = tilt
    for atom in atoms:
        atom['mb_atom_radius'] = atom.get('mb_atom_radius', _atom_radius(atom['mb_element'], old_settings))
    record['settings'] = settings
    update_junctions(owner, record)
    for slot, old, new in assignments:
        slot.link = 'OBJECT'
        slot.material = new
    write_record(owner, record)
    for old in {old for _, old, _ in assignments}:
        if old.users == 0 and old.get('mb_style_material'):
            bpy.data.materials.remove(old)


class MB_OT_restyle(bpy.types.Operator):
    bl_idname = 'molecular_blender.restyle'
    bl_label = 'Restyle molecule'
    bl_description = 'Change the selected molecule without rereading the source or recomputing fields'
    bl_options = {'REGISTER', 'UNDO'}

    atom_scale: FloatProperty(name='Atom size multiplier', default=1.0, min=0.01, max=100)
    bond_thickness: FloatProperty(name='Bond thickness', default=0.2, min=0.001, max=50)
    colors: EnumProperty(name='Palette', items=(('default', 'Default', ''), ('vmd', 'VMD', '')))
    surface_opacity: FloatProperty(name='Surface opacity', default=1.0, min=0, max=1)

    @classmethod
    def poll(cls, context):
        return context.mode == 'OBJECT' and resolve_owner(context.object) is not None

    def invoke(self, context, event):
        settings = read_record(resolve_owner(context.object))['settings']
        for name in ('atom_scale', 'bond_thickness', 'colors', 'surface_opacity'):
            setattr(self, name, settings[name])
        return context.window_manager.invoke_props_dialog(self)

    def execute(self, context):
        try:
            restyle_molecule(resolve_owner(context.object), atom_scale=self.atom_scale,
                             bond_thickness=self.bond_thickness, colors=self.colors,
                             surface_opacity=self.surface_opacity)
        except (ValueError, KeyError) as exc:
            self.report({'ERROR'}, str(exc))
            return {'CANCELLED'}
        return {'FINISHED'}


def register():
    bpy.utils.register_class(MB_OT_restyle)


def unregister():
    bpy.utils.unregister_class(MB_OT_restyle)
