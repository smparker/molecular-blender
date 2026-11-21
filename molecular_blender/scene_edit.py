"""Molecule-scoped appearance edits; scientific fields are never recomputed."""
import bpy
from bpy.props import BoolProperty, FloatProperty, FloatVectorProperty, EnumProperty

from .geometry_import import track_new
from .options import validated_options
from .periodictable import elements
from .scene_state import (atom_bond_style, read_record, write_record, resolve_owner, owned_objects,
                          normalize_duplicate_ids, ID, ROLE)
from .stylers import get_styler
from .ui_options import option_property
from .util import darken_color


def _atom_radius(element, settings):
    if settings['plot_style'] == 'sticks':
        return settings['bond_thickness']
    return elements[element].vdw * (1.0 if settings['plot_style'] == 'vdw' else 0.25)


def _private_curve(obj):
    if obj.data.users > 1:
        obj.data = track_new(obj.data.copy())


def apply_bond_profile(obj, bond, order):
    """Rebuild a profile with its saved orientation, including older scene records."""
    from .scene_geometry import set_bond_order
    if not obj.data.get('mb_trajectory_tilt') and 'tilt' not in bond and bond['direction'] is not None:
        set_bond_order(obj, order, bond['direction'])
        return
    # Older records with no saved orientation can have native Blender tilt edits.
    points = obj.data.splines[0].bezier_points
    tilts = [bond.get('tilt', point.tilt) for point in points]
    set_bond_order(obj, order)
    for point, tilt in zip(points, tilts):
        point.tilt = tilt


def space_filling_bond(obj, enabled):
    """Hide a bond for the atom envelope, remembering its own visibility."""
    key = 'mb_space_filling_visibility'
    if enabled and key not in obj:
        obj[key] = [int(obj.hide_viewport), int(obj.hide_render)]
        obj.hide_viewport = obj.hide_render = True
    elif not enabled and key in obj:
        obj.hide_viewport, obj.hide_render = map(bool, obj[key])
        del obj[key]


def update_junctions(owner, record):
    """Scale objects, preserving manual nonuniform scaling and linked atom meshes."""
    from .scene_geometry import _multiple_bond_profile
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
            elif node.type in ('BSDF_DIFFUSE', 'EMISSION'):
                node.inputs['Color'].default_value = color


def _surface_output(material):
    output = next((n for n in material.node_tree.nodes if n.type == 'OUTPUT_MATERIAL' and n.is_active_output), None)
    if output is None or not output.inputs['Surface'].links:
        raise ValueError('Surface material has no connected surface shader')
    return output


def _set_opacity(material, opacity):
    """Wrap the existing shader, retaining its original surface appearance."""
    material.use_nodes = True
    tree = material.node_tree
    mix = tree.nodes.get('MB Surface Opacity')
    if mix is None:
        output = _surface_output(material)
        original = output.inputs['Surface'].links[0].from_socket
        transparent = tree.nodes.new('ShaderNodeBsdfTransparent')
        mix = tree.nodes.new('ShaderNodeMixShader')
        mix.name = 'MB Surface Opacity'
        mix.label = 'Molecule surface opacity'
        tree.links.new(transparent.outputs[0], mix.inputs[1])
        tree.links.new(original, mix.inputs[2])
        tree.links.new(mix.outputs[0], output.inputs['Surface'])
    mix.inputs[0].default_value = opacity
    material.diffuse_color[3] = opacity
    if hasattr(material, 'surface_render_method'):
        material.surface_render_method = 'DITHERED'
    elif hasattr(material, 'blend_method'):
        material.blend_method = 'BLEND'


def restyle_molecule(owner, *, atom_scale=None, bond_thickness=None,
                     colors=None, surface_opacity=None, reset_style=False):
    normalize_duplicate_ids()
    record = read_record(owner)
    old_settings = record['settings']
    # Older scene records only retain their most recently saved settings.
    initial_style = record.setdefault('initial_style', atom_bond_style(old_settings))
    if reset_style:
        atom_scale = initial_style['atom_scale']
        bond_thickness = initial_style['bond_thickness']
        colors = initial_style['colors']
    atom_scale = old_settings['atom_scale'] if atom_scale is None else atom_scale
    bond_thickness = old_settings['bond_thickness'] if bond_thickness is None else bond_thickness
    colors = old_settings['colors'] if colors is None else colors
    edit_surface_opacity = surface_opacity is not None
    surface_opacity = old_settings['surface_opacity'] if surface_opacity is None else surface_opacity
    settings = validated_options(dict(old_settings, atom_scale=atom_scale,
                                      bond_thickness=bond_thickness, colors=colors,
                                      surface_opacity=surface_opacity))
    atoms = owned_objects(owner, 'atom')
    bonds = owned_objects(owner, 'bond')
    if bond_thickness <= 0 and (bonds or settings['plot_style'] == 'sticks'):
        raise ValueError('Bond thickness must be positive when restyling bonds')
    surfaces = owned_objects(owner, 'surface') if edit_surface_opacity else []
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
            accented = 'mb_accent_color' in obj
            if not is_surface and (symbol is None or accented) and not (reset_style and accented):
                continue
            for slot in obj.material_slots:
                material = slot.material
                if material is None:
                    continue
                key = (material.as_pointer(), symbol, is_surface)
                if key not in copies:
                    copy = track_new(material.copy())
                    copies[key] = copy
                    copy['mb_style_material'] = True
                    if is_surface:
                        _set_opacity(copy, surface_opacity)
                    else:
                        _set_color(copy, styler.element_color(elements[symbol]) if symbol else styler.bond_color)
                assignments.append((slot, material, copies[key]))
    except Exception:
        for material in copies.values():
            bpy.data.materials.remove(material)
        raise

    bond_records = {b['id']: b for b in record['bonds']}
    for obj in bonds:
        _private_curve(obj)
        # Preserve element-dependent relative radii, even after repeated edits.
        initial = obj.get('mb_radius', obj.data.bevel_depth)
        nominal = obj.get('mb_nominal_radius', old_settings['bond_thickness'])
        obj['mb_nominal_radius'] = nominal
        obj['mb_radius'] = initial
        obj.data.bevel_depth = bond_thickness * initial / nominal if nominal else bond_thickness
        apply_bond_profile(obj, bond_records[obj[ID]], obj.get('mb_order', 1))
    for atom in atoms:
        atom['mb_atom_radius'] = atom.get('mb_atom_radius', _atom_radius(atom['mb_element'], old_settings))
    record['settings'] = settings
    update_junctions(owner, record)
    for slot, old, new in assignments:
        slot.link = 'OBJECT'
        slot.material = new
    if reset_style:
        for obj in atoms + bonds:
            if 'mb_accent_color' in obj:
                del obj['mb_accent_color']
    write_record(owner, record)
    for old in {old for _, old, _ in assignments}:
        if old.users == 0 and old.get('mb_style_material'):
            bpy.data.materials.remove(old)


class MB_OT_reset_atom_bond_style(bpy.types.Operator):
    bl_idname = 'molecular_blender.reset_atom_bond_style'
    bl_label = 'Reset Atom and Bond Style'
    bl_description = 'Restore the saved starting atom size, bond thickness and palette, and clear fragment accents'
    bl_options = {'REGISTER', 'UNDO'}

    @classmethod
    def poll(cls, context):
        owner = resolve_owner(context.object)
        return context.mode == 'OBJECT' and owner is not None and read_record(owner)['settings']['object_type'] != 'wireframe'

    def execute(self, context):
        try:
            restyle_molecule(resolve_owner(context.object), reset_style=True)
        except (ValueError, RuntimeError) as exc:
            self.report({'ERROR'}, str(exc))
            return {'CANCELLED'}
        return {'FINISHED'}


class MB_OT_restyle(bpy.types.Operator):
    bl_idname = 'molecular_blender.restyle'
    bl_label = 'Restyle molecule'
    bl_description = 'Change the selected molecule without rereading the source or recomputing fields'
    bl_options = {'REGISTER', 'UNDO'}

    atom_scale: option_property('atom_scale', name='Atom size multiplier', min=0.01, max=100)
    bond_thickness: option_property('bond_thickness', name='Bond thickness', api_default=True, min=0.001)
    colors: option_property('colors', name='Palette')
    surface_opacity: option_property('surface_opacity', name='Surface opacity', min=0, max=1)

    @classmethod
    def poll(cls, context):
        return context.mode == 'OBJECT' and resolve_owner(context.object) is not None

    def invoke(self, context, event):
        settings = read_record(resolve_owner(context.object))['settings']
        for name in ('atom_scale', 'bond_thickness', 'colors', 'surface_opacity'):
            setattr(self, name, settings[name])
        self._initial_opacity = self.surface_opacity
        return context.window_manager.invoke_props_dialog(self)

    def execute(self, context):
        opacity = self.surface_opacity if self.properties.is_property_set('surface_opacity') else None
        if opacity == getattr(self, '_initial_opacity', None):
            opacity = None
        try:
            restyle_molecule(resolve_owner(context.object), atom_scale=self.atom_scale,
                             bond_thickness=self.bond_thickness, colors=self.colors,
                             surface_opacity=opacity)
        except (ValueError, KeyError) as exc:
            self.report({'ERROR'}, str(exc))
            return {'CANCELLED'}
        return {'FINISHED'}


def make_molecule_independent(owner):
    """Detach geometry and top-level materials shared outside this hierarchy."""
    objects = {o for o in [owner, *owned_objects(owner)] if o.type in {'MESH', 'CURVE', 'SURFACE', 'META', 'FONT'}}
    if any(o.library for o in objects):
        raise ValueError('Make library-linked objects local in Blender before editing their molecule')
    geometry = {o.data for o in objects}
    geometry.update(o["mb_esp_original_mesh"] for o in objects if "mb_esp_original_mesh" in o)
    users = bpy.data.user_map(subset=geometry)
    copies = {data: data.copy() for data in geometry if data.library or users[data] - objects}
    for obj in objects:
        if obj.data in copies:
            obj.data = copies[obj.data]
        if obj.get("mb_esp_original_mesh") in copies:
            obj["mb_esp_original_mesh"] = copies[obj["mb_esp_original_mesh"]]
    geometry = {copies.get(data, data) for data in geometry}

    # Generated bevel objects live outside collections and can be shared even
    # when the centerline curves have already been made single-user.
    profiles = {}
    for curve in geometry:
        profile = getattr(curve, 'bevel_object', None)
        if profile is not None and profile.get('molecular_blender_bond_profile'):
            profiles.setdefault(profile, []).append(curve)
    for profile, curves in profiles.items():
        if (profile.library or profile.users_collection
                or any(c not in geometry and c.bevel_object == profile for c in bpy.data.curves)):
            local = profile.copy()
            local.data = profile.data.copy()
            local.use_fake_user = True
            for curve in curves:
                curve.bevel_object = local
        elif profile.data.users > 1:
            profile.data = profile.data.copy()

    materials = {m for data in geometry for m in data.materials if m}
    materials.update(slot.material for obj in objects for slot in obj.material_slots if slot.material)
    materials.update(o["mb_esp_original_material"] for o in objects if "mb_esp_original_material" in o)
    users = bpy.data.user_map(subset=materials)
    local_users = objects | geometry
    copies = {m: m.copy() for m in materials if m.library or users[m] - local_users}
    for copy in copies.values():
        copy["mb_style_material"] = True
    for data in geometry:
        for i, material in enumerate(data.materials):
            if material in copies:
                data.materials[i] = copies[material]
    for obj in objects:
        if obj.get("mb_esp_original_material") in copies:
            obj["mb_esp_original_material"] = copies[obj["mb_esp_original_material"]]
        for slot in obj.material_slots:
            if slot.link == 'OBJECT' and slot.material in copies:
                slot.material = copies[slot.material]


class MB_OT_make_independent(bpy.types.Operator):
    bl_idname = 'molecular_blender.make_independent'
    bl_label = 'Make Molecule Independent'
    bl_description = 'Give this molecule private geometry and materials wherever they are shared with other objects'
    bl_options = {'UNDO'}

    @classmethod
    def poll(cls, context):
        return context.mode == 'OBJECT' and resolve_owner(context.object) is not None

    def execute(self, context):
        try:
            make_molecule_independent(resolve_owner(context.object))
        except (ValueError, RuntimeError) as exc:
            self.report({'ERROR'}, str(exc))
            return {'CANCELLED'}
        return {'FINISHED'}


class MB_OT_style_selected(bpy.types.Operator):
    bl_idname = 'molecular_blender.style_selected'
    bl_label = 'Style Selected'
    bl_description = 'Color selected atoms, selected bonds and bonds between selected atoms'
    bl_options = {'REGISTER', 'UNDO'}

    color: FloatVectorProperty(name='Accent color', subtype='COLOR', size=3,
                               min=0., max=1., default=(0.95, 0.35, 0.05))
    restore: BoolProperty(options={'HIDDEN', 'SKIP_SAVE'})

    @classmethod
    def poll(cls, context):
        return context.mode == 'OBJECT' and resolve_owner(context.object) is not None

    def invoke(self, context, event):
        if self.restore:
            return self.execute(context)
        return context.window_manager.invoke_props_dialog(self)

    def execute(self, context):
        owner = resolve_owner(context.object)
        selected = [obj for obj in context.selected_objects if resolve_owner(obj) == owner]
        atoms = [obj for obj in selected if obj.get(ROLE) == 'atom']
        atom_ids = {obj[ID] for obj in atoms}
        bond_ids = {obj[ID] for obj in selected if obj.get(ROLE) == 'bond'}
        targets = atoms + [obj for obj in owned_objects(owner, 'bond')
                           if obj[ID] in bond_ids or set(obj['mb_endpoints']) <= atom_ids]
        if not targets:
            self.report({'INFO'}, 'Select atoms or bonds to style')
            return {'CANCELLED'}
        styler = get_styler(read_record(owner)['settings'])
        copies = {}
        for obj in targets:
            if self.restore and 'mb_accent_color' not in obj:
                continue
            symbol = obj.get('mb_element')
            color = (*self.color, 1.) if not self.restore else (
                styler.element_color(elements[symbol]) if symbol else styler.bond_color)
            for slot in obj.material_slots:
                material = slot.material
                if material is None:
                    continue
                key = (material, color)
                if key not in copies:
                    copy = material.copy()
                    copy['mb_style_material'] = True
                    _set_color(copy, color)
                    copies[key] = copy
                slot.link = 'OBJECT'
                slot.material = copies[key]
            if self.restore:
                del obj['mb_accent_color']
            else:
                obj['mb_accent_color'] = tuple(self.color)
        for old in {mat for mat, color in copies}:
            if old.users == 0 and old.get('mb_style_material'):
                bpy.data.materials.remove(old)
        return {'FINISHED'}


def register():
    bpy.utils.register_class(MB_OT_restyle)
    bpy.utils.register_class(MB_OT_reset_atom_bond_style)
    bpy.utils.register_class(MB_OT_style_selected)
    bpy.utils.register_class(MB_OT_make_independent)


def unregister():
    bpy.utils.unregister_class(MB_OT_make_independent)
    bpy.utils.unregister_class(MB_OT_style_selected)
    bpy.utils.unregister_class(MB_OT_reset_atom_bond_style)
    bpy.utils.unregister_class(MB_OT_restyle)
