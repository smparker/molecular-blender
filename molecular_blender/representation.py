"""Switch single-frame representations using the displayed atoms and saved bond graph."""
import bpy
from bpy.props import EnumProperty

from .bond_edit import _atom_objects, _remove_curve, _remove_unused_curve, _validate_geometry
from .geometry_import import rollback_new_ids, track_new
from .periodictable import elements
from .plotter import set_bond_order
from .scene_edit import _atom_radius, space_filling_bond, update_junctions
from .scene_state import ID, normalize_duplicate_ids, owned_objects, read_record, resolve_owner, write_record
from .stylers import get_styler

STYLES = (('sticks', 'Sticks', 'Element-colored bonds with small atom junctions'),
          ('bs', 'Ball and stick', 'Element spheres with single-color bonds'),
          ('fixedbs', 'Fixed ball and stick', 'Ball and stick with uniform bond thickness'),
          ('vdw', 'Space-filling (van der Waals)', 'Show element van der Waals spheres and hide bonds'))


def can_switch(settings):
    return (settings['plot_type'] != 'animate' and settings['object_type'] != 'wireframe'
            and settings['plot_style'] in {item[0] for item in STYLES})


def switch_representation(owner, style):
    normalize_duplicate_ids()
    record = read_record(owner)
    settings = record['settings']
    if not can_switch(settings) or style not in {item[0] for item in STYLES}:
        raise ValueError('Choose a single-frame molecule with individually selectable atoms')
    if style == settings['plot_style']:
        return
    if style == 'sticks' and settings['bond_thickness'] <= 0:
        raise ValueError('Choose a positive bond thickness before switching to sticks')
    atoms = _atom_objects(owner)
    radii = {atom: atom.get('mb_atom_radius', _atom_radius(atom['mb_element'], settings))
             for atom in atoms.values()}
    for atom, radius in radii.items():
        if radius <= 0 or atom.get('mb_style_scale', atom.get('mb_base_scale', [1])[0]) <= 0:
            raise ValueError('Cannot resize an atom with zero saved radius or scale')
    curves = owned_objects(owner, 'bond')
    _validate_geometry(curves)
    old_materials = {slot.material for obj in curves for slot in obj.material_slots
                     if slot.material and slot.material.get('mb_style_material')}
    groups = {}
    bpy.context.view_layer.update()
    for obj in curves:
        for hook in obj.modifiers:
            if hook.type == 'HOOK' and hook.object in radii and hook.object.matrix_world.determinant() == 0:
                raise ValueError('The hooked atom has a singular transform')
        groups.setdefault(obj[ID], []).append(obj)
    bonds = {bond['id']: bond for bond in record['bonds']}
    materials = {}
    styler = get_styler(settings)

    def material_for(atom, key):
        if atom and 'mb_accent_color' not in atom and atom.active_material:
            return atom.active_material
        if key not in materials:
            name = key + ('_mat' if atom else '_bond_mat')
            mat = bpy.data.materials.get(name) if settings['recycle_materials'] else None
            if mat is None:
                mat = styler.atom_material(name, elements[key]) if atom else styler.bond_material(name, key)
                mat['mb_style_material'] = True
            materials[key] = mat
        return materials[key]

    # Prepare all replacement curve data before modifying the original objects.
    staged = []
    surplus = []
    prior_style = settings['plot_style']
    bond_style = record.get('bond_representation', 'fixedbs') if prior_style == 'vdw' else prior_style
    rebuild = style != 'vdw' and style != bond_style
    with rollback_new_ids():
        for identifier, parts in (groups.items() if rebuild else ()):
            parts.sort(key=lambda obj: obj.data.bevel_factor_start)
            endpoints = [atoms.get(i) for i in bonds[identifier]['atoms']]
            if None in endpoints or len(parts) > 2:
                raise ValueError('The displayed bond no longer has its original atom endpoints or curve parts')
            symbols = [atom['mb_element'] for atom in endpoints]
            element = min(symbols, key=lambda key: elements[key].vdw)
            key = 'universal' if settings['universal_bonds'] else element
            radius = settings['bond_thickness']
            if style == 'bs' and not settings['universal_bonds']:
                radius = 0.1 * elements[element].vdw
            split = style == 'sticks' and symbols[0] != symbols[1]
            count = 2 if split else 1
            for i in range(count):
                old = parts[min(i, len(parts) - 1)]
                copy = track_new(old.copy())
                copy.data = track_new(old.data.copy())
                curve = copy.data
                symbol = symbols[i] if style == 'sticks' else None
                mat = old.active_material if 'mb_accent_color' in old else material_for(
                    endpoints[i] if symbol else None, symbol or key)
                if curve.materials:
                    curve.materials[0] = mat
                else:
                    curve.materials.append(mat)
                copy.material_slots[0].link = 'DATA'
                curve.bevel_factor_start, curve.bevel_factor_end = ((i * 0.5, (i + 1) * 0.5) if split else (0, 1))
                if abs(curve.bevel_depth - radius) > 1e-7:
                    curve.bevel_depth = radius
                    if old.get('mb_order', 1) > 1:
                        tilts = [p.tilt for p in curve.splines[0].bezier_points]
                        set_bond_order(copy, old['mb_order'])
                        for p, tilt in zip(curve.splines[0].bezier_points, tilts):
                            p.tilt = tilt
                staged.append((old, copy, i < len(parts), symbol))
            surplus.extend((obj, parts[0]) for obj in parts[count:])

    for old, copy, existing, symbol in staged:
        if existing:
            old_curve = old.data
            old.data = copy.data
            old.material_slots[0].link = 'DATA'
            bpy.data.objects.remove(copy, do_unlink=True)
            _remove_unused_curve(old_curve)
            obj = old
        else:
            obj = copy
            for collection in old.users_collection:
                collection.objects.link(obj)
            if old.select_get():
                obj.select_set(True)
        if symbol:
            obj['mb_element'] = symbol
        elif 'mb_element' in obj:
            del obj['mb_element']
        obj['mb_radius'] = obj.data.bevel_depth
        obj['mb_nominal_radius'] = settings['bond_thickness']
    for obj, primary in surplus:
        if obj.select_get():
            primary.select_set(True)
        if bpy.context.view_layer.objects.active == obj:
            bpy.context.view_layer.objects.active = primary
        _remove_curve(obj)
    for atom, radius in radii.items():
        atom['mb_atom_radius'] = radius
    if style == 'vdw':
        record['bond_representation'] = prior_style
    elif prior_style == 'vdw':
        record.pop('bond_representation', None)
    for obj in owned_objects(owner, 'bond'):
        space_filling_bond(obj, style == 'vdw')
    settings['plot_style'] = style
    update_junctions(owner, record)
    write_record(owner, record)
    for mat in old_materials:
        if not mat.users:
            bpy.data.materials.remove(mat)


class MB_OT_representation(bpy.types.Operator):
    bl_idname = 'molecular_blender.change_representation'
    bl_label = 'Change Representation'
    bl_description = 'Switch atom and bond geometry using the current molecule'
    bl_options = {'REGISTER', 'UNDO'}

    style: EnumProperty(name='Representation', items=STYLES)

    @classmethod
    def poll(cls, context):
        owner = resolve_owner(context.object)
        return context.mode == 'OBJECT' and owner is not None and can_switch(read_record(owner)['settings'])

    def invoke(self, context, event):
        self.style = read_record(resolve_owner(context.object))['settings']['plot_style']
        return context.window_manager.invoke_props_dialog(self)

    def execute(self, context):
        try:
            switch_representation(resolve_owner(context.object), self.style)
        except (ValueError, KeyError, RuntimeError) as exc:
            self.report({'ERROR'}, str(exc))
            return {'CANCELLED'}
        return {'FINISHED'}
