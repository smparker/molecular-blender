"""Undoable scene topology editing, using saved IDs and ordinary curve bonds."""
import math

import bpy
from bpy.props import EnumProperty, FloatProperty

from .diagnostics import report_exception
from .bond_perception import perceive_bond_orders
from .scene_resources import discard_profile
from .periodictable import elements
from .scene_geometry import bind_bond_endpoints, configure_bond_curve
from .scene_state import (read_record, write_record, resolve_owner, owned_objects,
                          new_id, tag, ID, ROLE, normalize_duplicate_ids)
from .scene_edit import _private_curve, apply_bond_profile, space_filling_bond, update_junctions
from .stylers import get_styler


def _atom_objects(owner):
    atoms = {}
    for obj in owned_objects(owner, 'atom'):
        if obj[ID] in atoms:
            raise ValueError('Duplicate atom identity within this molecule; duplicate the whole molecule instead')
        atoms[obj[ID]] = obj
    return atoms


def selected_bond(context):
    """Resolve either two atoms or one logical bond (including both colored halves)."""
    selected = list(context.selected_objects)
    owners = {resolve_owner(obj) for obj in selected}
    if len(owners) != 1 or None in owners:
        raise ValueError('Select atoms or a bond from one molecule')
    owner = owners.pop()
    if len(selected) == 2 and all(obj.get(ROLE) == 'atom' for obj in selected):
        return owner, [obj[ID] for obj in selected], None
    if selected and all(obj.get(ROLE) == 'bond' for obj in selected):
        identifiers = {obj[ID] for obj in selected}
        if len(identifiers) == 1:
            return owner, None, identifiers.pop()
    raise ValueError('Select exactly two atoms or one bond')


def _find_bond(record, atom_ids=None, bond_id=None):
    if bond_id is not None:
        bond = next((b for b in record['bonds'] if b['id'] == bond_id), None)
        if bond is None:
            raise ValueError('The selected bond is no longer in this molecule')
        return bond
    if atom_ids is None or len(atom_ids) != 2 or len(set(atom_ids)) != 2:
        raise ValueError('A bond needs two distinct atoms')
    valid = {a['id'] for a in record['atoms']}
    if not set(atom_ids) <= valid:
        raise ValueError('Bond endpoints must belong to this molecule')
    return next((b for b in record['bonds'] if set(b['atoms']) == set(atom_ids)), None)


def _remove_curve(obj):
    curve = obj.data
    bpy.data.objects.remove(obj, do_unlink=True)
    _remove_unused_curve(curve)


def _remove_unused_curve(curve):
    if curve.users:
        return
    profile = curve.bevel_object
    bpy.data.curves.remove(curve)
    discard_profile(profile)


def _new_curves(owner, record, bond, atoms):
    """Build split/unsplit curves directly; bind hooks in the current parent frame."""
    settings = record['settings']
    first, second = (atoms[i] for i in bond['atoms'])
    bpy.context.view_layer.update()
    inverse = owner.matrix_world.inverted()
    coords = [inverse @ obj.matrix_world.translation for obj in (first, second)]
    if (coords[1] - coords[0]).length < 1e-8:
        raise ValueError('Cannot bond atoms at the same position')
    symbols = [obj['mb_element'] for obj in (first, second)]
    envelope = settings['plot_style'] == 'vdw' and settings['plot_type'] != 'animate'
    style = record.get('bond_representation', 'fixedbs') if envelope else settings['plot_style']
    colored = style == 'sticks'
    split = colored and symbols[0] != symbols[1]
    radius = settings['bond_thickness']
    if style == 'bs' and not settings['universal_bonds']:
        radius = 0.1 * min(elements[s].vdw for s in symbols)
    if radius <= 0:
        raise ValueError('Set a positive bond thickness before adding a bond')
    collection = owner.users_collection[0]
    objects = []
    material = None
    try:
        if not colored:
            material = get_styler(settings).bond_material(owner.name + '_manual_bond', 'universal')
            material['mb_style_material'] = True
        for index in range(2 if split else 1):
            curve = bpy.data.curves.new(owner.name + '_bond', 'CURVE')
            obj = bpy.data.objects.new(owner.name + '_bond', curve)
            objects.append(obj)
            curve['mb_generated'] = True
            collection.objects.link(obj)
            obj.parent = owner
            configure_bond_curve(curve, coords, radius, index if split else None)
            if colored:
                atom = (first, second)[index]
                if 'mb_accent_color' in atom:
                    base = get_styler(settings).atom_material(owner.name + '_bond_' + symbols[index],
                                                             elements[symbols[index]])
                    base['mb_style_material'] = True
                    curve.materials.append(base)
                elif atom.active_material is not None:
                    curve.materials.append(atom.active_material)
                obj['mb_element'] = symbols[index]
            else:
                curve.materials.append(material)
            apply_bond_profile(obj, bond, bond['order'])
            tag(obj, owner, 'bond', bond['id'])
            obj['mb_endpoints'] = bond['atoms']
            obj['mb_order'] = bond['order']
            obj['mb_radius'] = radius
            obj['mb_nominal_radius'] = settings['bond_thickness']
            bpy.context.view_layer.update()
            bind_bond_endpoints(obj, first, second)
            if envelope:
                space_filling_bond(obj, True)
        return objects
    except Exception:
        for obj in objects:
            _remove_curve(obj)
        if material is not None and material.users == 0:
            bpy.data.materials.remove(material)
        raise


def _validate_geometry(objects):
    for obj in objects:
        if obj.type != 'CURVE' or len(obj.data.splines) != 1 or len(obj.data.splines[0].bezier_points) != 2:
            raise ValueError('Bond editing requires straight two-point curves')
        points = obj.data.splines[0].bezier_points
        if (points[1].co - points[0].co).length < 1e-8:
            raise ValueError('Cannot edit a zero-length bond')


def _apply_order(objects, bond):
    for obj in objects:
        _private_curve(obj)
        apply_bond_profile(obj, bond, bond['order'])
        obj['mb_order'] = bond['order']


def set_manual_bond(owner, order, *, atom_ids=None, bond_id=None):
    if isinstance(order, bool) or not isinstance(order, int) or order not in (1, 2, 3):
        raise ValueError('Bond order must be 1, 2, or 3')
    normalize_duplicate_ids()
    record = read_record(owner)
    if record['settings']['object_type'] == 'wireframe':
        raise ValueError('Manual bonds require individually selectable atoms')
    atoms = _atom_objects(owner)
    bond = _find_bond(record, atom_ids, bond_id)
    endpoints = atom_ids if bond is None else bond['atoms']
    if not set(endpoints) <= atoms.keys():
        raise ValueError('Both endpoint atoms must be present in the scene')
    if bond is None:
        bond = dict(id=new_id(), atoms=list(endpoints), order=order,
                    provenance='manual', locked=True, direction=None)
        _new_curves(owner, record, bond, atoms)
        record['bonds'].append(bond)
    else:
        objects = [obj for obj in owned_objects(owner, 'bond') if obj[ID] == bond['id']]
        _validate_geometry(objects)
        bond.update(order=order, provenance='manual', locked=True)
        if objects:
            _apply_order(objects, bond)
        else:
            _new_curves(owner, record, bond, atoms)
    if 'frame_orders' in bond:
        bond['frame_orders'] = [order if prior else 0 for prior in bond['frame_orders']]
    record.pop('perception', None)
    update_junctions(owner, record)
    write_record(owner, record)
    return bond['id']


def delete_bond(owner, *, atom_ids=None, bond_id=None):
    record = read_record(owner)
    bond = _find_bond(record, atom_ids, bond_id)
    if bond is None:
        raise ValueError('There is no bond between the selected atoms')
    for obj in list(owned_objects(owner, 'bond')):
        if obj[ID] == bond['id']:
            _remove_curve(obj)
    record['bonds'].remove(bond)
    record.pop('perception', None)
    # Topology, including deleted edges, is authoritative for later perception.
    record.setdefault('deleted_bonds', []).append(dict(atoms=bond['atoms'], provenance='manual'))
    update_junctions(owner, record)
    write_record(owner, record)


def reperceive_bonds(owner):
    """Solve the saved graph with manual orders constrained, applying only success.

    Explicit source orders are also protected. Connectivity is never regenerated,
    so manual additions/deletions survive. Coordinates are read from displayed
    atoms in the molecule frame, not from their old import positions.
    """
    record = read_record(owner)
    if record['settings']['plot_type'] == 'animate':
        raise ValueError('Post-import perception currently requires a static molecule')
    objects = _atom_objects(owner)
    atom_records = record['atoms']
    if len(objects) != len(atom_records):
        raise ValueError('Perception requires all atoms, including hydrogens, present in the scene')
    indices = {a['id']: i for i, a in enumerate(atom_records)}
    edges = [tuple(indices[i] for i in b['atoms']) for b in record['bonds']]
    fixed = {i: b['order'] for i, b in enumerate(record['bonds'])
             if b['locked'] or b['provenance'] == 'explicit'}
    bpy.context.view_layer.update()
    inverse = owner.matrix_world.inverted()
    coords = [inverse @ objects[a['id']].matrix_world.translation for a in atom_records]
    symbols = [a['element'] for a in atom_records]
    lengths = [(coords[i] - coords[j]).length / (elements[symbols[i]].covalent + elements[symbols[j]].covalent)
               for i, j in edges]
    report = {}
    orders, charges = perceive_bond_orders(symbols, edges,
        charge=record['settings']['molecular_charge'], bond_lengths=lengths,
        max_iterations=record['settings']['bond_order_max_iterations'],
        atom_ids=[a['id'] for a in atom_records], fixed_orders=fixed, diagnostics=report)
    curves = owned_objects(owner, 'bond')
    _validate_geometry(curves)
    for index, (bond, order) in enumerate(zip(record['bonds'], orders)):
        if index not in fixed:
            bond.update(order=order, provenance='perceived')
        _apply_order([obj for obj in curves if obj[ID] == bond['id']], bond)
    record['perception'] = [report]
    for atom, charge in zip(atom_records, charges):
        atom['formal_charge'] = charge
    update_junctions(owner, record)
    write_record(owner, record)


class MB_OT_set_bond(bpy.types.Operator):
    bl_idname = 'molecular_blender.set_bond'
    bl_label = 'Set bond order'
    bl_description = 'Add or update a bond and lock its manual order against perception'
    bl_options = {'REGISTER', 'UNDO'}
    order: EnumProperty(name='Order', items=(('1', 'Single', ''), ('2', 'Double', ''), ('3', 'Triple', '')))

    @classmethod
    def poll(cls, context):
        return context.mode == 'OBJECT' and resolve_owner(context.object) is not None

    def execute(self, context):
        try:
            owner, atoms, bond = selected_bond(context)
            set_manual_bond(owner, int(self.order), atom_ids=atoms, bond_id=bond)
        except (ValueError, KeyError) as exc:
            report_exception(self, exc)
            return {'CANCELLED'}
        return {'FINISHED'}


def _selected_multiple_bond(context):
    owner, atoms, bond_id = selected_bond(context)
    record = read_record(owner)
    bond = _find_bond(record, atoms, bond_id)
    parts = [] if bond is None else [
        obj for obj in owned_objects(owner, 'bond') if obj[ID] == bond['id']]
    if bond and 'frame_orders' in bond:
        if not any(not obj.hide_viewport and obj.get('mb_order', 1) > 1 for obj in parts):
            raise ValueError('Select a bond with multiple strands at the current frame')
        parts = [obj for obj in parts if obj.get('mb_order', 1) > 1]
    if not parts or any(obj.get('mb_order', 1) < 2 for obj in parts):
        raise ValueError('Select a displayed double or triple bond, or its two atoms')
    _validate_geometry(parts)
    return owner, record, bond, parts


class MB_OT_set_bond_tilt(bpy.types.Operator):
    bl_idname = 'molecular_blender.set_bond_tilt'
    bl_label = 'Set strand tilt'
    bl_description = 'Rotate the strands of a double or triple bond around its centerline'
    bl_options = {'REGISTER', 'UNDO'}
    tilt: FloatProperty(name='Strand tilt', subtype='ANGLE', min=-math.pi, max=math.pi)

    @classmethod
    def poll(cls, context):
        if context.mode != 'OBJECT':
            return False
        try:
            _selected_multiple_bond(context)
        except (ValueError, KeyError) as exc:
            cls.poll_message_set(str(exc))
            return False
        return True

    def invoke(self, context, event):
        _, _, bond, parts = _selected_multiple_bond(context)
        curve = parts[0].data
        self.tilt = (bond.get('tilt_offset', 0.0) if curve.get('mb_trajectory_tilt') else
                     math.remainder(curve.splines[0].bezier_points[0].tilt, math.tau))
        return context.window_manager.invoke_props_dialog(self)

    def draw(self, context):
        _, _, _, parts = _selected_multiple_bond(context)
        label = 'Offset from animated plane' if parts[0].data.get('mb_trajectory_tilt') else 'Strand tilt'
        self.layout.prop(self, 'tilt', text=label)

    def execute(self, context):
        try:
            owner, record, bond, parts = _selected_multiple_bond(context)
            if not math.isfinite(self.tilt):
                raise ValueError('Strand tilt must be finite')
            for obj in parts:
                _private_curve(obj)
                curve = obj.data
                if curve.get('mb_trajectory_tilt') and curve.animation_data and curve.animation_data.action:
                    action = curve.animation_data.action
                    if action.users > 1:
                        action = curve.animation_data.action = action.copy()
                    delta = self.tilt - bond.get('tilt_offset', 0.0)
                    for point in curve.splines[0].bezier_points:
                        fcurve = action.fcurves.find(point.path_from_id('tilt'))
                        if fcurve:
                            for key in fcurve.keyframe_points:
                                key.co_ui.y += delta
                            fcurve.update()
                else:
                    for point in curve.splines[0].bezier_points:
                        point.tilt = self.tilt
            if parts[0].data.get('mb_trajectory_tilt'):
                bond['tilt_offset'] = self.tilt
                scene = context.scene
                scene.frame_set(scene.frame_current, subframe=scene.frame_subframe)
            else:
                bond['tilt'] = self.tilt
            write_record(owner, record)
        except (ValueError, KeyError) as exc:
            report_exception(self, exc)
            return {'CANCELLED'}
        return {'FINISHED'}


class MB_OT_delete_bond(bpy.types.Operator):
    bl_idname = 'molecular_blender.delete_bond'
    bl_label = 'Delete bond'
    bl_options = {'REGISTER', 'UNDO'}
    poll = classmethod(MB_OT_set_bond.poll.__func__)

    def execute(self, context):
        try:
            owner, atoms, bond = selected_bond(context)
            delete_bond(owner, atom_ids=atoms, bond_id=bond)
        except (ValueError, KeyError) as exc:
            report_exception(self, exc)
            return {'CANCELLED'}
        return {'FINISHED'}


class MB_OT_reperceive(bpy.types.Operator):
    bl_idname = 'molecular_blender.reperceive'
    bl_label = 'Perceive unlocked bonds'
    bl_description = 'Recalculate orders on the current graph, preserving manual and source bond orders'
    bl_options = {'REGISTER', 'UNDO'}
    poll = classmethod(MB_OT_set_bond.poll.__func__)

    def execute(self, context):
        try:
            reperceive_bonds(resolve_owner(context.object))
        except (ValueError, KeyError) as exc:
            report_exception(self, exc)
            return {'CANCELLED'}
        return {'FINISHED'}


CLASSES = (MB_OT_set_bond, MB_OT_set_bond_tilt, MB_OT_delete_bond, MB_OT_reperceive)


def register():
    for cls in CLASSES:
        bpy.utils.register_class(cls)


def unregister():
    for cls in reversed(CLASSES):
        bpy.utils.unregister_class(cls)
