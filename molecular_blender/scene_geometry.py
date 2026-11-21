# Copyright (C) 2014-2025 Shane Parker, Joshua Szekely
# SPDX-License-Identifier: GPL-3.0-or-later

"""Build Blender geometry for atoms, bonds, rings and gradient arrows."""
from collections import namedtuple
from math import cos, sin, tau

import bpy

from .bond_orientation import strand_tilt
from .geometry_import import track_new
from .scene_resources import discard_profile
from .options import get_frame_mapper
from .periodictable import elements
from .util import Timer, stopwatch, unique_name


Style = namedtuple(
    'Style', ['atom_size', 'bond_size', 'atom_scaling', 'bond_scaling', 'split_bond'])
STYLES = {
    'vdw': Style(atom_size="scaled", bond_size="none",
                 atom_scaling=1.0, bond_scaling=0.0, split_bond=False),
    'bs': Style(atom_size="scaled", bond_size="vdw",
                atom_scaling=0.25, bond_scaling=0.1, split_bond=False),
    'fixedbs': Style(atom_size="scaled", bond_size="fixed",
                     atom_scaling=0.25, bond_scaling=1.0, split_bond=False),
    'sticks': Style(atom_size="bond", bond_size="fixed",
                    atom_scaling=0.2, bond_scaling=1.0, split_bond=True)}


def configure_bond_curve(curve, coords, radius, part=None):
    """Set up the same two-point extrusion for imported and manually added bonds."""
    curve.dimensions = '3D'
    curve.resolution_u = 2
    curve.bevel_depth = radius
    curve.use_fill_caps = True
    if part is not None:
        curve.bevel_factor_start, curve.bevel_factor_end = part * 0.5, (part + 1) * 0.5
    spline = curve.splines.new('BEZIER')
    spline.bezier_points.add(1)
    for point, co in zip(spline.bezier_points, coords):
        point.co = co
        point.handle_left_type = point.handle_right_type = 'AUTO'


def bind_bond_endpoints(bond, first, second):
    """Bind both Bézier points and their handles in the evaluated rest pose."""
    for index, atom in enumerate((first, second)):
        hook = bond.modifiers.new('Bond endpoint', 'HOOK')
        hook.object = atom
        hook.vertex_indices_set(range(3 * index, 3 * index + 3))
        hook.center = bond.data.splines[0].bezier_points[index].co
        hook.matrix_inverse = atom.matrix_world.inverted() @ bond.matrix_world


def _parent_geometry(context, pairs):
    """Keep source-coordinate locations, animation channels and the original selection convention."""
    if not pairs:
        return
    for child, parent in pairs:
        child.parent = parent
        child.matrix_parent_inverse = parent.matrix_world.inverted_safe()
    context.view_layer.update()
    for obj in context.selected_objects:
        obj.select_set(False)
    child, parent = pairs[-1]
    child.select_set(True)
    context.view_layer.objects.active = parent


def _multiple_bond_profile(radius, order):
    """Strand radius, center spacing and total profile radius."""
    strand_radius = 0.7 * radius
    spacing = 1.7 * radius
    return strand_radius, spacing, spacing * (order - 1) / 2 + strand_radius


def set_bond_order(obj, order, direction=None):
    """Extrude parallel strands along an existing straight two-point bond curve.

    ``direction`` is the desired strand separation in the curve's local space
    (the molecule coordinate axes for bonds made by PlotMolecule). Its component
    along the bond is ignored. Omitting it uses Blender's untwisted profile axis.
    The centerline, hooks, materials and bevel factors remain on the bond curve.
    """
    if order not in (1, 2, 3):
        raise ValueError("Bond order must be 1, 2 or 3")
    curve = obj.data
    if len(curve.splines) != 1 or len(curve.splines[0].bezier_points) != 2:
        raise ValueError("Expected a straight two-point Bezier bond")
    points = curve.splines[0].bezier_points
    tangent = points[1].co - points[0].co
    if tangent.length < 1e-8:
        raise ValueError("Cannot orient a zero-length bond")
    tangent.normalize()
    tilt = strand_tilt(tangent, direction) if order != 1 else 0.0

    old_profile = curve.bevel_object
    if order == 1:
        curve.bevel_mode = 'ROUND'
        curve.bevel_object = None
    else:
        radius = curve.bevel_depth
        profile = track_new(bpy.data.curves.new(obj.name + "_profile", type='CURVE'))
        profile.dimensions = '2D'
        profile.resolution_u = 12
        # Disconnected circles in one bevel object make disconnected tubes, while
        # retaining a single centerline and its endpoint hooks. Each strand is
        # slightly smaller than a single bond, with a visible gap between them.
        strand_radius, spacing, _ = _multiple_bond_profile(radius, order)
        for strand in range(order):
            offset = (strand - (order - 1) / 2) * spacing
            circle = profile.splines.new('BEZIER')
            circle.bezier_points.add(3)
            circle.use_cyclic_u = True
            for p, (x, y) in zip(circle.bezier_points,
                                 ((1, 0), (0, 1), (-1, 0), (0, -1))):
                p.co = (offset + strand_radius * x, strand_radius * y, 0)
                p.handle_left_type = p.handle_right_type = 'AUTO'
        profile_obj = track_new(bpy.data.objects.new(obj.name + "_profile", profile))
        # Referenced by the bond, but deliberately not linked into the scene:
        # the bevel source itself should never appear as visible geometry.
        # Blender does not count bevel references as object users on save.
        profile_obj.use_fake_user = True
        profile_obj["molecular_blender_bond_profile"] = True
        curve.bevel_mode = 'OBJECT'
        curve.bevel_object = profile_obj
        curve.twist_mode = 'Z_UP'
        # These centerlines are straight. With extra longitudinal samples,
        # Blender's bevel-factor trimming can duplicate a ring at the midpoint.
        # The resulting collapsed faces corrupt smooth normals at color splits.
        curve.resolution_u = 1
    for p in points:
        p.tilt = tilt
    discard_profile(old_profile)


@stopwatch("PlotMolecule")
def PlotMolecule(context, molecule, options):
    """Plots the given molecule object to the specified context"""

    clock = Timer()

    object_type = options["object_type"]
    center_of_mass = molecule.center_of_mass()

    # update molecule name to make sure it's unique
    molecule.name = unique_name(molecule.name, bpy.data.objects.keys())

    top_collection = context.view_layer.active_layer_collection.collection
    molecule_collection = track_new(bpy.data.collections.new(molecule.name))
    atom_collection = track_new(bpy.data.collections.new("{}-atoms".format(molecule.name)))
    top_collection.children.link(molecule_collection)
    molecule_collection.children.link(atom_collection)

    # Creates empty to collect all the new objects
    bpy.ops.object.empty_add(type='PLAIN_AXES', location=center_of_mass)
    molecule_obj = track_new(context.object)
    top_collection.objects.unlink(molecule_obj)
    molecule_collection.objects.link(molecule_obj)

    # set name for new parent
    molecule_obj.name = molecule.name
    molecule.name = molecule_obj.name
    molecule.scene_owner = molecule_obj
    molecule.scene_collection = molecule_collection

    plot_style = STYLES[options["plot_style"]]

    bond_thickness = options["bond_thickness"]
    draw_bonds = bond_thickness != 0.0
    draw_charges = options["charges"] != "none"

    # Stick models use tiny spheres as junctions. Widen just those junctions
    # touched by multiple bonds so the outer strands still meet their atoms.
    junction_scales = {}
    if draw_bonds and options["draw_bond_orders"] and plot_style.atom_size == "bond":
        for bond in molecule.bonds:
            max_order = max(bond.frame_orders, default=bond.order)
            if max_order > 1 and not (bond.iatom.hidden or bond.jatom.hidden):
                _, _, extent = _multiple_bond_profile(bond_thickness, max_order)
                for atom in (bond.iatom, bond.jatom):
                    junction_scales[atom.index] = max(
                        junction_scales.get(atom.index, 1.0), extent / bond_thickness)

    # Preserve the first available suffix without rescanning earlier atoms.
    obj_keys = set(bpy.data.objects.keys())
    next_suffix = {}

    # list of objects that need to be processed after scene creation
    link_to = []
    to_hook = []
    to_parent = []

    clock.tick_print("preamble")

    charges_collection = None
    charge_geometry = {}
    element_objects = {}
    atom_objects = {}
    # Check to see if original atom already exists.
    # if yes, create translated linked duplicate, if no, create new object
    for atom in molecule.atoms:
        if atom.hidden:
            continue
        base_atom = molecule.name + "_" + atom.element.symbol
        suffix = next_suffix.get(base_atom, 0)
        while f"{base_atom}{suffix}" in obj_keys:
            suffix += 1
        atom.name = f"{base_atom}{suffix}"
        next_suffix[base_atom] = suffix + 1
        obj_keys.add(atom.name)
        if atom.element.symbol in element_objects:
            atom_copy = track_new(element_objects[atom.element.symbol].copy())
            atom_copy.location = atom.position
            atom_copy.name = atom.name

            # store object in link_to list
            link_to.append((atom_copy, atom_collection))
            to_parent.append((atom_copy, molecule_obj))

            # hold onto reference
            atom_obj = atom_copy
        else:  # Create the base atom from which all other of same element will be copied
            radius = plot_style.atom_scaling * \
                atom.element.vdw if (plot_style.atom_size !=
                                     "bond") else bond_thickness

            if object_type == "nurbs":
                bpy.ops.surface.primitive_nurbs_surface_sphere_add(
                    radius=radius, location=atom.position)
            elif object_type == "meta":
                bpy.ops.object.metaball_add(
                    type='BALL', radius=radius, location=atom.position)
            else:
                bpy.ops.mesh.primitive_uv_sphere_add(
                    location=atom.position, radius=radius)

            atom_obj = track_new(context.object)
            track_new(atom_obj.data)
            bpy.ops.object.shade_smooth()

            atom_obj.name = atom.name
            atom_obj.data.name = atom.name
            element_objects[atom.element.symbol] = atom_obj
            top_collection.objects.unlink(atom_obj)

            atom_obj.data.materials.append(
                bpy.data.materials[molecule.materials[atom.element.symbol]])

            link_to.append((atom_obj, atom_collection))

            # make parent active
            to_parent.append((atom_obj, molecule_obj))

        atom.name = atom_obj.name
        atom_objects[atom] = atom_obj
        if plot_style.atom_size == "bond":
            atom_obj.scale = (junction_scales.get(atom.index, 1.0),) * 3

        if draw_charges:  # set a sphere on top of atom that will show charge
            if charges_collection is None:
                charges_collection = track_new(bpy.data.collections.new(f"{molecule.name}-charges"))
                molecule_collection.children.link(charges_collection)
            for sign, material, charge in (('plus', molecule.pluscharge_mat, atom.charge),
                                           ('neg', molecule.negcharge_mat, -atom.charge)):
                key = (atom_obj.data, material)
                if key not in charge_geometry:
                    geometry = track_new(atom_obj.data.copy())
                    geometry.materials[0] = bpy.data.materials[material]
                    charge_geometry[key] = geometry
                overlay = track_new(bpy.data.objects.new(f"{atom.name}_{sign}", charge_geometry[key]))
                setattr(atom, sign + '_charge', overlay.name)
                overlay.location = atom.position
                overlay.scale = molecule.scale(max(0.0, charge))
                link_to.append((overlay, charges_collection))
                to_parent.append((overlay, atom_obj))

    clock.tick_print("atom creation")

    if draw_bonds:  # Add the bonds
        bond_collection = track_new(bpy.data.collections.new("{}-bonds".format(molecule.name)))
        molecule_collection.children.link(bond_collection)
        for bond in molecule.bonds:  # make curves for bonds
            iatom = bond.iatom
            jatom = bond.jatom

            if iatom.hidden or jatom.hidden:
                continue

            rad = bond_thickness
            if not options['universal_bonds'] and plot_style.bond_size == "vdw":
                rad = elements[bond.style].vdw * plot_style.bond_scaling

            to_split = plot_style.split_bond and iatom.element != jatom.element
            orders = sorted(set(bond.frame_orders or (bond.order,)) - {0}) if options['draw_bond_orders'] else [1]
            base_names = bond.make_names(molecule.name, to_split)
            bond.name = []
            for order in orders:
                suffix = f'_order{order}' if len(orders) > 1 else ''
                names = [unique_name(base + suffix, obj_keys) for base in base_names]
                obj_keys.update(names)
                for i, name in enumerate(names):
                    atom = (iatom, jatom)[i]
                    material = (molecule.materials[atom.element.symbol] if plot_style.split_bond
                                else molecule.bond_materials[bond.style])
                    curve = track_new(bpy.data.curves.new(name, type='CURVE'))
                    curve.materials.append(bpy.data.materials[material])
                    coords = (end.position - molecule_obj.location for end in (iatom, jatom))
                    configure_bond_curve(curve, coords, rad, i if to_split else None)
                    obj = track_new(bpy.data.objects.new(name, curve))
                    if order > 1:
                        set_bond_order(obj, order, bond.direction)
                    obj.parent = molecule_obj
                    link_to.append((obj, bond_collection))
                    if bond.frame_orders:
                        obj['mb_profile_order'] = order
                    bond.name.append(obj.name)
                    to_hook.append((obj, atom_objects[iatom], atom_objects[jatom]))

    clock.tick_print("bond creation")

    # finalize by linking in objects and updating scene
    for ob, coll in link_to:
        coll.objects.link(ob)

    clock.tick_print("link objects")

    # Hooks need evaluated endpoint matrices, and charge overlays are parented
    # to newly created atoms. Plain atoms only need the existing root transform.
    if (options["hook_atoms"] and to_hook) or draw_charges or not to_parent:
        context.view_layer.update()

    clock.tick_print("update scene")

    if options["hook_atoms"]:
        for bond, first, second in to_hook:
            bind_bond_endpoints(bond, first, second)

        clock.tick_print("set hooks")

    _parent_geometry(context, to_parent)
    clock.tick_print("set parentage")
    if options["gradient"]:
        plot_gradients(molecule, options)

    return


def _gradient_arrow_mesh(name):
    """Closed unit arrow along +Z. The tip is exactly one vector length from the atom."""
    sides = 16
    circle = [(cos(tau * i / sides), sin(tau * i / sides)) for i in range(sides)]
    vertices = [(r*x, r*y, z) for r, z in ((1, 0), (1, .73), (2.5, .73)) for x, y in circle]
    tip = len(vertices)
    vertices.append((0, 0, 1))
    faces = [tuple(reversed(range(sides)))]
    for i in range(sides):
        j = (i + 1) % sides
        for base in (0, sides):
            faces.append((base+i, base+j, base+sides+j, base+sides+i))
        faces.append((2*sides+i, 2*sides+j, tip))
    mesh = track_new(bpy.data.meshes.new(name))
    mesh.from_pydata(vertices, [], faces)
    mesh.update()
    return mesh


def plot_gradients(molecule, options):
    mesh = None
    owner = molecule.scene_owner
    frame_for = get_frame_mapper(options, molecule.frame_numbers)
    animate = options['plot_type'] == 'animate' and options['object_type'] != 'wireframe'
    interpolation = {'hold': 'CONSTANT', 'linear': 'LINEAR'}.get(options['trajectory_interpolation'])
    for atom in molecule.atoms:
        samples = atom.trajectory if animate else (atom,)
        if atom.hidden or not any(sample.gradient.length > 1e-30 for sample in samples):
            continue
        if mesh is None:
            mesh = _gradient_arrow_mesh(molecule.name + '_gradient')
        name = atom.name or f"{molecule.name}_{atom.element.symbol}{atom.index}"
        arrow = track_new(bpy.data.objects.new(name + '_gradient', mesh))
        molecule.scene_collection.objects.link(arrow)
        arrow.parent = owner
        arrow.matrix_parent_inverse = owner.matrix_world.inverted_safe()
        arrow.rotation_mode = 'QUATERNION'
        for i, sample in enumerate(samples):
            arrow.location = sample.position
            length = sample.gradient.length
            radius = .04 * min(length, 1.)
            arrow.scale = (radius, radius, length)
            if length > 1e-30:
                rotation = sample.gradient.to_track_quat('Z', 'Y')
                rotation.make_compatible(arrow.rotation_quaternion)
                arrow.rotation_quaternion = rotation
            if animate:
                for path in ('location', 'rotation_quaternion', 'scale'):
                    arrow.keyframe_insert(data_path=path, frame=frame_for(i))
        if animate and interpolation:
            for curve in arrow.animation_data.action.fcurves:
                for key in curve.keyframe_points:
                    key.interpolation = interpolation


def PlotWireFrame(context, molecule, options):
    """driver for wireframe plotter"""
    for item in context.selectable_objects:
        item.select_set(False)

    verts = [atom.position for atom in molecule.atoms]
    indices = {atom.index: i for i, atom in enumerate(molecule.atoms)}
    bonds = [(indices[b.iatom.index], indices[b.jatom.index]) for b in molecule.bonds]

    molec_mesh = track_new(bpy.data.meshes.new(molecule.name))
    molec_mesh.from_pydata(verts, bonds, [])
    molec_mesh.update()
    molec_obj = track_new(bpy.data.objects.new(molecule.name, molec_mesh))
    molec_obj.data = molec_mesh
    context.collection.objects.link(molec_obj)
    molec_obj.select_set(True)
    context.view_layer.objects.active = molec_obj
    bpy.ops.object.convert(target='CURVE')
    track_new(context.object)
    track_new(context.object.data)
    context.object.data.fill_mode = 'FULL'
    context.object.data.render_resolution_u = 12
    context.object.data.bevel_depth = 0.1
    context.object.data.bevel_resolution = 4
    bpy.ops.object.shade_smooth()
    molecule.name = context.object.name
    molecule.scene_owner = context.object
    collection = track_new(bpy.data.collections.new(molecule.name))
    context.collection.children.link(collection)
    collection.objects.link(context.object)
    context.collection.objects.unlink(context.object)
    molecule.scene_collection = collection
    if options["gradient"]:
        plot_gradients(molecule, options)


@stopwatch("draw rings")
def plot_rings(context, molecule, options):
    """Given a list of atoms and bonds, determine where aromatic rings are and plot them"""
    planarCycles = molecule.rings

    ringmat = bpy.data.materials[molecule.materials["ring"]]
    ring_collection = track_new(bpy.data.collections.new('{}-rings'.format(molecule.name)))
    molecule.scene_collection.children.link(ring_collection)

    to_parent = []
    for cycle in planarCycles:
        objname = unique_name(molecule.name + "_ring", bpy.data.objects.keys())
        meshname = objname + "_mesh"
        ringMesh = track_new(bpy.data.meshes.new(meshname))
        verts = [molecule.atoms[iatom].position for iatom in cycle]
        bonds = [[imesh, (imesh + 1) % len(cycle)]
                 for imesh in range(len(cycle))]

        ringMesh.from_pydata(verts, bonds, [range(len(cycle))])
        ringMesh.update()
        ringObj = track_new(bpy.data.objects.new(objname, ringMesh))
        ringObj.data.materials.append(ringmat)
        ringObj.data = ringMesh
        ring_collection.objects.link(ringObj)
        to_parent.append((ringObj, molecule.scene_owner))
    _parent_geometry(context, to_parent)
