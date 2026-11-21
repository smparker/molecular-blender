# -*- coding: utf-8 -*-
#
#  Molecular Blender
#  Filename: plotter.py
#  Copyright (C) 2014-2025 Shane Parker, Joshua Szekely
#
#  This file is part of Molecular Blender.
#
#  Molecular Blender is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 3, or (at your option)
#  any later version.
#
#  Molecular Blender is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU Library General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with Molecular Blender; see COPYING.
#  If not, see <http://www.gnu.org/licenses/>.
#

"""Methods to handle Blender interactions for Molecular Blender"""

from collections import namedtuple
from itertools import chain
import re

import bpy
import mathutils
from mathutils.kdtree import KDTree
import bmesh
import numpy as np

from .geometry_import import geometry_import, track_new
from .bond_orientation import strand_tilt, trajectory_tilts
from .molecule import Molecule
from .periodictable import elements
from .stylers import get_styler
from .isosurfaces import cube_isosurface, molden_isosurface

from .util import stopwatch, Timer, unique_name
from .importers import molecule_from_file


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
    if (old_profile is not None
            and old_profile.users == int(old_profile.use_fake_user)
            and old_profile.get("molecular_blender_bond_profile")
            and not any(c.bevel_object == old_profile for c in bpy.data.curves)):
        old_data = old_profile.data
        bpy.data.objects.remove(old_profile)
        if old_data.users == 0:
            bpy.data.curves.remove(old_data)


def plot_prep(molecule, options):
    """Preprocess options"""
    process_options(options, cube=molecule.volume is not None)
    if molecule.volume is None and molecule.orbitals is None:
        options["isosurfaces"] = False
    if options["charges"] != "none":
        molecule.chgfac = options["charge_factor"]
        molecule.chgoff = options["charge_offset"]
    return options

@stopwatch("make_materials")
def make_materials(molecule, styler, options):
    """Make all materials to be used by plotter"""

    make_atom_materials(molecule, styler, options)

    if molecule.bonds:
        make_bond_materials(molecule, styler, options)

    if molecule.rings:
        make_ring_materials(molecule, styler, options)

@stopwatch("make_atom_materials")
def make_atom_materials(molecule, styler, options):
    """Given a molecule object, creates required materials and populates materials dictionary"""

    # set of new atoms to make materials for
    atom_base_set = set([i.element.symbol for i in molecule.atoms])

    for a in atom_base_set:
        atom = a + "_mat"
        if options["recycle_materials"]:
            if atom in bpy.data.materials.keys():
                molecule.materials[a] = atom
                continue

        # get unique name for new material
        atom = unique_name(atom, bpy.data.materials.keys())
        molecule.materials[a] = atom  # store unique materials
        styler.atom_material(atom, elements[a])

    if options["charges"] != "none":  # make materials used for showing charges
        positive, negative = styler.charge_material("pluscharge_mat", "negcharge_mat", elements[a])
        molecule.pluscharge_mat = positive.name
        molecule.negcharge_mat = negative.name
    return


@stopwatch("make_bond_materials")
def make_bond_materials(molecule, styler, options):
    """Given a molecule object, creates the corresponding materials for the bonds"""

    # set of unique bond types
    bond_base_set = set([i.style for i in molecule.bonds])

    for b in bond_base_set:
        bond = b + "_bond_mat"
        if options["recycle_materials"]:
            if bond in bpy.data.materials.keys():
                molecule.bond_materials[b] = bond
                continue

        # unique name for bond material
        bond = unique_name(bond, bpy.data.materials.keys())
        molecule.bond_materials[b] = bond
        styler.bond_material(bond, b)
    return


@stopwatch("make_ring_materials")
def make_ring_materials(molecule, styler, options):
    """Given a molecule object, creates the corresponding materials for rings"""
    ringmat = "ring_mat"
    if not (options["recycle_materials"] and ringmat in bpy.data.materials.keys()):
        ringmat = unique_name(ringmat, bpy.data.materials.keys())
        styler.ring_material(ringmat)
    molecule.materials["ring"] = ringmat
    return


@stopwatch("make_iso_materials")
def make_iso_materials(molecule, styler, vertex_sets, options):
    """
    Given a molecule and an isosurface list, creates the materials
    for the isosurfaces and adds keys under "material" to each vertex set
    """

    orbitals = set([ v["name"] for v in vertex_sets ])

    # returns dict of materials created
    out = {}

    for o in orbitals:
        all_surfaces = [ v for v in vertex_sets if v["name"] == o ]

        # list of positive and negative isovalue surfaces
        positive_surfaces = [ v for v in all_surfaces if v["isovalue"] >= 0.0 ]
        negative_surfaces = [ v for v in all_surfaces if v["isovalue"] < 0.0 ]

        # sort surface lists by magnitude
        positive_surfaces.sort(key=lambda v: abs(v["isovalue"]), reverse=True)
        negative_surfaces.sort(key=lambda v: abs(v["isovalue"]), reverse=True)

        for surf_list, mat_basename in [ (positive_surfaces, "orbital_mat_plus"), (negative_surfaces, "orbital_mat_negative") ]:
            if not surf_list:
                continue

            nsurf = len(surf_list)
            vset = surf_list.pop(0)
            iv = vertex_sets.index(vset)
            matname = "{orbital}_{type}_inner".format(orbital=o, type=mat_basename)
            if not (options["recycle_materials"] and matname in bpy.data.materials.keys()):
                matname = unique_name(matname, bpy.data.materials.keys())
                if nsurf > 1:
                    # if more than one surface, first one is "inner"
                    material = styler.isosurface_material(matname)
                else:
                    # but the outer is so cool, so use that if there's only one
                    material = styler.outer_isosurface_material(matname)
                material['mb_surface_material'] = True
            else:
                material = bpy.data.materials.get(matname)
            out[matname] = material
            vertex_sets[iv]["material"] = material

            # any remaining elements get the "outer" material
            for surf in surf_list:
                matname = "{orbital}_{type}_outer".format(orbital=o, type=mat_basename)
                iv = vertex_sets.index(surf)
                if not (options["recycle_materials"] and matname in bpy.data.materials.keys()):
                    matname = unique_name(matname, bpy.data.materials.keys())
                    material = styler.outer_isosurface_material(matname)
                    material['mb_surface_material'] = True
                else:
                    material = bpy.data.materials.get(matname)
                out[matname] = material
                vertex_sets[iv]["material"] = material

    return out


@stopwatch("PlotMolecule")
def PlotMolecule(context, molecule, options):
    """Plots the given molecule object to the specified context"""
    global namedtuple

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
    molecule.scene_owner = molecule_obj
    molecule.scene_collection = molecule_collection

    Style = namedtuple(
        'Style', ['atom_size', 'bond_size', 'atom_scaling', 'bond_scaling', 'split_bond'])
    style_dict = {
        'vdw': Style(atom_size="scaled", bond_size="none",
                     atom_scaling=1.0, bond_scaling=0.0, split_bond=False),
        'bs': Style(atom_size="scaled", bond_size="vdw",
                    atom_scaling=0.25, bond_scaling=0.1, split_bond=False),
        'fixedbs': Style(atom_size="scaled", bond_size="fixed",
                         atom_scaling=0.25, bond_scaling=1.0, split_bond=False),
        'sticks': Style(atom_size="bond", bond_size="fixed",
                        atom_scaling=0.2, bond_scaling=1.0, split_bond=True)}
    plot_style = style_dict[options["plot_style"]]

    bond_thickness = options["bond_thickness"]
    draw_bonds = bond_thickness != 0.0
    draw_charges = options["charges"] != "none"

    # Stick models use tiny spheres as junctions. Widen just those junctions
    # touched by multiple bonds so the outer strands still meet their atoms.
    junction_scales = {}
    if draw_bonds and options["draw_bond_orders"] and plot_style.atom_size == "bond":
        for bond in molecule.bonds:
            if bond.order > 1 and not (bond.iatom.hidden or bond.jatom.hidden):
                _, _, extent = _multiple_bond_profile(bond_thickness, bond.order)
                for atom in (bond.iatom, bond.jatom):
                    junction_scales[atom.index] = max(
                        junction_scales.get(atom.index, 1.0), extent / bond_thickness)

    # local list of keys to keep making unique names
    obj_keys = bpy.data.objects.keys().copy()

    # list of objects that need to be processed after scene creation
    link_to = []
    to_hook = []
    to_parent = []

    clock.tick_print("preamble")

    element_objects = {}
    # Check to see if original atom already exists.
    # if yes, create translated linked duplicate, if no, create new object
    for atom in molecule.atoms:
        if atom.hidden:
            continue
        base_atom = molecule.name + "_" + atom.element.symbol
        atom.name = unique_name(base_atom, obj_keys, 0)
        if atom.element.symbol in element_objects:
            atom_copy = track_new(element_objects[atom.element.symbol].copy())
            atom_copy.location = atom.position
            atom_copy.name = atom.name

            # add new name to local list of keys
            obj_keys.append(atom.name)

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

            obj_keys.append(atom.name)

            atom_obj.data.materials.append(
                bpy.data.materials[molecule.materials[atom.element.symbol]])

            link_to.append((atom_obj, atom_collection))

            # make parent active
            to_parent.append((atom_obj, molecule_obj))

        if plot_style.atom_size == "bond":
            atom_obj.scale = (junction_scales.get(atom.index, 1.0),) * 3

        if draw_charges:  # set a sphere on top of atom that will show charge
            charges_collection = track_new(bpy.data.collections.new("{}-charges".format(molecule.name)))
            molecule_collection.children.link(charges_collection)

            plus_obj = track_new(atom_obj.copy())
            plus_obj.data = track_new(plus_obj.data.copy())
            plus_obj.data.materials[0] = bpy.data.materials[molecule.pluscharge_mat]
            plus_obj.name = atom.name + "_plus"
            atom.plus_charge = plus_obj.name
            link_to.append((plus_obj, charges_collection))

            neg_obj = track_new(atom_obj.copy())
            neg_obj.data = track_new(neg_obj.data.copy())
            neg_obj.data.materials[0] = bpy.data.materials[molecule.negcharge_mat]
            neg_obj.name = atom.name + "_neg"
            atom.neg_charge = neg_obj.name
            link_to.append((neg_obj, charges_collection))

            # initialize scale
            pc = max(0.0, atom.charge)
            plus_obj.scale = molecule.scale(pc)

            to_parent.append((plus_obj, atom_obj))

            # initialize scale
            nc = -min(0.0, atom.charge)
            neg_obj.scale = molecule.scale(nc)

            to_parent.append((neg_obj, atom_obj))

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
            bond.name = [ unique_name(b, obj_keys) for b in bond.make_names(
                molecule.name, to_split) ]
            obj_keys.extend(bond.name)
            if to_split:
                coords = [iatom.position, jatom.position]

                curves = [track_new(bpy.data.curves.new(bname, type='CURVE')) for bname in bond.name]
                for a, c in zip([iatom, jatom], curves):
                    c.materials.append(
                        bpy.data.materials[molecule.materials[a.element.symbol]])

                for c, cname, bounds in zip(curves, bond.name, [(0.0, 0.5), (0.5, 1.0)]):
                    c.dimensions = '3D'
                    c.resolution_u = 2

                    bondline = c.splines.new('BEZIER')
                    bondline.bezier_points.add(len(coords) - 1)

                    for i, pnt in enumerate(coords):
                        p = bondline.bezier_points[i]
                        p.co = pnt - molecule_obj.location
                        p.handle_right_type = p.handle_left_type = 'AUTO'

                    c_obj = track_new(bpy.data.objects.new(cname, c))
                    c_obj.data.bevel_depth = rad
                    c_obj.data.use_fill_caps = True
                    c_obj.data.bevel_factor_start, c_obj.data.bevel_factor_end = bounds
                    if options["draw_bond_orders"] and bond.order != 1:
                        set_bond_order(c_obj, bond.order, bond.direction)
                    c_obj.parent = bpy.data.objects[molecule.name]

                    link_to.append((c_obj, bond_collection))
                    to_hook.append((iatom.name, jatom.name, cname))
            else:
                assert len(bond.name) == 1
                # Only one object should be associated with the bonds so only one name
                bondname = bond.name[0]
                # create curve object
                coords = [iatom.position, jatom.position]
                curve = track_new(bpy.data.curves.new(bondname, type='CURVE'))
                bond_mat = bpy.data.materials[molecule.materials[iatom.element.symbol]] \
                    if plot_style.split_bond else \
                    bpy.data.materials[molecule.bond_materials[bond.style]]
                curve.materials.append(bond_mat)

                curve.dimensions = '3D'
                curve.resolution_u = 2
                bondline = curve.splines.new('BEZIER')
                bondline.bezier_points.add(len(coords) - 1)
                for i, pnt in enumerate(coords):
                    p = bondline.bezier_points[i]
                    p.co = pnt - molecule_obj.location
                    p.handle_right_type = p.handle_left_type = 'AUTO'
                curveOB = track_new(bpy.data.objects.new(bondname, curve))
                curveOB.data.bevel_depth = rad
                curveOB.data.use_fill_caps = True
                if options["draw_bond_orders"] and bond.order != 1:
                    set_bond_order(curveOB, bond.order, bond.direction)
                curveOB.parent = bpy.data.objects[molecule.name]

                link_to.append((curveOB, bond_collection))

                to_hook.append((iatom.name, jatom.name, bondname))

    clock.tick_print("bond creation")

    # finalize by linking in objects and updating scene
    for ob, coll in link_to:
        coll.objects.link(ob)

    clock.tick_print("link objects")

    # update scene
    context.view_layer.update()

    clock.tick_print("update scene")

    if options["hook_atoms"]:
        for iatom, jatom, bond in to_hook:
            bind_bond_endpoints(bpy.data.objects[bond], bpy.data.objects[iatom], bpy.data.objects[jatom])

        clock.tick_print("set hooks")

    if options["gradient"]:
        for atom in molecule.atoms:
            if atom.hidden:
                continue
            if atom.gradient.length > 1.0e-30:
                curve = track_new(bpy.data.curves.new(atom.name + "_gradient", type='CURVE'))
                curve.dimensions = '3D'
                curve.resolution_u = 2
                gradline = curve.splines.new('BEZIER')
                gradline.bezier_points.add(1)
                gradline.bezier_points[0].co = atom.position
                gradline.bezier_points[1].co = atom.position + atom.gradient
                for p in gradline.bezier_points:
                    p.handle_right_type = p.handle_left_type = 'AUTO'
                curveOB = track_new(bpy.data.objects.new(atom.name + "_gradient", curve))

                curveOB.data.use_fill_caps = True
                context.collection.objects.link(curveOB)
                to_parent.append((curveOB, molecule_obj))

    _parent_geometry(context, to_parent)
    clock.tick_print("set parentage")

    return


def PlotWireFrame(context, molecule, _options):
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

def AnimateMolecule(context, molecule, options):
    """Keyframes position of atoms and opacity of bonds"""
    frame_map = get_frame_mapper(options, molecule.frame_numbers)
    interpolation = {'hold': 'CONSTANT', 'linear': 'LINEAR'}.get(options.get('trajectory_interpolation'))

    # check to make sure trajectory information is stored for at least the first atom
    if not molecule.atoms[0].trajectory:
        raise Exception(
            "Trajectory must be read before calling AnimateMolecule")
    if all([atom.name == '' for atom in molecule.atoms]):
        raise Exception(
            "No names found for any atom in Molecule. Was PlotMolecule or PlotAtoms called?")
    for item in context.selected_objects:
        item.select_set(False)
    for atom in filter(lambda x: x.name != '', molecule.atoms):
        atom_obj = bpy.data.objects.get(atom.name)
        if atom_obj is None:
            raise Exception(
                "Atom " + atom.name + " not found. Was PlotMolecule or PlotAtoms called properly?")

        if atom.hidden:
            continue

        pluscharge = bpy.data.objects[atom.plus_charge] if (
            options["charges"] != "none") else ""
        negcharge = bpy.data.objects[atom.neg_charge] if (
            options["charges"] != "none") else ""

        for (iframe, isnap) in enumerate(atom.trajectory):
            atom_obj.location = isnap.position
            atom_obj.keyframe_insert(
                data_path='location', frame=frame_map(iframe))

            if options["charges"] != "none":
                if options["charges"] == "scale":
                    pc, nc = max(0.0, isnap.charge), -min(0.0, isnap.charge)
                    pluscharge.scale = molecule.scale(pc)
                    negcharge.scale = molecule.scale(nc)
                    for x in [pluscharge, negcharge]:
                        x.keyframe_insert(data_path='scale',
                                          frame=frame_map(iframe))

        if interpolation:
            animated = [(atom_obj, 'location')]
            if options['charges'] == 'scale':
                animated.extend(((pluscharge, 'scale'), (negcharge, 'scale')))
            for obj, path in animated:
                for curve in obj.animation_data.action.fcurves:
                    if curve.data_path == path:
                        for point in curve.keyframe_points:
                            point.interpolation = interpolation

    if options["draw_bond_orders"] and options["hook_atoms"]:
        _animate_strand_tilts(molecule, frame_map, interpolation)

    if options["animate_bonds"] == "dynamic" and molecule.explicit_bonds is None:
        bondmask = molecule.bond_mask(options)
        for bond in molecule.bonds:
            i, j = bond.iatom.index, bond.jatom.index
            for bond_obj in [ bpy.data.objects[b] for b in bond.name ]:
                for (iframe, mask) in enumerate(bondmask[(i, j)]):
                    bond_obj.hide_viewport = not mask
                    bond_obj.keyframe_insert(
                        data_path="hide_viewport", frame=frame_map(iframe))
                    bond_obj.hide_render = not mask
                    bond_obj.keyframe_insert(
                        data_path="hide_render", frame=frame_map(iframe))
    return


def _animate_strand_tilts(molecule, frame_map, interpolation):
    bonds = [bond for bond in molecule.bonds if bond.order > 1 and bond.name]
    if not bonds:
        return
    neighbors = {atom: [] for atom in molecule.atoms}
    for bond in molecule.bonds:
        neighbors[bond.iatom].append(bond.jatom)
        neighbors[bond.jatom].append(bond.iatom)
    for bond in bonds:
        curves = [bpy.data.objects[name].data for name in bond.name]
        for i, tilt in enumerate(trajectory_tilts(bond, neighbors)):
            for curve in curves:
                for point in curve.splines[0].bezier_points:
                    point.tilt = tilt
                    point.keyframe_insert(data_path='tilt', frame=frame_map(i))
        for curve in curves:
            curve['mb_trajectory_tilt'] = True
            if interpolation:
                for fcurve in curve.animation_data.action.fcurves:
                    for key in fcurve.keyframe_points:
                        key.interpolation = interpolation


def field_weld_distance(resolution):
    # Intersections are shared to float32 precision. A fixed 0.01 Å weld
    # also merges distinct vertices on fine grids and can tear a closed mesh.
    return resolution * 1e-4


def index_field_triangles(triangles, resolution):
    """Share float32 intersections and close sub-grid seams before allocating a Blender mesh."""
    points = np.ascontiguousarray(triangles, dtype=np.float32).reshape(-1, 3)
    if not len(points):
        return points, np.empty((0, 3), dtype=np.int32)
    # Sort complete XYZ triples as one key. Most shared edges match exactly.
    keys = points.view(np.dtype((np.void, points.strides[0]))).ravel()
    _, first, corners = np.unique(keys, return_index=True, return_inverse=True)
    vertices = points[first]
    faces = corners.reshape(-1, 3)

    # Adjacent marching cells can interpolate the same edge in opposite directions.
    # Merge residual roundoff by distance, including points across rounding-bin boundaries.
    tree = KDTree(len(vertices))
    for index, point in enumerate(vertices):
        tree.insert(point, index)
    tree.balance()
    representatives = np.arange(len(vertices), dtype=np.int32)
    tolerance = field_weld_distance(resolution)
    for index, point in enumerate(vertices):
        if representatives[index] == index:
            for _, other, _ in tree.find_range(point, tolerance):
                if other > index and representatives[other] == other:
                    representatives[other] = index
    faces = representatives[faces]
    a, b, c = faces.T
    faces = faces[(a != b) & (b != c) & (c != a)]
    # Collapsing a tiny sliver can also make two faces share all three corners.
    face_keys = np.sort(faces, axis=1).view(np.dtype((np.void, 3 * faces.itemsize))).ravel()
    _, first_face = np.unique(face_keys, return_index=True)
    faces = faces[np.sort(first_face)]
    used, indices = np.unique(faces, return_inverse=True)
    return vertices[used], indices.reshape(-1, 3).astype(np.int32)


def field_geometry(definition, resolution):
    if 'faces' in definition:
        return definition['vertices'], definition['faces']
    return index_field_triangles(definition['triangles'], resolution)


def _write_indexed_triangles(mesh, vertices, faces):
    """Fill Blender's vertex, loop and polygon arrays without a second mesh or weld."""
    mesh.vertices.add(len(vertices))
    mesh.vertices.foreach_set('co', vertices.ravel())
    mesh.loops.add(faces.size)
    mesh.loops.foreach_set('vertex_index', faces.ravel())
    mesh.polygons.add(len(faces))
    mesh.polygons.foreach_set('loop_start', np.arange(0, faces.size, 3, dtype=np.int32))
    if hasattr(mesh, 'shade_flat'):
        mesh.shade_flat()
    mesh.update(calc_edges=True)


def create_surface_mesh_data(name, verts, faces, weld_distance):
    me = bpy.data.meshes.new(name)
    try:
        if weld_distance is None:
            _write_indexed_triangles(me, verts, faces)
            return me
        me.from_pydata(verts, [], faces)
        me.update()
        bm = bmesh.new()
        try:
            bm.from_mesh(me)
            bmesh.ops.remove_doubles(bm, verts=bm.verts, dist=weld_distance)
            bm.to_mesh(me)
        finally:
            bm.free()
    except Exception:
        bpy.data.meshes.remove(me)
        raise
    return me


def create_mesh(name, verts, faces, material, context, remesh=8, weld_distance=0.01):
    me = create_surface_mesh_data(name, verts, faces, weld_distance)
    ob = bpy.data.objects.new(name, me)  # create a new object
    ob.data = me          # link the mesh data to the object
    if ob.data.materials:
        ob.data.materials[0] = material
    else:
        ob.data.materials.append(material)

    if remesh > 0:
        mod = ob.modifiers.new('Remesh', 'REMESH')
        mod.mode = 'SMOOTH'
        mod.octree_depth = remesh
        mod.scale = 0.99
        mod.use_smooth_shade = True
        mod.use_remove_disconnected = False

    return ob


def create_geometry(verts):
    """Build geometry out of vertices"""
    faces = []
    faceoffset = 0
    for ver in verts:
        if len(ver) == 4:
            faces.append((faceoffset + 0, faceoffset + 1,
                          faceoffset + 2, faceoffset + 3))
            faceoffset += 4
        elif len(ver) == 3:
            faces.append((faceoffset + 0, faceoffset + 1, faceoffset + 2))
            faceoffset += 3
    return list(chain.from_iterable(verts)), faces


def compute_surface_vertex_sets(molecule, options, progress_state=None, *, indexed=False):
    """Compute isosurface triangle data, synchronously or in an isolated worker process.

    progress_state is an optional dict with keys 'completed' and 'total' that this function
    updates as isosurfaces are computed, for use by a modal progress bar.
    With indexed=True, Molden fields carry shared vertices and triangle indices.
    """
    if molecule.volume is None and molecule.orbitals is None:
        return [], {}

    vertex_sets = []
    vertex_trajectory_map = {}
    fmap = get_frame_mapper(options, getattr(molecule, "frame_numbers", None))

    # Count total tasks for progress tracking
    if progress_state is not None:
        if molecule.volume is not None:
            total = 1 + (len(molecule.volume_trajectory) if molecule.volume_trajectory else 0)
        elif molecule.orbitals is not None:
            total = (1 if options["volume"] in ("density", "alpha", "beta", "spin") else
                     len([x.strip() for x in options["orbital"].split(',')]))
            total *= 1 + len(molecule.orbitals_trajectory or [])
        else:
            total = 1
        progress_state["total"] = max(total, 1)
        progress_state["completed"] = 0

    def tick():
        if progress_state is not None:
            progress_state["completed"] = progress_state.get("completed", 0) + 1

    if molecule.volume is not None:  # volumetric data was read
        frames = [molecule.volume, *(molecule.volume_trajectory or ())]
        for vol in frames:
            vol.select_dataset(options.get("cube_dataset", 0))
        levels = options["isovalues"]
        if options["cumulative"]:
            reference = frames[0]
            if reference.linear_combination:
                reference = next((vol for vol in frames if np.any(vol.data)), None)
                if reference is None:
                    for _ in frames:
                        tick()
                    return [], {}
            levels = reference.isovalue_containing_proportion(levels, options.get("cube_field_type", "auto"))
        for i, vol in enumerate(frames):
            name = "cube-{:5d}".format(i) if i else "cube"
            surfaces = cube_isosurface(vol.data, vol.origin, vol.axes, levels, name=name)
            if vol.linear_combination:
                surfaces = [s for s in surfaces if len(s['triangles'])]
            for surface in surfaces:
                surface['series'] = 'cube'
            if len(frames) > 1:
                vertex_trajectory_map[fmap(i)] = [s['name'] for s in surfaces]
            vertex_sets.extend(surfaces)
            tick()

    elif molecule.orbitals is not None:  # orbital data was read
        orbitals = molecule.orbitals
        channel = options.get('orbital_spin', 'all')
        channel = None if channel == 'all' else channel
        orbnames = ([options["volume"]] if options["volume"] in ('density', 'alpha', 'beta', 'spin') else
                    [x.strip().lower() for x in options['orbital'].split(',')])
        orblist = [name if orbitals.is_density(name) else orbitals.resolve_orbital(name, channel)
                   for name in orbnames]

        isovals = options["isovalues"]
        resolution = options["resolution"]

        def field_surface(model, name, number, label, levels=None):
            if model.is_density(name):
                if not np.any(model.density_matrix(name)[0]):
                    return [], levels
                field = model.get_density(name)
            else:
                field = model.get_orbital(model.resolve_orbital(name, channel))
            if levels is None:
                levels = field.isovalue_containing_proportion(isovals) if options['cumulative'] else isovals
            method = 'adaptive' if model.is_density(name) else 'tiled'
            surfaces = molden_isosurface(field, levels, resolution, label, method=method)
            if indexed:
                for surface in surfaces:
                    surface['vertices'], surface['faces'] = index_field_triangles(surface.pop('triangles'), resolution)
            return surfaces, levels

        for orbname, orbnumber in zip(orbnames, orblist):
            vset, orbital_isovals = field_surface(orbitals, orbname, orbnumber, orbname)
            tick()
            if vset:
                vertex_sets.extend(vset)

            if molecule.orbitals_trajectory is not None:
                f0 = fmap(0)
                vertex_trajectory_map.setdefault(f0, []).extend([v["name"] for v in vset])
                for i, orbital_snap in enumerate(molecule.orbitals_trajectory, 1):
                    snap_name = f"{orbname:s}-{i:4d}"
                    vsets_snap, _ = field_surface(orbital_snap, orbname, orbnumber, snap_name, orbital_isovals)
                    for v in vsets_snap:
                        v['series'] = orbname
                    tick()
                    fi = fmap(i)
                    vertex_trajectory_map.setdefault(fi, []).extend([v["name"] for v in vsets_snap])
                    vertex_sets.extend(vsets_snap)

    return vertex_sets, vertex_trajectory_map


def create_surface_meshes(molecule, styler, context, options, vertex_sets, vertex_trajectory_map,
                          *, select_surfaces=True):
    """Create Blender mesh objects from precomputed isosurface data. Must run on main thread."""
    if not vertex_sets:
        return

    molecule_collection = getattr(molecule, "scene_collection", None) or bpy.data.collections[molecule.name]
    make_iso_materials(molecule, styler, vertex_sets, options)
    mol_obj = getattr(molecule, "scene_owner", None) or bpy.data.objects[molecule.name]

    # Vertices are in source coordinates. The imported empty starts at the
    # center of mass, while a wireframe root starts at the source origin.
    parent_inverse = mathutils.Matrix.Identity(4)
    if options["object_type"] != "wireframe":
        parent_inverse = mathutils.Matrix.Translation(-molecule.center_of_mass())
    function_names = dict.fromkeys(v["name"] for v in vertex_sets)
    for func in function_names:
        fsets = [v for v in vertex_sets if v["name"] == func]
        fmeshes = []
        for v in fsets:
            name = "{0}_{1}_{2}".format(molecule.name, v["name"], v["isovalue"])
            if molecule.volume is not None:
                verts, faces = create_geometry(v["triangles"])
                weld_distance = .01
            else:
                verts, faces = field_geometry(v, options['resolution'])
                weld_distance = None
            vmesh = create_mesh(name, verts, faces, v["material"], context,
                                remesh=options["remesh"], weld_distance=weld_distance)
            v["mesh"] = vmesh
            fmeshes.append(vmesh)

        function_collection = bpy.data.collections.new(f"{molecule.name}-isosurface-{func}")
        function_collection['mb_surface_collection'] = True
        molecule_collection.children.link(function_collection)
        for m in fmeshes:
            function_collection.objects.link(m)
            m.parent = mol_obj
            m.matrix_parent_inverse = parent_inverse
            if select_surfaces:
                m.select_set(True)

    if select_surfaces:
        mol_obj.select_set(True)
        context.view_layer.objects.active = mol_obj

    if hasattr(molecule, 'scene_owner'):
        from .scene_state import persist_surfaces
        persist_surfaces(mol_obj, vertex_sets, options)
        if options.get('surface_opacity', 1.0) != 1.0:
            from .scene_edit import _set_opacity
            for definition in vertex_sets:
                obj = definition['mesh']
                material = obj.active_material.copy()
                material['mb_style_material'] = True
                obj.material_slots[0].link = 'OBJECT'
                obj.material_slots[0].material = material
                _set_opacity(material, options['surface_opacity'])

    if vertex_trajectory_map:
        for frame in vertex_trajectory_map:
            for v in vertex_sets:
                obj = v["mesh"]
                vn = v["name"]
                hide = vn not in vertex_trajectory_map[frame]
                obj.hide_viewport = hide
                obj.hide_render = hide
                obj.keyframe_insert(data_path="hide_viewport", frame=frame)
                obj.keyframe_insert(data_path="hide_render", frame=frame)

    if hasattr(molecule, 'scene_owner'):
        from .scene_surfaces import SOLO_GROUP, apply_solo, surface_groups
        if SOLO_GROUP in mol_obj:
            from .scene_state import read_record
            apply_solo(mol_obj, surface_groups(mol_obj, read_record(mol_obj)))


@stopwatch("draw surfaces")
def draw_surfaces(molecule, styler, context, options):
    """Draws isosurfaces"""
    vertex_sets, vertex_trajectory_map = compute_surface_vertex_sets(molecule, options, indexed=True)
    create_surface_meshes(molecule, styler, context, options, vertex_sets, vertex_trajectory_map)

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


def default_options(options):
    """add defaults to options in case called outside of blender"""
    from .options import DEFAULTS, molden_convention
    options['molden_convention'] = molden_convention(options)
    for name, value in DEFAULTS.items():
        options.setdefault(name, value)
    return options


def process_options(options, *, cube=False):
    """postprocess choices that might interact with each other"""

    hooking = {"on": True, "off": False,
               "auto": options["plot_type"] == "animate" }
    options["hook_atoms"] = hooking[options["hook_atoms"]]

    options["isosurfaces"] = "isovalues" in options and options["isovalues"] != ""
    if options["isosurfaces"]:
        isovals = options["isovalues"].split(',')
        isovals = [float(x) for x in isovals]
        if not isovals:
            raise Exception("Isovalues must be specified")

        both_signs = options["volume"] in ("orbital", "spin")
        if cube:
            both_signs = options.get("cube_both_signs", both_signs)
        if both_signs:
            # automatically include positive and negative
            iso_set = set()
            for iso in isovals:
                iso_set.add(iso)
                iso_set.add(-iso)
            isovals = list(iso_set)

        options["isovalues"] = isovals

    return options


def get_frame_mapper(options, frame_numbers=None):
    """Map a snapshot position to a Blender frame, including explicit manifest timing."""
    stride = options["keystride"]
    def f(i):
        return (frame_numbers[i] if frame_numbers is not None else i) * stride + 1
    return f


def BlendMoleculeSetup(context, filename, **options):
    """Load molecule and create atoms/bonds. Returns (molecule, styler, options).

    Does not draw isosurfaces - call draw_surfaces separately if needed.
    """
    from .options import validated_options
    options = validated_options(options)
    name = filename.rsplit('.', 1)[0].rsplit('/')[-1]
    input_data = molecule_from_file(filename, options)
    molecule = Molecule.from_dict(name, input_data)
    molecule.surface_input = input_data
    molecule.determine_bonding(options)
    styler = get_styler(options)
    molecule.settings = dict(options)

    with geometry_import(context):
        make_materials(molecule, styler, options)
        modes = input_data.get('vibrations')
        if options['object_type'] == 'wireframe' or input_data.get('frame_numbers'):
            modes = None
        if modes and options['hook_atoms'] == 'auto':
            options['hook_atoms'] = 'on'
        plot_prep(molecule, options)

        if options["object_type"] == "wireframe":
            PlotWireFrame(context, molecule, options)
            if options["plot_type"] == "animate":
                print("animation of wireframes not currently supported!")
        else:
            PlotMolecule(context, molecule, options)
            if options["plot_type"] == "animate":
                AnimateMolecule(context, molecule, options)

        if options["find_aromatic"]:
            plot_rings(context, molecule, options)

        from .scene_state import persist_molecule
        persist_molecule(molecule, filename, molecule.settings)
        if modes:
            from .vibrations import store_modes
            store_modes(molecule.scene_owner, modes)
        if molecule.settings.get('atom_scale', 1.0) != 1.0:
            from .scene_edit import restyle_molecule
            restyle_molecule(molecule.scene_owner, atom_scale=molecule.settings['atom_scale'],
                             bond_thickness=molecule.settings['bond_thickness'],
                             colors=molecule.settings['colors'],
                             surface_opacity=molecule.settings.get('surface_opacity', 1.0))
    return molecule, styler, options


def BlendMolecule(context, filename, **options):
    """basic driver that calls the appropriate plot functions"""
    molecule, styler, options = BlendMoleculeSetup(context, filename, **options)
    if options["isosurfaces"]:
        draw_surfaces(molecule, styler, context, options)
