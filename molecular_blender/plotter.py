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
import bmesh

from .molecule import Molecule
from .periodictable import elements
from .stylers import get_styler
from .isosurfaces import cube_isosurface, molden_isosurface

from .util import stopwatch, Timer, unique_name
from .importers import molecule_from_file

def plot_prep(molecule, options):
    """Preprocess options"""
    process_options(options)
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
        pc = unique_name("pluscharge", bpy.data.materials.keys()) + "_mat"
        nc = unique_name("negcharge", bpy.data.materials.keys()) + "_mat"

        molecule.pluscharge_mat = pc
        molecule.negcharge_mat = nc

        styler.charge_material(pc, nc, elements[a])
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
        molecule.materials["ring"] = ringmat
        styler.ring_material(ringmat)
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
    molecule_collection = bpy.data.collections.new(molecule.name)
    atom_collection = bpy.data.collections.new("{}-atoms".format(molecule.name))
    top_collection.children.link(molecule_collection)
    molecule_collection.children.link(atom_collection)

    # Creates empty to collect all the new objects
    bpy.ops.object.empty_add(type='PLAIN_AXES', location=center_of_mass)
    molecule_obj = context.object
    top_collection.objects.unlink(molecule_obj)
    molecule_collection.objects.link(molecule_obj)

    # set name for new parent
    molecule_obj.name = molecule.name

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

    # local list of keys to keep making unique names
    obj_keys = bpy.data.objects.keys().copy()

    # list of objects that need to be processed after scene creation
    link_to = []
    to_hook = []
    to_parent = []

    clock.tick_print("preamble")

    # Check to see if original atom already exists.
    # if yes, create translated linked duplicate, if no, create new object
    for atom in molecule.atoms:
        if atom.hidden:
            continue
        # Unselect Everything
        bpy.ops.object.select_all(action='DESELECT')
        base_atom = molecule.name + "_" + atom.element.symbol
        ref_atom = base_atom + "0"
        atom_obj = ""
        if ref_atom in obj_keys:  # duplicate existing atom
            # name new object
            atom.name = unique_name(base_atom, obj_keys, 0)

            # Create the linked duplicate object
            atom_obj = bpy.data.objects[ref_atom]

            atom_copy = atom_obj.copy()
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
            atom.name = ref_atom
            bpy.ops.object.select_all(action='DESELECT')

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

            bpy.ops.object.shade_smooth()

            atom_obj = context.object

            atom_obj.name = ref_atom
            atom_obj.data.name = ref_atom
            top_collection.objects.unlink(atom_obj)

            obj_keys.append(ref_atom)

            atom_obj.data.materials.append(
                bpy.data.materials[molecule.materials[atom.element.symbol]])

            link_to.append((atom_obj, atom_collection))

            # make parent active
            to_parent.append((atom_obj, molecule_obj))

        if draw_charges:  # set a sphere on top of atom that will show charge
            charges_collection = bpy.data.collections.new("{}-charges".format(molecule.name))
            molecule_collection.children.link(charges_collection)

            plus_obj = atom_obj.copy()
            plus_obj.data = plus_obj.data.copy()
            plus_obj.data.materials[0] = bpy.data.materials[molecule.pluscharge_mat]
            plus_obj.name = atom.name + "_plus"
            atom.plus_charge = plus_obj.name
            link_to.append((plus_obj, charges_collection))

            neg_obj = atom_obj.copy()
            neg_obj.data = neg_obj.data.copy()
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
        bond_collection = bpy.data.collections.new("{}-bonds".format(molecule.name))
        molecule_collection.children.link(bond_collection)
        for bond in molecule.bonds:  # make curves for bonds
            # deselect all
            bpy.ops.object.select_all(action='DESELECT')

            iatom = bond.iatom
            jatom = bond.jatom

            if iatom.hidden or jatom.hidden:
                continue

            rad = bond_thickness
            if not options['universal_bonds'] and plot_style.bond_size == "vdw":
                rad = elements[i].vdw * plot_style.bond_scaling

            to_split = plot_style.split_bond and iatom.element != jatom.element
            bond.name = [ unique_name(b, obj_keys) for b in bond.make_names(
                molecule.name, to_split) ]
            obj_keys.extend(bond.name)
            if to_split:
                coords = [iatom.position, jatom.position]

                curves = [bpy.data.curves.new(
                    bname, type='CURVE') for bname in bond.name]
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

                    c_obj = bpy.data.objects.new(cname, c)
                    c_obj.data.bevel_depth = rad
                    c_obj.data.use_fill_caps = True
                    c_obj.data.bevel_factor_start, c_obj.data.bevel_factor_end = bounds
                    c_obj.parent = bpy.data.objects[molecule.name]

                    link_to.append((c_obj, bond_collection))
                    to_hook.append((iatom.name, jatom.name, cname))
            else:
                assert len(bond.name) == 1
                # Only one object should be associated with the bonds so only one name
                bondname = bond.name[0]
                # create curve object
                coords = [iatom.position, jatom.position]
                curve = bpy.data.curves.new(bondname, type='CURVE')
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
                curveOB = bpy.data.objects.new(bondname, curve)
                curveOB.data.bevel_depth = rad
                curveOB.data.use_fill_caps = True
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
            # Hook to atom1
            bpy.ops.object.select_all(action='DESELECT')
            bpy.data.objects[iatom].select_set(True)
            bpy.data.objects[bond].select_set(True)
            context.view_layer.objects.active = bpy.data.objects[bond]
            bpy.ops.object.mode_set(mode='EDIT')
            bpy.ops.curve.de_select_first()
            bpy.ops.object.hook_add_selob()
            bpy.ops.object.mode_set(mode='OBJECT')

            # Hook to atom2
            bpy.ops.object.select_all(action='DESELECT')
            bpy.data.objects[jatom].select_set(True)
            bpy.data.objects[bond].select_set(True)
            context.view_layer.objects.active = bpy.data.objects[bond]
            bpy.ops.object.mode_set(mode='EDIT')
            bpy.ops.curve.de_select_first()
            bpy.ops.curve.de_select_last()
            bpy.ops.object.hook_add_selob()
            bpy.ops.object.mode_set(mode='OBJECT')

        clock.tick_print("set hooks")

    for child, parent in to_parent:
        bpy.ops.object.select_all(action='DESELECT')
        child.select_set(True)
        context.view_layer.objects.active = parent
        bpy.ops.object.parent_set(type='OBJECT', keep_transform=False)

    clock.tick_print("set parentage")

    if options["gradient"]:
        for atom in molecule.atoms:
            if atom.hidden:
                continue
            if atom.gradient.length > 1.0e-30:
                bpy.ops.object.select_all(action='DESELECT')
                curve = bpy.data.curves.new(
                    atom.name + "_gradient", type='CURVE')
                curve.dimensions = '3D'
                curve.resolution_u = 2
                gradline = curve.splines.new('BEZIER')
                gradline.bezier_points.add(1)
                gradline.bezier_points[0].co = atom.position
                gradline.bezier_points[1].co = atom.position + atom.gradient
                for p in gradline.bezier_points:
                    p.handle_right_type = p.handle_left_type = 'AUTO'
                curveOB = bpy.data.objects.new(atom.name + "_gradient", curve)

                curveOB.data.use_fill_caps = True
                context.collection.objects.link(curveOB)
                context.view_layer.objects.active = curveOB
                curveOB.select_set(True)
                context.view_layer.objects.active = bpy.data.objects[molecule.name]
                bpy.ops.object.parent_set(type='OBJECT', keep_transform=False)

    return


def PlotWireFrame(context, molecule, _options):
    """driver for wireframe plotter"""
    for item in context.selectable_objects:
        item.select_set(False)

    verts = [atom.position for atom in molecule.atoms]
    bonds = [(b.iatom.index, b.jatom.index) for b in molecule.bonds]

    molec_mesh = bpy.data.meshes.new(molecule.name)
    molec_mesh.from_pydata(verts, bonds, [])
    molec_mesh.update()
    molec_obj = bpy.data.objects.new(molecule.name, molec_mesh)
    molec_obj.data = molec_mesh
    context.collection.objects.link(molec_obj)
    molec_obj.select_set(True)
    context.view_layer.objects.active = bpy.data.objects[molecule.name]
    bpy.ops.object.convert(target='CURVE')
    context.object.data.fill_mode = 'FULL'
    context.object.data.render_resolution_u = 12
    context.object.data.bevel_depth = 0.1
    context.object.data.bevel_resolution = 4
    bpy.ops.object.shade_smooth()

def AnimateMolecule(context, molecule, options):
    """Keyframes position of atoms and opacity of bonds"""
    frame_map = get_frame_mapper(options)

    # check to make sure trajectory information is stored for at least the first atom
    if not molecule.atoms[0].trajectory:
        raise Exception(
            "Trajectory must be read before calling AnimateMolecule")
    if all([atom.name == '' for atom in molecule.atoms]):
        raise Exception(
            "No names found for any atom in Molecule. Was PlotMolecule or PlotAtoms called?")
    for atom in filter(lambda x: x.name != '', molecule.atoms):
        if atom.name not in bpy.data.objects.keys():
            raise Exception(
                "Atom " + atom.name + " not found. Was PlotMolecule or PlotAtoms called properly?")

        if atom.hidden:
            continue

        # deselect everything for good measure
        for item in context.selectable_objects:
            item.select_set(False)

        atom_obj = bpy.data.objects[atom.name]
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

    if options["animate_bonds"] == "dynamic":
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


def create_mesh(name, verts, faces, material, context, remesh=8):
    """Some black magic to make a mesh with the given name, verts, and faces"""
    me = bpy.data.meshes.new(name)  # create a new mesh
    me.from_pydata(verts, [], faces)
    me.update()      # update the mesh with the new data

    bm = bmesh.new()
    bm.from_mesh(me)
    bmesh.ops.remove_doubles(bm, verts=bm.verts, dist=0.01)
    bm.to_mesh(me)

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


@stopwatch("draw surfaces")
def draw_surfaces(molecule, styler, context, options):
    """Draws isosurfaces"""

    # quick return if no data
    if molecule.volume is None and molecule.orbitals is None:
        return

    # master list of all sets of vertices for plotting
    vertex_sets = []
    # map set names to the frame at which they are active
    vertex_trajectory_map = {}

    fmap = get_frame_mapper(options)

    wm = context.window_manager

    # set up collections
    molecule_collection = bpy.data.collections[molecule.name]

    if molecule.volume is not None:  # volumetric data was read
        vol = molecule.volume
        isovals = options["isovalues"]
        if options["cumulative"]:
            isovals = vol.isovalue_containing_proportion(isovals)

        # marching cubes to make the surface
        vsets = cube_isosurface(vol.data, vol.origin, vol.axes, isovals, name="cube", wm=wm)
        if vsets:
            vertex_sets.extend(vsets)
        if molecule.volume_trajectory is not None:
            f0 = fmap(0)
            if f0 not in vertex_trajectory_map:
                vertex_trajectory_map[f0] = []
            vertex_trajectory_map[f0].extend( [ v["name"] for v in vsets ] )

            for i, vol in enumerate(molecule.volume_trajectory, 1):
                name = "cube-{:5d}".format(i)
                vsets = cube_isosurface(vol.data, vol.origin, vol.axes, isovals, name=name, wm=wm)
                fi = fmap(i)
                if fi not in vertex_trajectory_map:
                    vertex_trajectory_map[fi] = []
                vertex_trajectory_map[fi].extend( [ v["name"] for v in vsets ] )
                vertex_sets.extend(vsets)
    elif molecule.orbitals is not None:  # orbital data was read
        orbitals = molecule.orbitals
        homo = orbitals.homo()
        lumo = homo+1
        orbnames = [ x.strip() for x in options["orbital"].split(',') ]
        orblist = []
        for o in orbnames:
            neworb = None

            if o == "homo":
                neworb = homo
            elif o == "lumo":
                neworb = lumo

            mh = re.search(r"homo\s*-\s*([0-9]+)", o)
            if mh:
                neworb = homo - int(mh.group(1))
            ml = re.search(r"lumo\s*\+\s*([0-9]+)", o)
            if ml:
                neworb = lumo + int(ml.group(1))

            if orbitals.is_density(o):
                neworb = o

            if neworb is None:
                try:
                    neworb = int(o)
                except ValueError:
                    raise Exception(f"Could not understand orbital {o}! Looking for {orbitals.densities.keys()}?")
                    pass

            orblist.append(neworb)

        isovals = options["isovalues"]
        resolution = options["resolution"]

        for orbname, orbnumber in zip(orbnames, orblist):
            if orbitals.is_density(orbname):
                orb = orbitals.get_density(orbname)
            else:
                orb = orbitals.get_orbital(orbnumber)
            if options["cumulative"]:
                orbital_isovals = orb.isovalue_containing_proportion(isovals)
            else:
                orbital_isovals = isovals
            vset = molden_isosurface(orb, orbital_isovals, resolution, orbname, wm)
            if vset:
                vertex_sets.extend(vset)

            if molecule.orbitals_trajectory is not None:
                f0 = fmap(0)
                if f0 not in vertex_trajectory_map:
                    vertex_trajectory_map[f0] = []
                vertex_trajectory_map[f0].extend( [ v["name"] for v in vsets ] )

            if molecule.orbitals_trajectory is not None:
                for i, orbital_snap in enumerate(molecule.volume_trajectory, 1):
                    name = f"{orbname:s}-{i:4d}"
                    if orbname == "density":
                        orb = orbital_snap.get_density()
                    else:
                        orb = orbital_snap.get_orbital(orbnumber)
                    vsets = molecule_isosurface(orb, orbital_isovals, resolution, name, wm)

                    fi = fmap(i)
                    if fi not in vertex_trajectory_map:
                        vertex_trajectory_map[fi] = []
                    vertex_trajectory_map[fi].extend( [ v["name"] for v in vsets ] )

                    vertex_sets.extend(vsets)

    meshes = []
    materials = make_iso_materials(molecule, styler, vertex_sets, options)
    mol_obj = bpy.data.objects[molecule.name]

    # collect meshes by name
    function_names = set([ v["name"] for v in vertex_sets ])
    for func in function_names:
        fsets = [ v for v in vertex_sets if v["name"] == func ]
        fmeshes = []
        for v in fsets:
            # add surfaces to the scene
            name = "{0}_{1}_{2}".format(molecule.name, v["name"], v["isovalue"])
            verts, faces = create_geometry(v["triangles"])
            vmesh = create_mesh(name, verts, faces, v["material"], context, remesh=options["remesh"])
            v["mesh"] = vmesh
            fmeshes.append(vmesh)

        # create collection for this function
        function_collection = bpy.data.collections.new(f"{molecule.name}-isosurface-{func}")
        molecule_collection.children.link(function_collection)
        # link to molecule collection
        for m in fmeshes:
            function_collection.objects.link(m)
            m.select_set(True)

        meshes.extend(fmeshes)

    mol_obj.select_set(True)
    # parent them
    context.view_layer.objects.active = mol_obj
    bpy.ops.object.parent_set(type='OBJECT', keep_transform=False)

    if vertex_trajectory_map: # keyframe visibility for surfaces
        for frame in vertex_trajectory_map:
            for v in vertex_sets:
                obj = v["mesh"]
                vn = v["name"]
                hide = vn not in vertex_trajectory_map[frame]
                obj.hide_viewport = hide
                obj.hide_render = hide
                obj.keyframe_insert(data_path="hide_viewport", frame=frame)
                obj.keyframe_insert(data_path="hide_render", frame=frame)

@stopwatch("draw rings")
def plot_rings(context, molecule, options):
    """Given a list of atoms and bonds, determine where aromatic rings are and plot them"""
    planarCycles = molecule.rings

    ringmat = bpy.data.materials[molecule.materials["ring"]]
    ring_collection = bpy.data.collections.new('{}-rings'.format(molecule.name))
    bpy.data.collections[molecule.name].children.link(ring_collection)

    to_parent = []
    for cycle in planarCycles:
        objname = unique_name(molecule.name + "_ring", bpy.data.objects.keys())
        meshname = objname + "_mesh"
        ringMesh = bpy.data.meshes.new(meshname)
        verts = [molecule.atoms[iatom].position for iatom in cycle]
        bonds = [[imesh, (imesh + 1) % len(cycle)]
                 for imesh in range(len(cycle))]

        ringMesh.from_pydata(verts, bonds, [range(len(cycle))])
        ringMesh.update()
        ringObj = bpy.data.objects.new(objname, ringMesh)
        ringObj.data.materials.append(ringmat)
        ringObj.data = ringMesh
        ring_collection.objects.link(ringObj)
        to_parent.append(ringObj)

    molecule_obj = bpy.data.objects[molecule.name]
    for child in to_parent:
        bpy.ops.object.select_all(action='DESELECT')
        child.select_set(True)
        context.view_layer.objects.active = molecule_obj
        bpy.ops.object.parent_set(type="OBJECT", keep_transform=False)


def default_options(options):
    """add defaults to options in case called outside of blender"""
    # force defaults in case called outside of the blender UI
    defaults = {"bonds": True,
                "bond_thickness": 0.2,
                "hook_atoms": "auto",
                "plot_style": "sticks",
                "plot_type": "auto",
                "object_type": "mesh",
                "keystride": 1,
                "animate_bonds": "staticfirst",
                "universal_bonds": True,
                "ignore_hydrogen": True,
                "gradient": False,
                "charges": "none",
                "charges_offset": 0.0,
                "charges_factor": 1.0,
                "find_aromatic": False,
                "recycle_materials": False,
                "isovalues": "0.25,0.50",
                "cumulative": True,
                "volume": "orbital",
                "orbital": "homo",
                "remesh" : 6,
                "resolution": 0.05,
                "colors": "default"
               }

    for d in defaults:
        if d not in options:
            options[d] = defaults[d]

    return options


def process_options(options):
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

        if options["volume"] == "orbital":
            # automatically include positive and negative
            iso_set = set()
            for iso in isovals:
                iso_set.add(iso)
                iso_set.add(-iso)
            isovals = list(iso_set)

        options["isovalues"] = isovals

    return options


def get_frame_mapper(options):
    """ returns a function that maps snapshot number to frame number """
    stride = options["keystride"]
    def f(i):
        return i*stride + 1
    return f


def BlendMolecule(context, filename, **options):
    """basic driver that calls the appropriate plot functions"""

    options = default_options(options)
    name = filename.rsplit('.', 1)[0].rsplit('/')[-1]
    molecule = Molecule.from_dict(name, molecule_from_file(filename, options))
    molecule.determine_bonding(options)
    styler = get_styler(options)

    make_materials(molecule, styler, options)
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

    if options["isosurfaces"]:
        draw_surfaces(molecule, styler, context, options)
