# -*- coding: utf-8 -*-
#
#  Molecular Blender
#  Filename: plotter.py
#  Copyright (C) 2014 Shane Parker, Joshua Szekely
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
from .marching_cube import cube_isosurface, molden_isosurface

from .util import stopwatch, Timer, unique_name
from .importers import molecule_from_file


def plot_prep(molecule, options):
    """Preprocess options"""
    if options["charges"] != "none":
        molecule.chgfac = options["charge_factor"]
        molecule.chgoff = options["charge_offset"]

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
    """Given a molecule and an isosurface list, creates the materials for the isosurfaces"""

    orbitals = set([ v["name"] for v in vertex_sets ])

    out = {}
    for o in orbitals:
        for mat in ["orbital_mat_plus", "orbital_mat_minus" ]:
            key = "%s_%s" % (o, mat)
            if not (options["recycle_materials"] and mat in bpy.data.materials.keys()):
                imat = unique_name(mat, bpy.data.materials.keys())
                material = styler.isosurface_material(imat)
                out[key] = material
            else:
                out[key] = bpy.data.materials.get(mat)

    return out


@stopwatch("PlotMolecule")
def PlotMolecule(context, molecule, options):
    """Plots the given molecule object to the specified context"""
    global namedtuple

    clock = Timer()

    object_type = options["object_type"]
    center_of_mass = molecule.center_of_mass()

    # Creates empty to collect all the new objects
    bpy.ops.object.empty_add(type='PLAIN_AXES', location=center_of_mass)

    # update molecule name to make sure it's unique
    molecule.name = unique_name(molecule.name, bpy.data.objects.keys())
    # set name for new parent
    context.object.name = molecule.name

    # hold onto reference to parent
    molecule_obj = bpy.data.objects[molecule.name]

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

    # local list of keys to keep making unique names
    obj_keys = bpy.data.objects.keys().copy()

    # list of objects that need to be processed after scene creation
    to_link = []
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
        base_atom = molecule.name + "_" + atom.element.name
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

            # store object in to_link list
            to_link.append(atom_copy)

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
                    location=atom.position, size=radius)
            context.object.name = ref_atom
            context.object.data.name = ref_atom
            obj_keys.append(ref_atom)

            # hold onto reference
            atom_obj = context.object

            context.object.data.materials.append(
                bpy.data.materials[molecule.materials[atom.element.symbol]])
            bpy.ops.object.shade_smooth()

            # make parent active
            context.scene.objects.active = molecule_obj
            bpy.ops.object.parent_set(type='OBJECT', keep_transform=False)

        if options["charges"] != "none":  # set a sphere on top of atom that will show charge
            plus_obj = atom_obj.copy()
            plus_obj.data = plus_obj.data.copy()
            plus_obj.data.materials[0] = bpy.data.materials[molecule.pluscharge_mat]
            plus_obj.name = atom.name + "_plus"
            atom.plus_charge = plus_obj.name
            to_link.append(plus_obj)

            neg_obj = atom_obj.copy()
            neg_obj.data = neg_obj.data.copy()
            neg_obj.data.materials[0] = bpy.data.materials[molecule.negcharge_mat]
            neg_obj.name = atom.name + "_neg"
            atom.neg_charge = neg_obj.name
            to_link.append(neg_obj)

            # initialize scale
            pc = max(0.0, atom.charge)
            plus_obj.scale = molecule.scale(pc)

            to_parent.append((plus_obj, atom_obj))

            # initialize scale
            nc = -min(0.0, atom.charge)
            neg_obj.scale = molecule.scale(nc)

            to_parent.append((neg_obj, atom_obj))

    clock.tick_print("atom creation")

    if plot_style.bond_size != "none":  # Add the bonds
        # Make circles the correct size for the bonds
        bevnames = {}

        bond_bevels = set([i.style for i in molecule.bonds])
        # set up bevel types
        for i in bond_bevels:
            rad = bond_thickness
            if not options['universal_bonds'] and plot_style.bond_size == "vdw":
                rad = elements[i].vdw * plot_style.bond_scaling

            bpy.ops.curve.primitive_bezier_circle_add(
                radius=rad, location=(0, 0, 0))
            bevelname = unique_name(
                molecule.name + "_" + i + "_bond", bpy.data.objects.keys())
            bevnames[i] = bevelname
            context.object.name = bevelname

            context.scene.objects.active = molecule_obj
            bpy.ops.object.parent_set(type='OBJECT', keep_transform=False)

        for bond in molecule.bonds:  # make curves for bonds
            # deselect all
            bpy.ops.object.select_all(action='DESELECT')

            iatom = bond.iatom
            jatom = bond.jatom

            if iatom.hidden or jatom.hidden:
                continue

            bond.bevelname = bevnames[bond.style]

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
                    c_obj.data.bevel_object = bpy.data.objects[bond.bevelname]
                    c_obj.data.use_fill_caps = True
                    c_obj.data.bevel_factor_start, c_obj.data.bevel_factor_end = bounds
                    c_obj.parent = bpy.data.objects[molecule.name]

                    to_link.append(c_obj)
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
                curveOB.data.bevel_object = bpy.data.objects[bond.bevelname]
                curveOB.data.use_fill_caps = True
                curveOB.parent = bpy.data.objects[molecule.name]

                to_link.append(curveOB)

                to_hook.append((iatom.name, jatom.name, bondname))

    clock.tick_print("bond creation")

    # finalize by linking in objects and updating scene
    scn = context.scene
    for ob in to_link:  # link in objects
        scn.objects.link(ob)

    clock.tick_print("link objects")

    # update scene
    scn.update()

    clock.tick_print("update scene")

    if options["hook_atoms"]:
        for iatom, jatom, bond in to_hook:
            # Hook to atom1
            bpy.ops.object.select_all(action='DESELECT')
            bpy.data.objects[iatom].select = True
            bpy.data.objects[bond].select = True
            context.scene.objects.active = bpy.data.objects[bond]
            bpy.ops.object.mode_set(mode='EDIT')
            bpy.ops.curve.de_select_first()
            bpy.ops.object.hook_add_selob()
            bpy.ops.object.mode_set(mode='OBJECT')

            # Hook to atom2
            bpy.ops.object.select_all(action='DESELECT')
            bpy.data.objects[jatom].select = True
            bpy.data.objects[bond].select = True
            context.scene.objects.active = bpy.data.objects[bond]
            bpy.ops.object.mode_set(mode='EDIT')
            bpy.ops.curve.de_select_first()
            bpy.ops.curve.de_select_last()
            bpy.ops.object.hook_add_selob()
            bpy.ops.object.mode_set(mode='OBJECT')

        clock.tick_print("set hooks")

    for child, parent in to_parent:
        bpy.ops.object.select_all(action='DESELECT')
        child.select = True
        context.scene.objects.active = parent
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
                context.scene.objects.link(curveOB)
                context.scene.objects.active = curveOB
                curveOB.select = True
                context.scene.objects.active = bpy.data.objects[molecule.name]
                bpy.ops.object.parent_set(type='OBJECT', keep_transform=False)

    return


def PlotWireFrame(context, molecule, _options):
    """driver for wireframe plotter"""
    for item in context.selectable_objects:
        item.select = False

    verts = [atom.position for atom in molecule.atoms]
    bonds = [(b.iatom.index, b.jatom.index) for b in molecule.bonds]

    molec_mesh = bpy.data.meshes.new(molecule.name)
    molec_mesh.from_pydata(verts, bonds, [])
    molec_mesh.update()
    molec_obj = bpy.data.objects.new(molecule.name, molec_mesh)
    molec_obj.data = molec_mesh
    context.scene.objects.link(molec_obj)
    molec_obj.select = True
    context.scene.objects.active = bpy.data.objects[molecule.name]
    bpy.ops.object.convert(target='CURVE')
    context.object.data.fill_mode = 'FULL'
    context.object.data.render_resolution_u = 12
    context.object.data.bevel_depth = 0.1
    context.object.data.bevel_resolution = 4
    bpy.ops.object.shade_smooth()


def AnimateMolecule(context, molecule, options):
    """Keyframes position of atoms and opacity of bonds"""
    kstride = options["keystride"]

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
            item.select = False

        atom_obj = bpy.data.objects[atom.name]
        pluscharge = bpy.data.objects[atom.plus_charge] if (
            options["charges"] != "none") else ""
        negcharge = bpy.data.objects[atom.neg_charge] if (
            options["charges"] != "none") else ""

        for (iframe, isnap) in enumerate(atom.trajectory):
            atom_obj.location = isnap.position
            atom_obj.keyframe_insert(
                data_path='location', frame=iframe * kstride + 1)

            if options["charges"] != "none":
                if options["charges"] == "scale":
                    pc, nc = max(0.0, isnap.charge), -min(0.0, isnap.charge)
                    pluscharge.scale = molecule.scale(pc)
                    negcharge.scale = molecule.scale(nc)
                    for x in [pluscharge, negcharge]:
                        x.keyframe_insert(data_path='scale',
                                          frame=iframe * kstride + 1)

    if options["animate_bonds"] == "dynamic":
        bondmask = molecule.bond_mask(options)
        for bond in molecule.bonds:
            i, j = bond.iatom.index, bond.jatom.index
            for bond_obj in [ bpy.data.objects[b] for b in bond.name ]:
                for (iframe, mask) in enumerate(bondmask[(i, j)]):
                    bond_obj.hide = not mask
                    bond_obj.keyframe_insert(
                        data_path="hide", frame=iframe * kstride + 1)
                    bond_obj.hide_render = not mask
                    bond_obj.keyframe_insert(
                        data_path="hide_render", frame=iframe * kstride + 1)
    return


def create_mesh(name, verts, faces, material, context, remesh=True):
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

    if remesh:
        mod = ob.modifiers.new('Remesh', 'REMESH')
        mod.mode = 'SMOOTH'
        mod.octree_depth = 7
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
    vertex_sets = []
    wm = context.window_manager
    if molecule.volume is not None:  # volumetric data was read
        vol = molecule.volume
        isovals = options["isovalues"]

        # marching cubes to make the surface
        vsets = cube_isosurface(vol.data, vol.origin, vol.axes, isovals, wm)
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

            try:
                neworb = int(o)
            except ValueError:
                pass

            if neworb is None:
                raise Exception("Could not understand orbital list!")
            orblist.append(neworb)

        isovals = options["isovalues"]
        resolution = options["resolution"]

        for orbname, orbnumber in zip(orbnames, orblist):
            orb = orbitals.get_orbital(orbnumber)
            vset = molden_isosurface(orb, isovals, resolution, orbname, wm)
            vertex_sets.extend(vset)

    meshes = []
    materials = make_iso_materials(molecule, styler, vertex_sets, options)

    for v in vertex_sets:
        # add surfaces to the scene
        name = "{0}_{1}_{2}".format(molecule.name, v["name"], v["isovalue"])
        matname = "{0}_orbital_mat_{1}".format(v["name"], "plus" if v["isovalue"] > 0 else "minus")
        verts, faces = create_geometry(v["triangles"])
        meshes.append(create_mesh(name, verts, faces, materials[matname], context, options["remesh"]))

    mol_obj = bpy.data.objects[molecule.name]
    for m in meshes:
        context.scene.objects.link(m)
        m.select = True

    mol_obj.select = True
    # parent them
    context.scene.objects.active = mol_obj
    bpy.ops.object.parent_set(type='OBJECT', keep_transform=False)


@stopwatch("draw rings")
def plot_rings(context, molecule, options):
    """Given a list of atoms and bonds, determine where aromatic rings are and plot them"""
    planarCycles = molecule.rings

    ringmat = bpy.data.materials[molecule.materials["ring"]]

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
        context.scene.objects.link(ringObj)
        to_parent.append(ringObj)

    molecule_obj = bpy.data.objects[molecule.name]
    for child in to_parent:
        bpy.ops.object.select_all(action='DESELECT')
        child.select = True
        context.scene.objects.active = molecule_obj
        bpy.ops.object.parent_set(type="OBJECT", keep_transform=False)


def process_options(filename, options):
    """postprocess choices that might interact with each other"""

    # force defaults in case called outside of the blender UI
    defaults = {"bonds": True,
                "bond_thickness": 0.2,
                "hook_atoms": "auto",
                "plot_style": "sticks",
                "plot_type": "frame",
                "object_type": "mesh",
                "keystride": 2,
                "animate_bonds": "staticfirst",
                "universal_bonds": True,
                "ignore_hydrogen": True,
                "gradient": False,
                "charges": "none",
                "charges_offset": 0.0,
                "charges_factor": 1.0,
                "find_aromatic": False,
                "recycle_materials": True,
                "isovalues": "",
                "volume": "orbital",
                "orbital": 0,
                "remesh" : True,
                "resolution": 0.5,
                "colors": "default"
               }

    for d in defaults:
        if d not in options:
            options[d] = defaults[d]

    hooking = {"on": True, "off": False,
               "auto": options["plot_type"] == "animate"}
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


def BlendMolecule(context, filename, **options):
    """basic driver that calls the appropriate plot functions"""

    options = process_options(filename, options)
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
