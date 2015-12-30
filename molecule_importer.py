#
#  Molecular Blender
#  Filename: molecular_blender.py
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


from .periodictable import elements,element,symbols
from .find_planar_rings import plotRings

# import pybel if it is available
try: import pybel
except ImportError: pass
import sys
import math

from collections import namedtuple

import bpy
import mathutils

class Atom():
    """Collect information on single atom"""
    def __init__(self, symbol, position, index, name = ""):
        self.el = elements[symbol]
        self.position = mathutils.Vector(position)
        self.index = index
        self.name = name
        self.trajectory = []

class Bond():
    """Join two atoms in a bond"""
    def __init__(self, iatom, jatom, name = "", bevel = ""):
        self.iatom = iatom
        self.jatom = jatom
        self.name = name
        self.bevel = bevel
        self.threshold = 1.2*(iatom.el.covalent + jatom.el.covalent)

    def is_bonded(self, distance = None):
        dist = distance if (distance is not None) else (self.iatom.position - self.jatom.position).length
        return dist <= self.threshold

    def make_name(self, basename):
        """builds name for a molecule from basename and connected atoms"""
        iname = self.iatom.name.split('_')[-1]
        jname = self.jatom.name.split('_')[-1]
        return basename + "_" + iname + "-" + jname

class Molecule():
    """Atoms and bonds form a molecule"""
    def __init__(self, name, atoms):
        self.name = name
        self.atoms = atoms
        self.bonds = []
        self.materials = {}

    def COM(self):
        """Computes center of mass from atoms list"""
        total_mass = float(0.0)
        out = mathutils.Vector((0.0, 0.0, 0.0))
        for i in self.atoms:
            out += i.position * i.el.mass
            total_mass += i.el.mass
        out /= total_mass
        return out

    def determine_bonding(self, options):
        """builds list of bonded pairs"""
        atoms = self.atoms
        natoms = len(atoms)
        if (options["bonds"]):
            for i in range(natoms):
                for j in range(i):
                    # make a Bond object no matter what, then use it to check whether to keep
                    bond = Bond(atoms[i], atoms[j])
                    if (options["plot_type"] == "frame" or options["animate_bonds"] == "staticfirst"):
                        # only use the first frame to look for bonded pairs
                        vec = atoms[i].position - atoms[j].position
                        if bond.is_bonded(vec.length):
                            self.bonds.append(bond)
                    elif (options["animate_bonds"] in [ "staticall", "dynamic" ]):
                        # search through entire trajectory for a bond
                        for va, vb in zip(atoms[i].trajectory, atoms[j].trajectory):
                            if (bond.is_bonded((va - vb).length)):
                                self.bonds.append(bond)
                                break


    def bond_mask(self, options):
        """Construct mask determining whether a bond should be drawn at each frame"""
        outmask = { }
        natoms = len(self.atoms)
        for bond in self.bonds:
            iatom, jatom = bond.iatom, bond.jatom

            pairmask = []
            for va, vb in zip(iatom.trajectory, jatom.trajectory):
                vec = va - vb
                pairmask.append(bond.is_bonded(vec.length))

            outmask[(iatom.index,jatom.index)] = pairmask
        return outmask

def ImportXYZ(filename, options):
    """Read in xyz file and return list of atoms"""
    out = []
    fh = open(filename,"r")
    raw = fh.readlines()
    fh.close()

    if (options["plot_type"] == "frame"):
        natoms = int(raw[0])
        raw.pop(0) # first line contains number of atoms
        raw.pop(0) # second line is a comment
        if (natoms != len(raw)):
            raise Exception("Improperly formatted xyz file!")

        index = 0
        for line in raw:
            tmp = line.split()
            symb = str(tmp[0]).lower()
            position = ( float(tmp[1]), float(tmp[2]), float(tmp[3]) )
            out.append(Atom(symb, position, index))
            index += 1
        assert(index == natoms)
    elif (options["plot_type"] == "animate"):
        natoms = int(raw[0])
        if (len(raw) % (natoms+2) != 0):
            raise Exception("Trajectory file has the wrong number of lines. Should be a multiple of natoms+2.")
        nframes = int(len(raw)/(natoms+2))

        raw.pop(0)
        raw.pop(0)
        index = 0
        for line in range(natoms):
            tmp = raw[0].split()
            symb = str(tmp[0]).lower()
            position = ( float(tmp[1]), float(tmp[2]), float(tmp[3]) )
            new_atom = Atom(symb, position, index)
            index += 1
            new_atom.trajectory.append(mathutils.Vector(position))
            out.append(new_atom)
            raw.pop(0)

        assert(index == natoms)

        for ifrm in range(nframes-1):
            frame_atoms = int(raw[0])
            if (frame_atoms != natoms):
                raise Exception("All frames in trajectory must have the same number of atoms.")

            raw.pop(0)
            raw.pop(0)
            for i in range(natoms):
                tmp = raw[0].split()
                symb = str(tmp[0]).lower()
                if (symb != out[i].el.symbol):
                    raise Exception("The order of the atoms must be the same for each frame in the animation.")
                position = ( float(tmp[1]), float(tmp[2]), float(tmp[3]) )
                out[i].trajectory.append(mathutils.Vector(position))
                raw.pop(0)

    return out

def ImportPDB(filename):
    """Read PDB using Babel, if available"""
    return BabelImport(filename, "pdb")

def BabelImport(filename, filetype):
    """Implementation for reading a pdb"""
    if ('pybel' not in sys.modules.keys()):
        raise Exception("Importing files in " + filetype + " format requires PyBel")

    mol = pybel.readfile(filetype, filename).next()

    out = []
    for atom in mol:
        out.append(Atom(symbols[atom.atomicnum].lower(), atom.coord))

    return out

def unique_name(name, existing_names, starting_suffix = None):
    """If name is not in existing_names, returns name. Otherwise, returns name + "0", 1, 2, etc."""
    testname = name if starting_suffix is None else "%s%d" % (name, starting_suffix)
    if testname in existing_names:
        n = 0 if starting_suffix is None else starting_suffix+1
        while (True):
            testname = "%s%d" % (name, n)
            if testname in existing_names:
                n += 1
            else:
                return testname
    else:
        return testname

def MakeMaterials(molecule, options):
    """Given a molecule object, creates required materials and populates materials dictionary in molecule"""

    # set of new atoms to make materials for
    atom_base_set = set([ i.el.symbol for i in molecule.atoms])

    for a in atom_base_set:
        atom = a + "_mat"
        if (options["recycle_materials"]):
            if atom in bpy.data.materials.keys():
                molecule.materials[a] = atom
                continue

        atom = unique_name(atom, bpy.data.materials.keys()) # get unique name for new material
        molecule.materials[a] = atom # store unique materials

        bpy.data.materials.new(atom) # Creates new material
        bpy.data.materials[atom].diffuse_color = mathutils.Color(elements[a].color) # Sets color from atom dictionary
    return

def PlotMolecule(context, molecule, options):
    """Plots the given molecule object to the specified context"""
    global pt
    global namedtuple

    object_type = options["object_type"]
    center_of_mass = molecule.COM()

    # Creates empty to collect all the new objects
    bpy.ops.object.empty_add(type='PLAIN_AXES',location=center_of_mass)

    # update molecule name to make sure it's unique
    molecule.name = unique_name(molecule.name, bpy.data.objects.keys())

    # set name for new parent
    context.object.name = molecule.name

    Style = namedtuple('Style', ['atom_size', 'bond_size', 'atom_scaling', 'bond_scaling'])
    style_dict = { 'vdw'       : Style( atom_size="scaled", bond_size="none", atom_scaling=1.0, bond_scaling=0.0 ),
                   'bs'        : Style( atom_size="scaled", bond_size="vdw", atom_scaling=0.25, bond_scaling=0.1 ),
                   'fixedbs'   : Style( atom_size="scaled", bond_size="fixed", atom_scaling=0.25, bond_scaling=1.0 ),
                   'sticks'     : Style( atom_size="bond", bond_size="fixed", atom_scaling=0.2, bond_scaling=1.0 ) }
    plot_style = style_dict[options["plot_style"]]

    bond_thickness = options["bond_thickness"]

    # Check to see if original atom already exists, if yes, create translated linked duplicate, if no, create new object
    for atom in molecule.atoms:
        #Unselect Everything
        bpy.ops.object.select_all(action='DESELECT')
        base_atom = molecule.name + "_" + atom.el.name
        ref_atom = base_atom + "0"
        if ref_atom in bpy.data.objects.keys():
            # name new object
            atom.name = unique_name(base_atom, bpy.data.objects.keys(), 0)

            # Create the linked duplicate object
            bpy.data.objects[ref_atom].select = True #Set active object to base
            context.scene.objects.active = bpy.data.objects[ref_atom]
            translation_vector = atom.position - context.object.location
            bpy.ops.object.duplicate_move_linked(
                OBJECT_OT_duplicate={"linked":False, "mode":'TRANSLATION'},
                TRANSFORM_OT_translate={
                    "value":translation_vector,
                    "constraint_axis":(False, False, False),
                    "constraint_orientation":'GLOBAL',
                    "mirror":False,
                    "proportional":'DISABLED',
                    "proportional_edit_falloff":'SMOOTH',
                    "proportional_size":1,
                    "snap":False,
                    "snap_target":'CLOSEST',
                    "snap_point":(0, 0, 0),
                    "snap_align":False,
                    "snap_normal":(0, 0, 0),
                    "texture_space":False,
                    "remove_on_cancel":False,
                    "release_confirm":False}
                )

            context.object.name = atom.name
            context.object.data.name = atom.name
        else:   #Create the base atom from which all other of same element will be copied
            atom.name = ref_atom
            bpy.ops.object.select_all(action='DESELECT')

            radius = plot_style.atom_scaling * atom.el.vdw if (plot_style.atom_size != "bond") else bond_thickness

            if object_type == "nurbs":
                bpy.ops.surface.primitive_nurbs_surface_sphere_add(radius=radius,location=atom.position)
            elif object_type == "meta":
                bpy.ops.object.metaball_add(type='BALL',radius=radius,location=atom.position)
            else:
                bpy.ops.mesh.primitive_uv_sphere_add(location=atom.position, size=radius)
            context.object.name = ref_atom
            context.object.data.name = ref_atom
            context.object.data.materials.append(bpy.data.materials[molecule.materials[atom.el.symbol]])
            bpy.ops.object.shade_smooth()

            #make parent active
            context.scene.objects.active=bpy.data.objects[molecule.name]
            bpy.ops.object.parent_set(type='OBJECT',keep_transform=False)

    if len(molecule.bonds) > 0: # Add the bonds
        # Make circles the correct size for the bonds
        bond_bevels = set([i.el.symbol for i in molecule.atoms])
        bevnames = {}
        for i in bond_bevels:
            rad = elements[i].vdw if (plot_style.bond_size == "vdw") else bond_thickness
            rad *= plot_style.bond_scaling
            bpy.ops.curve.primitive_bezier_circle_add(radius=rad,location=(0,0,0))
            bevelname = unique_name(molecule.name + "_" + i + "_bond", bpy.data.objects.keys())
            bevnames[i] = bevelname
            context.object.name = bevelname
            context.scene.objects.active=bpy.data.objects[molecule.name]
            bpy.ops.object.parent_set(type='OBJECT',keep_transform=False)

        for bond in molecule.bonds: # make curves for bonds
            #deselect all
            bpy.ops.object.select_all(action='DESELECT')

            iatom = bond.iatom
            jatom = bond.jatom

            bevtype = (iatom if iatom.el.vdw > jatom.el.vdw else jatom).el.symbol
            bond.bevel = bevnames[bevtype]
            bond.name = unique_name(bond.make_name(molecule.name), bpy.data.objects.keys())

            #create curve object
            coords = [iatom.position, jatom.position]
            curve = bpy.data.curves.new(bond.name,type='CURVE')

            curve.dimensions = '3D'
            curve.resolution_u = 2
            bondline = curve.splines.new('BEZIER')
            bondline.bezier_points.add(len(coords)-1)
            for i, pnt in enumerate(coords):
                p = bondline.bezier_points[i]
                p.co = pnt
                p.handle_right_type = p.handle_left_type = 'AUTO'
            curveOB = bpy.data.objects.new(bond.name, curve)
            curveOB.data.bevel_object = bpy.data.objects[bond.bevel]
            curveOB.data.use_fill_caps = True
            scn = context.scene
            scn.objects.link(curveOB)
            scn.objects.active = curveOB
            curveOB.select = True
            context.scene.objects.active=bpy.data.objects[molecule.name]
            bpy.ops.object.parent_set(type='OBJECT',keep_transform=False)

            #Hook to atom1
            bpy.ops.object.select_all(action='DESELECT')
            bpy.data.objects[iatom.name].select=True
            bpy.data.objects[bond.name].select=True
            context.scene.objects.active = bpy.data.objects[bond.name]
            bpy.ops.object.mode_set(mode='EDIT')
            bpy.ops.curve.de_select_first()
            bpy.ops.object.hook_add_selob()
            bpy.ops.object.mode_set(mode='OBJECT')

            #Hook to atom2
            bpy.ops.object.select_all(action='DESELECT')
            bpy.data.objects[jatom.name].select=True
            bpy.data.objects[bond.name].select=True
            context.scene.objects.active = bpy.data.objects[bond.name]
            bpy.ops.object.mode_set(mode='EDIT')
            bpy.ops.curve.de_select_first()
            bpy.ops.curve.de_select_last()
            bpy.ops.object.hook_add_selob()
            bpy.ops.object.mode_set(mode='OBJECT')
    return

def PlotWireFrame(context, molecule, options):
    """driver for wireframe plotter"""
    for item in context.selectable_objects:
        item.select = False

    verts = [atom.position for atom in molecule.atoms]
    bonds = [(b.iatom.index, b.jatom.index) for b in molecule.bonds]

    molec_mesh = bpy.data.meshes.new(molecule.name)
    molec_mesh.from_pydata(verts,bonds,[])
    molec_mesh.update()
    molec_obj = bpy.data.objects.new(molecule.name,molec_mesh)
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
    if (len(molecule.atoms[0].trajectory) == 0):
        raise Exception("Trajectory must be read before calling AnimateMolecule")
    for atom in molecule.atoms:
        if (atom.name == ''):
            raise Exception("No name found for an atom. Was PlotMolecule or PlotAtoms called properly?")
        if atom.name not in bpy.data.objects.keys():
            raise Exception("Atom " + atom.name + " not found. Was PlotMolecule or PlotAtoms called properly?")

        # deselect everything for good measure
        for item in context.selectable_objects:
            item.select = False

        atom_obj = bpy.data.objects[atom.name]

        for (iframe, position) in enumerate(atom.trajectory):
            atom_obj.location = position
            atom_obj.keyframe_insert(data_path='location', frame = iframe*kstride + 1)

    if options["animate_bonds"] == "dynamic":
        bondmask = molecule.bond_mask(options)
        for bond in molecule.bonds:
            i, j = bond.iatom.index, bond.jatom.index
            bond_obj = bpy.data.objects[bond.name]
            for (iframe, mask) in enumerate(bondmask[(i,j)]):
                bond_obj.hide = not mask
                bond_obj.keyframe_insert(data_path="hide", frame=iframe*kstride+1)
                bond_obj.hide_render = not mask
                bond_obj.keyframe_insert(data_path="hide_render", frame=iframe*kstride+1)
    return

def BlendMolecule(context, filename, **options):
    """basic driver that calls the appropriate plot functions"""
    name = filename.rsplit('.', 1)[0].rsplit('/')[-1]
    atoms = ImportXYZ(filename, options)
    molecule = Molecule(name, atoms)
    molecule.determine_bonding(options)

    MakeMaterials(molecule, options)
    if (options["object_type"] == "wireframe"):
        PlotWireFrame(context, molecule, options)
        if (options["plot_type"] == "animate"):
            print("animation of wireframes not currently supported!")
    else:
        PlotMolecule(context, molecule, options)
        if (options["plot_type"] == "animate"):
            AnimateMolecule(context, molecule, options)

    if (options["find_aromatic"]):
        plotRings(context, molecule, options)
