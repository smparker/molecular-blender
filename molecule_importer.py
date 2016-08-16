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


from .periodictable import element,symbols,generate_table
from .find_planar_rings import plotRings

# import pybel if it is available
try: import pybel
except ImportError: pass
import sys
import math

from collections import namedtuple

import time
import bpy
import mathutils

elements = {} # container for elements

# decorator to measure time in a function using blender timer
def stopwatch(routine):
    def stopwatch_dec(func):
        def wrapper(*args, **kwargs):
            start = time.time()
            out = func(*args, **kwargs)
            end = time.time()
            print("%.4f sec elapsed in routine %s" % ((end-start), routine))
            return out
        return wrapper
    return stopwatch_dec

class Timer(object):
    def __init__(self):
        self.last_time = time.time()

    def tick(self):
        lt, self.last_time = self.last_time, time.time()
        return self.last_time - lt

    def tick_print(self, label):
        print("  %40s: %.4f sec" % (label, self.tick()))

class Snapshot(object):
    def __init__(self, position, charge, gradient):
        self.position = mathutils.Vector(position)
        self.charge = float(charge)
        self.gradient = mathutils.Vector(position)

class Atom(object):
    """Collect information on single atom"""
    def __init__(self, symbol, position, index, name = "", charge = 0.0, gradient = [0.0, 0.0, 0.0]):
        self.el = elements[symbol]
        self.position = mathutils.Vector(position)
        self.index = index
        self.name = name
        self.charge = charge
        self.gradient = mathutils.Vector(gradient)
        self.trajectory = []

class Bond(object):
    """Join two atoms in a bond"""
    def __init__(self, iatom, jatom, style = "", name = ""):
        self.iatom = iatom
        self.jatom = jatom
        self.name = name
        self.bevelname = ""
        if style == "":
            self.style = (jatom if iatom.el.vdw > jatom.el.vdw else iatom).el.symbol
        else:
            self.style = style
        self.threshold = 1.2*(iatom.el.covalent + jatom.el.covalent)

    def is_bonded(self, distance = None):
        dist = distance if (distance is not None) else (self.iatom.position - self.jatom.position).length
        return dist <= self.threshold

    def make_name(self, basename):
        """builds name for a molecule from basename and connected atoms"""
        iname = self.iatom.name.split('_')[-1]
        jname = self.jatom.name.split('_')[-1]
        return basename + "_" + iname + "-" + jname

class Molecule(object):
    """Atoms and bonds form a molecule"""
    def __init__(self, name, atoms):
        self.name = name
        self.atoms = atoms
        self.bonds = []
        self.materials = {}
        self.bond_materials = {}
        self.chgoff = 1.0
        self.chgfac = 1.0

    def COM(self):
        """Computes center of mass from atoms list"""
        total_mass = float(0.0)
        out = mathutils.Vector((0.0, 0.0, 0.0))
        for i in self.atoms:
            out += i.position * i.el.mass
            total_mass += i.el.mass
        out /= total_mass
        return out

    @stopwatch("determine_bonding")
    def determine_bonding(self, options):
        """builds list of bonded pairs"""
        atoms = self.atoms
        natoms = len(atoms)
        bondstyle = "universal" if options["universal_bonds"] else ""
        if (options["bonds"]):
            for i in range(natoms):
                for j in range(i):
                    # make a Bond object no matter what, then use it to check whether to keep
                    bond = Bond(atoms[i], atoms[j], bondstyle)
                    if (options["plot_type"] == "frame" or options["animate_bonds"] == "staticfirst"):
                        # only use the first frame to look for bonded pairs
                        vec = atoms[i].position - atoms[j].position
                        if bond.is_bonded(vec.length):
                            self.bonds.append(bond)
                    elif (options["animate_bonds"] in [ "staticall", "dynamic" ]):
                        # search through entire trajectory for a bond
                        for va, vb in zip(atoms[i].trajectory, atoms[j].trajectory):
                            if (bond.is_bonded((va.position - vb.position).length)):
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
                vec = va.position - vb.position
                pairmask.append(bond.is_bonded(vec.length))

            outmask[(iatom.index,jatom.index)] = pairmask
        return outmask

    def scale(self, chg):
        thr = 1.0 - self.chgoff
        out = self.chgfac * (chg - thr) + 1.0 if chg > thr else 0.0
        return mathutils.Vector((out,out,out))

@stopwatch("ImportXYZ")
def ImportXYZ(filename, options):
    """Read in xyz file and return list of atoms"""
    out = []
    fh = open(filename,"r")

    if (options["plot_type"] == "frame"):
        # first line contains number of atoms
        natoms = int(fh.readline().split()[0])
        # second line is a comment
        fh.readline()

        index = 0
        for line in fh:
            # Expecting:
            #   <symbol> <x> <y> <z> [<charge> [<vx> <vy> <vz>]]
            # e.g.
            #   h 0.0 0.0 1.0 0.0 1.0 2.0 3.0
            tmp = line.split()
            symb = str(tmp[0]).lower()
            position = [ float(x) for x in tmp[1:4] ]
            charge = float(tmp[4]) if len(tmp) >= 5 else 0.0
            vec = [ float(x) for x in tmp[5:8] ] if len(tmp) >= 8 else [ 0.0 for x in range(3) ]
            out.append(Atom(symb, position, index, charge=charge, gradient=vec))
            index += 1

        assert(index == natoms)
    elif (options["plot_type"] == "animate"):
        # first line contains number of atoms
        natoms = int(fh.readline().split()[0])
        # second line is a comment
        fh.readline()

        for iatom in range(natoms):
            tmp = fh.readline().split()
            symb = str(tmp[0]).lower()
            position = [ float(x) for x in tmp[1:4] ]
            charge = float(tmp[4]) if len(tmp) >= 5 else 0.0
            vec = [ float(x) for x in tmp[5:8] ] if len(tmp) >= 8 else [ 0.0 ] * 3
            new_atom = Atom(symb, position, iatom, charge=charge, gradient=vec)

            new_atom.trajectory.append(Snapshot(position, charge, vec))
            out.append(new_atom)

        while (True):
            line = fh.readline()
            if (line == ""): break
            frame_atoms = int(line.split()[0]) # natoms
            if (frame_atoms != natoms):
                raise Exception("All frames in trajectory must have the same number of atoms.")

            fh.readline() # comment line
            for i in range(natoms):
                tmp = fh.readline().split()
                symb = str(tmp[0]).lower()
                if (symb != out[i].el.symbol):
                    raise Exception("The order of the atoms must be the same for each frame in the animation.")
                position = [ float(x) for x in tmp[1:4] ]
                charge = float(tmp[4]) if len(tmp) >= 5 else 0.0
                vec = [ float(x) for x in tmp[5:8] ] if len(tmp) >= 8 else [ 0.0 for x in range(3) ]
                out[i].trajectory.append(Snapshot(position, charge, vec))

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

def plot_prep(molecule, options):
    if (options["charges"] != "none"):
        molecule.chgfac = options["charge_factor"]
        molecule.chgoff = options["charge_offset"]

@stopwatch("make_atom_materials")
def make_atom_materials(molecule, options):
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

    if (options["charges"] != "none"): # make materials used for showing charges
        pc = unique_name("pluscharge", bpy.data.materials.keys()) + "_mat"
        nc = unique_name("negcharge", bpy.data.materials.keys()) + "_mat"

        molecule.pluscharge_mat = pc
        molecule.negcharge_mat = nc

        for ch, col in [ (pc, (0,0,1)), (nc, (1,0,0)) ]:
            bpy.data.materials.new(ch)
            bpy.data.materials[ch].diffuse_color = mathutils.Color(col)
    return

@stopwatch("make_bond_materials")
def make_bond_materials(molecule, options):
    """Given a molecule object, creates the corresponding materials for the bonds"""

    # set of unique bond types
    bond_base_set = set([i.style for i in molecule.bonds])

    for b in bond_base_set:
        bond = b + "_bond_mat"
        if (options["recycle_materials"]):
            if bond in bpy.data.materials.keys():
                molecule.bond_materials[b] = bond
                continue

        bond = unique_name(bond, bpy.data.materials.keys()) # unique name for bond material
        molecule.bond_materials[b] = bond

        bpy.data.materials.new(bond)
        bpy.data.materials[bond].diffuse_color = mathutils.Color((0.9,0.9,0.9))
    return

@stopwatch("PlotMolecule")
def PlotMolecule(context, molecule, options):
    """Plots the given molecule object to the specified context"""
    global pt
    global namedtuple

    clock = Timer()

    object_type = options["object_type"]
    center_of_mass = molecule.COM()

    # Creates empty to collect all the new objects
    bpy.ops.object.empty_add(type='PLAIN_AXES',location=center_of_mass)

    # update molecule name to make sure it's unique
    molecule.name = unique_name(molecule.name, bpy.data.objects.keys())
    # set name for new parent
    context.object.name = molecule.name

    # hold onto reference to parent
    molecule_obj = bpy.data.objects[molecule.name]

    Style = namedtuple('Style', ['atom_size', 'bond_size', 'atom_scaling', 'bond_scaling'])
    style_dict = { 'vdw'       : Style( atom_size="scaled", bond_size="none", atom_scaling=1.0, bond_scaling=0.0 ),
                   'bs'        : Style( atom_size="scaled", bond_size="vdw", atom_scaling=0.25, bond_scaling=0.1 ),
                   'fixedbs'   : Style( atom_size="scaled", bond_size="fixed", atom_scaling=0.25, bond_scaling=1.0 ),
                   'sticks'     : Style( atom_size="bond", bond_size="fixed", atom_scaling=0.2, bond_scaling=1.0 ) }
    plot_style = style_dict[options["plot_style"]]

    bond_thickness = options["bond_thickness"]

    # local list of keys to keep making unique names
    obj_keys = bpy.data.objects.keys().copy()

    # list of objects that need to be processed after scene creation
    to_link = []
    to_hook = []
    to_parent = []

    clock.tick_print("preamble")

    # Check to see if original atom already exists, if yes, create translated linked duplicate, if no, create new object
    for atom in molecule.atoms:
        # Unselect Everything
        bpy.ops.object.select_all(action='DESELECT')
        base_atom = molecule.name + "_" + atom.el.name
        ref_atom = base_atom + "0"
        atom_obj = ""
        if ref_atom in obj_keys: # duplicate existing atom
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
            obj_keys.append(ref_atom)

            # hold onto reference
            atom_obj = context.object

            context.object.data.materials.append(bpy.data.materials[molecule.materials[atom.el.symbol]])
            bpy.ops.object.shade_smooth()

            # make parent active
            context.scene.objects.active = molecule_obj
            bpy.ops.object.parent_set(type='OBJECT',keep_transform=False)

        if (options["charges"] != "none"): # set a sphere on top of atom that will show charge
            plus_obj = atom_obj.copy()
            plus_obj.data = plus_obj.data.copy()
            plus_obj.data.materials[0] = bpy.data.materials[molecule.pluscharge_mat]
            plus_obj.name = atom.name + "_plus"
            to_link.append(plus_obj)

            neg_obj = atom_obj.copy()
            neg_obj.data = neg_obj.data.copy()
            neg_obj.data.materials[0] = bpy.data.materials[molecule.negcharge_mat]
            neg_obj.name = atom.name + "_neg"
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

    if len(molecule.bonds) > 0: # Add the bonds
        # Make circles the correct size for the bonds
        bevnames = {}

        bond_bevels = set([i.style for i in molecule.bonds])
        # set up bevel types
        for i in bond_bevels:
            rad = bond_thickness
            if (not options['universal_bonds'] and plot_style.bond_size == "vdw"):
                rad = elements[i].vdw * bond_scaling

            bpy.ops.curve.primitive_bezier_circle_add(radius=rad,location=(0,0,0))
            bevelname = unique_name(molecule.name + "_" + i + "_bond", bpy.data.objects.keys())
            bevnames[i] = bevelname
            context.object.name = bevelname

            context.scene.objects.active = molecule_obj
            bpy.ops.object.parent_set(type='OBJECT',keep_transform=False)

        make_bond_materials(molecule, options)

        for bond in molecule.bonds: # make curves for bonds
            # deselect all
            bpy.ops.object.select_all(action='DESELECT')

            iatom = bond.iatom
            jatom = bond.jatom

            bevtype = (iatom if iatom.el.vdw > jatom.el.vdw else jatom).el.symbol
            bond.bevelname = bevnames[bond.style]
            bond.name = unique_name(bond.make_name(molecule.name), obj_keys)
            obj_keys.append(bond.name)

            #create curve object
            coords = [iatom.position, jatom.position]
            curve = bpy.data.curves.new(bond.name,type='CURVE')
            curve.materials.append(bpy.data.materials[molecule.bond_materials[bond.style]])

            curve.dimensions = '3D'
            curve.resolution_u = 2
            bondline = curve.splines.new('BEZIER')
            bondline.bezier_points.add(len(coords)-1)
            for i, pnt in enumerate(coords):
                p = bondline.bezier_points[i]
                p.co = pnt - molecule_obj.location
                p.handle_right_type = p.handle_left_type = 'AUTO'
            curveOB = bpy.data.objects.new(bond.name, curve)
            curveOB.data.bevel_object = bpy.data.objects[bond.bevelname]
            curveOB.data.use_fill_caps = True
            curveOB.parent = bpy.data.objects[molecule.name]

            to_link.append(curveOB)

            to_hook.append((iatom.name, jatom.name, bond.name))

    clock.tick_print("bond creation")

    # finalize by linking in objects and updating scene
    scn = context.scene
    for ob in to_link: # link in objects
        scn.objects.link(ob)

    clock.tick_print("link objects")

    # update scene
    scn.update()

    clock.tick_print("update scene")

    # add hooks
    for iatom, jatom, bond in to_hook:
        # Hook to atom1
        bpy.ops.object.select_all(action='DESELECT')
        bpy.data.objects[iatom].select=True
        bpy.data.objects[bond].select=True
        context.scene.objects.active = bpy.data.objects[bond]
        bpy.ops.object.mode_set(mode='EDIT')
        bpy.ops.curve.de_select_first()
        bpy.ops.object.hook_add_selob()
        bpy.ops.object.mode_set(mode='OBJECT')

        # Hook to atom2
        bpy.ops.object.select_all(action='DESELECT')
        bpy.data.objects[jatom].select=True
        bpy.data.objects[bond].select=True
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

    if (options["gradient"]):
        for atom in molecule.atoms:
            if (atom.gradient.length > 1.0e-30):
                bpy.ops.object.select_all(action='DESELECT')
                curve = bpy.data.curves.new(atom.name+"_gradient", type='CURVE')
                curve.dimensions = '3D'
                curve.resolution_u = 2
                gradline = curve.splines.new('BEZIER')
                gradline.bezier_points.add(1)
                gradline.bezier_points[0].co = atom.position
                gradline.bezier_points[1].co = atom.position + atom.gradient
                for p in gradline.bezier_points:
                    p.handle_right_type = p.handle_left_type = 'AUTO'
                curveOB = bpy.data.objects.new(atom.name+"_gradient", curve)

                curveOB.data.use_fill_caps = True
                context.scene.objects.link(curveOB)
                context.scene.objects.active = curveOB
                curveOB.select = True
                context.scene.objects.active = bpy.data.objects[molecule.name]
                bpy.ops.object.parent_set(type='OBJECT', keep_transform=False)

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
        pluscharge = bpy.data.objects[atom.plus_charge] if (options["charges"]!="none") else ""
        negcharge  = bpy.data.objects[atom.neg_charge]  if (options["charges"]!="none") else ""

        for (iframe, isnap) in enumerate(atom.trajectory):
            atom_obj.location = isnap.position
            atom_obj.keyframe_insert(data_path='location', frame = iframe*kstride + 1)

            if (options["charges"] != "none"):
                if (options["charges"] == "scale"):
                    pc, nc = max(0.0, isnap.charge), -min(0.0, isnap.charge)
                    pluscharge.scale = molecule.scale(pc)
                    negcharge.scale = molecule.scale(nc)
                    for x in [ pluscharge, negcharge ]:
                        x.keyframe_insert(data_path='scale', frame=iframe*kstride+1)

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
    global elements # first redefine elements list
    elements = generate_table(options["colors"])
    name = filename.rsplit('.', 1)[0].rsplit('/')[-1]
    atoms = ImportXYZ(filename, options)
    molecule = Molecule(name, atoms)
    molecule.determine_bonding(options)

    make_atom_materials(molecule, options)
    plot_prep(molecule, options)

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
