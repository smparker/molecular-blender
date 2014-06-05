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

import molecular_blender.periodictable as pt

# import pybel if it is available
try: import pybel
except ImportError: pass
import sys
import math

from collections import namedtuple

import bpy
import mathutils

class Atom():
    el = pt.element(1.0,1.0,(1.0,1.0,1.0),1.0,"","")
    position = mathutils.Vector((0.0, 0.0, 0.0))
    name = "Name"   #I think this will be useful in determining which atom corredsponds to which object for animations, it will be set in the PlotAtoms subroutine
    trajectory = []

    def __init__(self, symbol, position):
        self.el = pt.elements[symbol]
        self.position = mathutils.Vector(position)
        self.name = ""
        self.trajectory = []

#Read in xyz file and return array of atoms
def ImportXYZ(filename, options):
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

        for line in raw:
            tmp = line.split()
            symb = str(tmp[0]).lower()
            position = ( float(tmp[1]), float(tmp[2]), float(tmp[3]) )
            out.append(Atom(symb, position))
    elif (options["plot_type"] == "animate"):
        natoms = int(raw[0])
        if (len(raw) % (natoms+2) != 0):
            raise Exception("Trajectory file has the wrong number of lines. Should be a multiple of natoms+2.")
        nframes = int(len(raw)/(natoms+2))

        raw.pop(0)
        raw.pop(0)
        for line in range(natoms):
            tmp = raw[0].split()
            symb = str(tmp[0]).lower()
            position = ( float(tmp[1]), float(tmp[2]), float(tmp[3]) )
            new_atom = Atom(symb, position)
            new_atom.trajectory.append(mathutils.Vector(position))
            out.append(new_atom)
            raw.pop(0)

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

#Read in pdb file
def ImportPDB(filename):
    return BabelImport(filename, "pdb")

#Use openbabel to import
def BabelImport(filename, filetype):
    if ('pybel' not in sys.modules.keys()):
        raise Exception("Importing files in " + filetype + " format requires PyBel")

    mol = pybel.readfile(filetype, filename).next()

    out = []
    for atom in mol:
        out.append(Atom(pt.symbols[atom.atomicnum].lower(), atom.coord))

    return out

#Compute center of mass
# returns tuple
def ComputeCOM(atoms):
    total_mass = float(0.0)
    out = mathutils.Vector((0.0, 0.0, 0.0))
    for i in atoms:
        out += i.position * i.el.mass
        total_mass += i.el.mass
    out /= total_mass
    return out

#Given an array of atoms, returns list of atom pairs that are bonded
# Uses the average of the two vdw radii as a bond cutoff. This may not be ideal
def ComputeBonds(atoms):
    out = []
    natoms = len(atoms)
    for i in range(natoms):
        for j in range(i):
            vec = atoms[i].position - atoms[j].position
            if (vec.length <= 1.2*(atoms[i].el.covalent + atoms[j].el.covalent)):
                out.append( (i,j) )
    return out

#Given a set of strings containing all elements in the molecule, creates required materials
def MakeMaterials(atoms):
    atom_base_set = set()
    for i in atoms:
        atom_base_set.add(i.el.symbol)

    for atom in atom_base_set:
        bpy.data.materials.new(atom)                #Creates new material
        bpy.data.materials[atom].diffuse_color = mathutils.Color(pt.elements[atom].color) #Sets color from atom dictionary
    return

#Given list of types Atom(), plots in scene
def PlotAtoms(context, atom_list, objtype="mesh"):
    #Check to see if original atom already exists, if yes, create translated linked duplicate, if no, create new object
    for atom in atom_list:
        #Unselect Everything
        for item in context.selectable_objects:
            item.select = False
        base_atom = ''.join([atom.el.name, "0"])
        if base_atom in bpy.data.objects.keys():
            #create name of new object
            sufx = 1
            while 1:
                atom_name = ''.join([atom.el.name, str(sufx)])
                if atom_name not in bpy.data.objects.keys():
                    break
                sufx += 1
            atom.name = atom_name
            #Create the linked duplicate object
            bpy.data.objects[base_atom].select = True #Set active object to base
            context.scene.objects.active = bpy.data.objects[base_atom]
            translation_vector = tuple([x-y for x,y in zip(atom.position,tuple(context.object.location))])
            bpy.ops.object.duplicate_move_linked(\
                OBJECT_OT_duplicate={"linked":False, "mode":'TRANSLATION'},\
                TRANSFORM_OT_translate={\
                    "value":translation_vector, \
                    "constraint_axis":(False, False, False), \
                    "constraint_orientation":'GLOBAL', \
                    "mirror":False, \
                    "proportional":'DISABLED', \
                    "proportional_edit_falloff":'SMOOTH', \
                    "proportional_size":1, \
                    "snap":False, \
                    "snap_target":'CLOSEST', \
                    "snap_point":(0, 0, 0), \
                    "snap_align":False, \
                    "snap_normal":(0, 0, 0), \
                    "texture_space":False, \
                    "remove_on_cancel":False, \
                    "release_confirm":False}\
                )

            context.object.name = atom_name
            context.object.data.name = atom_name
        else:   #Create the base atom from which all other of same element will be copied
            atom.name = base_atom
            for item in context.selectable_objects:
                item.select = False
            if objtype.lower() == "nurbs":
                bpy.ops.surface.primitive_nurbs_surface_sphere_add(radius=atom.el.vdw,location=atom.position)
            elif objtype.lower() == "meta":
                bpy.ops.object.metaball_add(type='BALL',radius=atom.el.vdw,location=atom.position)
            else:
                bpy.ops.mesh.primitive_uv_sphere_add(location=atom.position, size=atom.el.vdw)
            context.object.name = base_atom
            context.object.data.name = base_atom
            context.object.data.materials.append(bpy.data.materials[atom.el.symbol])
            bpy.ops.object.shade_smooth()
    return

#Given list of types Atom(), plots as a single molecule, all atoms as parent to an empty, with bonds
def PlotMolecule(context, atom_list, name, bonds, options):
    #Make parent
    global pt
    global namedtuple
    object_type = options["object_type"]
    center_of_mass = ComputeCOM(atom_list)
    bpy.ops.object.empty_add(type='PLAIN_AXES',location=center_of_mass)
    context.object.name = name

    Style = namedtuple('Style', ['atom_size', 'bond_size', 'atom_scaling', 'bond_scaling'])
    style_dict = { 'vdw'       : Style( atom_size="scaled", bond_size="none", atom_scaling=1.0, bond_scaling=0.0 ),
                   'bs'        : Style( atom_size="scaled", bond_size="vdw", atom_scaling=0.25, bond_scaling=0.1 ),
                   'fixedbs'   : Style( atom_size="scaled", bond_size="fixed", atom_scaling=0.25, bond_scaling=1.0 ),
                   'sticks'     : Style( atom_size="bond", bond_size="fixed", atom_scaling=0.2, bond_scaling=1.0 ) }
    plot_style = style_dict[options["plot_style"]]

    bond_thickness = options["bond_thickness"]

    #Check to see if original atom already exists, if yes, create translated linked duplicate, if no, create new object
    for atom in atom_list:
        #Unselect Everything
        bpy.ops.object.select_all(action='DESELECT')
        base_atom = ''.join([atom.el.name, "0"])
        if base_atom in bpy.data.objects.keys():
            #create name of new object
            sufx = 1
            while 1:
                atom_name = ''.join([atom.el.name, str(sufx)])
                if atom_name not in bpy.data.objects.keys():
                    break
                sufx += 1
            atom.name = atom_name
            #Create the linked duplicate object
            bpy.data.objects[base_atom].select = True #Set active object to base
            context.scene.objects.active = bpy.data.objects[base_atom]
            translation_vector = tuple([x-y for x,y in zip(atom.position,tuple(context.object.location))])
            bpy.ops.object.duplicate_move_linked(\
                OBJECT_OT_duplicate={"linked":False, "mode":'TRANSLATION'},\
                TRANSFORM_OT_translate={\
                    "value":translation_vector, \
                    "constraint_axis":(False, False, False), \
                    "constraint_orientation":'GLOBAL', \
                    "mirror":False, \
                    "proportional":'DISABLED', \
                    "proportional_edit_falloff":'SMOOTH', \
                    "proportional_size":1, \
                    "snap":False, \
                    "snap_target":'CLOSEST', \
                    "snap_point":(0, 0, 0), \
                    "snap_align":False, \
                    "snap_normal":(0, 0, 0), \
                    "texture_space":False, \
                    "remove_on_cancel":False, \
                    "release_confirm":False}\
                )

            context.object.name = atom_name
            context.object.data.name = atom_name
        else:   #Create the base atom from which all other of same element will be copied
            atom.name = base_atom
            bpy.ops.object.select_all(action='DESELECT')

            radius = plot_style.atom_scaling * atom.el.vdw if (plot_style.atom_size != "bond") else bond_thickness

            if object_type == "nurbs":
                bpy.ops.surface.primitive_nurbs_surface_sphere_add(radius=radius,location=atom.position)
            elif object_type == "meta":
                bpy.ops.object.metaball_add(type='BALL',radius=radius,location=atom.position)
            else:
                bpy.ops.mesh.primitive_uv_sphere_add(location=atom.position, size=radius)
            context.object.name = base_atom
            context.object.data.name = base_atom
            context.object.data.materials.append(bpy.data.materials[atom.el.symbol])
            bpy.ops.object.shade_smooth()
            #make parent active
            context.scene.objects.active=bpy.data.objects[name]
            bpy.ops.object.parent_set(type='OBJECT',keep_transform=False)
    #Add the bonds if able
    if len(bonds) > 0:
        #Make circles the correct size for the bonds
        bond_bevels = set()
        for i in atom_list:
            bond_bevels.add(i.el.symbol)
        for i in bond_bevels:
            rad = pt.elements[i].vdw if (plot_style.bond_size == "vdw") else bond_thickness
            rad *= plot_style.bond_scaling
            bpy.ops.curve.primitive_bezier_circle_add(radius=rad,location=(0,0,0))
            context.object.name = i + "_bond"
            context.scene.objects.active=bpy.data.objects[name]
            bpy.ops.object.parent_set(type='OBJECT',keep_transform=False)
        #make curves for bonds
        for (x,y) in bonds:
            #deselect all
            bpy.ops.object.select_all(action='DESELECT')
            A1 = atom_list[x]
            A2 = atom_list[y]
            bond_size = A2.el.symbol if (A1.el.vdw > A2.el.vdw) else A1.el.symbol
            bond_name = '-'.join([A1.name, A2.name])
            #create curve object
            coords = [A1.position,A2.position]
            curve = bpy.data.curves.new(bond_name,type='CURVE')
            curve.dimensions = '3D'
            curve.resolution_u = 2
            bondline = curve.splines.new('BEZIER')
            bondline.bezier_points.add(len(coords)-1)
            for i, pnt in enumerate(coords):
                p = bondline.bezier_points[i]
                p.co = pnt
                p.handle_right_type = p.handle_left_type = 'AUTO'
            curveOB = bpy.data.objects.new(bond_name, curve)
            curveOB.data.bevel_object = bpy.data.objects[bond_size + "_bond"]
            curveOB.data.use_fill_caps = True
            scn = context.scene
            scn.objects.link(curveOB)
            scn.objects.active = curveOB
            curveOB.select = True
            context.scene.objects.active=bpy.data.objects[name]
            bpy.ops.object.parent_set(type='OBJECT',keep_transform=False)
            #Hook to atom1
            bpy.ops.object.select_all(action='DESELECT')
            bpy.data.objects[A1.name].select=True
            bpy.data.objects[bond_name].select=True
            context.scene.objects.active = bpy.data.objects[bond_name]
            bpy.ops.object.mode_set(mode='EDIT')
            bpy.ops.curve.de_select_first()
            bpy.ops.object.hook_add_selob()
            bpy.ops.object.mode_set(mode='OBJECT')
            #Hook to atom2
            bpy.ops.object.select_all(action='DESELECT')
            bpy.data.objects[A2.name].select=True
            bpy.data.objects[bond_name].select=True
            context.scene.objects.active = bpy.data.objects[bond_name]
            bpy.ops.object.mode_set(mode='EDIT')
            bpy.ops.curve.de_select_first()
            bpy.ops.curve.de_select_last()
            bpy.ops.object.hook_add_selob()
            bpy.ops.object.mode_set(mode='OBJECT')
    return

def AnimateMolecule(context, atom_list, options):
    kstride = options["keystride"]
    # check to make sure trajectory information is stored for at least the first atom
    if (len(atom_list[0].trajectory) == 0):
        raise Exception("Trajectory must be read before calling AnimateMolecule")
    for atom in atom_list:
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
    return

def PlotBonds(context, bond_list):
    for (Atom1, Atom2) in bond_list:
        #Unselect everything first to be safe
        for item in context.selectable_objects:
            item.select = False
        bond_vector = tuple([x-y for x,y in zip(Atom2.position,Atom1.position)])
        bond_length = math.sqrt(math.pow(bond_vector[0],2)+math.pow(bond_vector[1],2)+math.pow(bond_vector[2],2))
        bond_center = tuple([0.5*(x+y) for x,y in zip(Atom2.position,Atom1.position)])
        bond_radius = 0.1*min([Atom1.el.vdw,Atom2.el.vdw])
        theta = math.acos(bond_vector[2]/bond_length)
        phi = math.atan2(bond_vector[1],bond_vector[0])
        bpy.ops.mesh.primitive_cylinder_add(radius=bond_radius,depth=bond_length,location=bond_center,rotation=(0,theta,phi))
        bond_name = '-'.join([Atom1.name, Atom2.name])
        context.object.name = bond_name
        context.object.data.name = bond_name

def BlendMolecule(context, filename, **options):
    name = filename.rsplit('.', 1)[0].rsplit('/')[-1]
    atoms = ImportXYZ(filename, options)
    MakeMaterials(atoms)
    bonds = ComputeBonds(atoms) if options["bonds"] else []
    if (options["object_type"] == "wireframe"):
        PlotWireFrame(context, atoms, name, bonds, options)
    else:
        PlotMolecule(context, atoms, name, bonds, options)
    if (options["plot_type"] == "animate"):
        AnimateMolecule(context, atoms, options)

def PlotWireFrame(context, atom_list, name, bonds, options):
    for item in context.selectable_objects:
                item.select = False
    verts = []
    for atom in atom_list:
        verts.append(atom.position)
    molec_mesh = bpy.data.meshes.new(name)
    molec_mesh.from_pydata(verts,bonds,[])
    molec_mesh.update()
    molec_obj = bpy.data.objects.new(name,molec_mesh)
    molec_obj.data = molec_mesh
    context.scene.objects.link(molec_obj)
    molec_obj.select = True
    context.scene.objects.active = bpy.data.objects[name]
    bpy.ops.object.convert(target='CURVE')
    context.object.data.fill_mode = 'FULL'
    context.object.data.render_resolution_u = 12
    context.object.data.bevel_depth = 0.1
    context.object.data.bevel_resolution = 4
    bpy.ops.object.shade_smooth()
