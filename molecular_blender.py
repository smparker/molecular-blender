import periodictable as pt

# import pybel if it is available
try: import pybel
except ImportError: pass
import sys
import math
import bpy
import mathutils

class Atom():
    el = pt.element(1.0,(1.0,1.0,1.0),1.0,"","")
    position = (0.0, 0.0, 0.0)
    name = "Name"   #I think this will be useful in determining which atom corredsponds to which object for animations, it will be set in the PlotAtoms subroutine

    def __init__(self, symbol, position):
        self.el = pt.elements[symbol]
        self.position = position

#Read in xyz file and return array of atoms
def ImportXYZ(filename):
    out = []
    fh = open(filename,"r")
    raw = fh.readlines()
    fh.close()

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

#Given an array of atoms, returns set of unique elements
def FormBaseSet(atoms):
    out = set()
    for i in atoms:
        out.add(i.el.symbol)
    return out

#Given an array of atoms, returns list of atom pairs that are bonded
# Uses the average of the two vdw radii as a bond cutoff. This may not be ideal
def ComputeBonds(atoms):
    out = []

    natoms = len(atoms)
    for i in range(natoms):
        for j in range(i):
            ixyz = i.position
            jxyz = j.position
            vec = (ixyz[0] - jxyz[0], ixyz[1] - jxyz[1], ixyz[2] - jxyz[2])
            distance = math.sqrt(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2])
            if (distance <= (0.5*(atoms[i].el.vdw + atoms[j].el.vdw))):
                out.append( (i,j) )
    return out

#Given a set of strings containing all elements in the molecule, creates required materials
def MakeMaterials(atom_base_set):
    for atom in atom_base_set:
        bpy.data.materials.new(atom)                #Creates new material
        bpy.data.materials[atom].diffuse_color = mathutils.Color(pt.elements[atom].color) #Sets color from atom dictionary
    return

#Given list of types Atom(), plots in scene
def PlotAtoms(atom_list, objtype="mesh"):
    #Check to see if original atom already exists, if yes, create translated linked duplicate, if no, create new object
    for atom in atom_list:
        base_atom = "0".join(atom.el.symbol)
        if base_atom in bpy.data.objects.keys():
            #create name of new object
            sufx = 0
            atom_name = str(atom.el.symbol).join(str(sufx))
            while atom_name in bpy.data.objects.keys():
                sufx += 1
                atom_name = str(atom.el.symbol).join(str(sufx))
            atom.name = atom_name
            #Create the linked duplicate object
            bpy.context.scene.objects.active = bpy.data.objects[base_atom] #Set active object to base
            translation_vector = tuple([y-x for x,y in zip(atom.position,tuple(bpy.context.object.location))])
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
            bpy.context.object.name = atom_name
            bpy.context.object.data.name = base_atom
        else:   #Create the base atom from which all other of same element will be copied
            atom.name = base_atom
            bpy.ops.mesh.primitive_uv_sphere_add(location=atom.position, size=atom.el.vdw)
            bpy.context.object.name = base_atom
            bpy.context.object.data.name = base_atom
            bpy.context.object.data.materials.append(bpy.data.materials[atom.el.symbol])
            bpy.ops.object.shade_smooth()
    return
