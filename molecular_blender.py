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
    position = mathutils.Vector((0.0, 0.0, 0.0))
    name = "Name"   #I think this will be useful in determining which atom corredsponds to which object for animations, it will be set in the PlotAtoms subroutine
    trajectory = []

    def __init__(self, symbol, position):
        self.el = pt.elements[symbol]
        self.position = mathutils.Vector(position)
        self.name = ""
        self.trajectory = []

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

# Reads a concatenated list of XYZ files into a trajectory
#   requires the same number of atoms in the same order for the entire trajectory
#   returns an atom list with the trajectory field filled
def ImportXYZTrajectory(filename):
    atom_list = []
    fh = open(filename, "r")
    raw = fh.readlines()
    fh.close()

    # number of lines in fh should be a multiple of (natoms+2)
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
        atom_list.append(new_atom)
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
            if (symb != atom_list[i].el.symbol):
                raise Exception("The order of the atoms must be the same for each frame in the animation.")
            position = ( float(tmp[1]), float(tmp[2]), float(tmp[3]) )
            atom_list[i].trajectory.append(mathutils.Vector(position))
            raw.pop(0)

    return atom_list

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
            ixyz = atoms[i].position
            jxyz = atoms[j].position
            vec = (ixyz[0] - jxyz[0], ixyz[1] - jxyz[1], ixyz[2] - jxyz[2])
            distance = math.sqrt(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2])
            if (distance <= (0.5*(atoms[i].el.vdw + atoms[j].el.vdw))):
                out.append( (i,j) )
    return out

#Given a set of strings containing all elements in the molecule, creates required materials
def MakeMaterials(atoms):
    atom_base_set = FormBaseSet(atoms)
    for atom in atom_base_set:
        bpy.data.materials.new(atom)                #Creates new material
        bpy.data.materials[atom].diffuse_color = mathutils.Color(pt.elements[atom].color) #Sets color from atom dictionary
    return

#Given list of types Atom(), plots in scene
def PlotAtoms(atom_list, objtype="mesh"):
    #Check to see if original atom already exists, if yes, create translated linked duplicate, if no, create new object
    for atom in atom_list:
        #Unselect Everything
        for item in bpy.context.selectable_objects:
            item.select = False
        base_atom = ''.join([atom.el.symbol, "0"])
        if base_atom in bpy.data.objects.keys():
            #create name of new object
            sufx = 1
            while 1:
                atom_name = ''.join([atom.el.symbol, str(sufx)])
                if atom_name not in bpy.data.objects.keys():
                    break
                sufx += 1
            atom.name = atom_name
            #Create the linked duplicate object
            bpy.data.objects[base_atom].select = True #Set active object to base
            bpy.context.scene.objects.active = bpy.data.objects[base_atom]
            translation_vector = tuple([x-y for x,y in zip(atom.position,tuple(bpy.context.object.location))])
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
            bpy.context.object.data.name = atom_name
        else:   #Create the base atom from which all other of same element will be copied
            atom.name = base_atom
            for item in bpy.context.selectable_objects:
                item.select = False
            if objtype.lower() == "nurbs":
                bpy.ops.surface.primitive_nurbs_surface_sphere_add(radius=atom.el.vdw,location=atom.position)
            elif objtype.lower() == "metaballs":
                bpy.ops.object.metaball_add(type='BALL',radius=atom.el.vdw,location=atom.position)
            else:
                bpy.ops.mesh.primitive_uv_sphere_add(location=atom.position, size=atom.el.vdw)
            bpy.context.object.name = base_atom
            bpy.context.object.data.name = base_atom
            bpy.context.object.data.materials.append(bpy.data.materials[atom.el.symbol])
            bpy.ops.object.shade_smooth()
    return

#Given list of types Atom(), plots as a single molecule, all atoms as parent to an empty, with bonds
def PlotMolecule(atom_list, objtype="mesh", name="molecule", bonds=[]):
    #Make parent
    global pt
    center_of_mass = ComputeCOM(atom_list)
    bpy.ops.object.empty_add(type='PLAIN_AXES',location=center_of_mass)
    bpy.context.object.name = name
    #Check to see if original atom already exists, if yes, create translated linked duplicate, if no, create new object
    for atom in atom_list:
        if len(bonds) > 0:
            scale_factor = 0.2
        else:
            scale_factor = 1.0
        #Unselect Everything
        for item in bpy.context.selectable_objects:
            item.select = False
        base_atom = ''.join([atom.el.symbol, "0"])
        if base_atom in bpy.data.objects.keys():
            #create name of new object
            sufx = 1
            while 1:
                atom_name = ''.join([atom.el.symbol, str(sufx)])
                if atom_name not in bpy.data.objects.keys():
                    break
                sufx += 1
            atom.name = atom_name
            #Create the linked duplicate object
            bpy.data.objects[base_atom].select = True #Set active object to base
            bpy.context.scene.objects.active = bpy.data.objects[base_atom]
            translation_vector = tuple([x-y for x,y in zip(atom.position,tuple(bpy.context.object.location))])
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
            bpy.context.object.data.name = atom_name
        else:   #Create the base atom from which all other of same element will be copied
            atom.name = base_atom
            for item in bpy.context.selectable_objects:
                item.select = False
            if objtype.lower() == "nurbs":
                bpy.ops.surface.primitive_nurbs_surface_sphere_add(radius=scale_factor*atom.el.vdw,location=atom.position)
            elif objtype.lower() == "metaballs":
                bpy.ops.object.metaball_add(type='BALL',radius=scale_factor*atom.el.vdw,location=atom.position)
            else:
                bpy.ops.mesh.primitive_uv_sphere_add(location=atom.position, size=scale_factor*atom.el.vdw)
            bpy.context.object.name = base_atom
            bpy.context.object.data.name = base_atom
            bpy.context.object.data.materials.append(bpy.data.materials[atom.el.symbol])
            bpy.ops.object.shade_smooth()
            #make parent active
            bpy.context.scene.objects.active=bpy.data.objects[name]
            bpy.ops.object.parent_set(type='OBJECT',keep_transform=False)
    #Add the bonds if able
    if len(bonds) > 0:
        #Make circles the correct size for the bonds
        bond_bevels = set()
        for i in atom_list:
            bond_bevels.add(i.el.symbol)
        for i in bond_bevels:
            rad = 0.1*pt.elements[i].vdw
            bpy.ops.curve.primitive_bezier_circle_add(radius=rad,location=(0,0,0))
            bpy.context.object.name = i + "_bond"
            bpy.context.scene.objects.active=bpy.data.objects[name]
            bpy.ops.object.parent_set(type='OBJECT',keep_transform=False)
        #make curves for bonds
        for (x,y) in bonds:
            #deselect all
            for item in bpy.context.selectable_objects:
                item.select = False
            A1 = atom_list[x]
            A2 = atom_list[y]
            if A1.el.vdw > A2.el.vdw:
                bond_size = A2.el.symbol
            else:
                bond_size = A1.el.symbol
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
            scn = bpy.context.scene
            scn.objects.link(curveOB)
            scn.objects.active = curveOB
            curveOB.select = True
            bpy.context.scene.objects.active=bpy.data.objects[name]
            bpy.ops.object.parent_set(type='OBJECT',keep_transform=False)
            #Hook to atom1
            for item in bpy.context.selectable_objects:
                item.select = False
            bpy.data.objects[A1.name].select=True
            bpy.data.objects[bond_name].select=True
            bpy.context.scene.objects.active = bpy.data.objects[bond_name]
            bpy.ops.object.mode_set(mode='EDIT')
            bpy.ops.curve.de_select_first()
            bpy.ops.object.hook_add_selob()
            bpy.ops.object.mode_set(mode='OBJECT')
            #Hook to atom2
            for item in bpy.context.selectable_objects:
                item.select = False
            bpy.data.objects[A2.name].select=True
            bpy.data.objects[bond_name].select=True
            bpy.context.scene.objects.active = bpy.data.objects[bond_name]
            bpy.ops.object.mode_set(mode='EDIT')
            bpy.ops.curve.de_select_first()
            bpy.ops.curve.de_select_last()
            bpy.ops.object.hook_add_selob()
            bpy.ops.object.mode_set(mode='OBJECT')
    return

def AnimateMolecule(atom_list, kstride = 5):
    # check to make sure trajectory information is stored for at least the first atom
    if (len(atom_list[0].trajectory) == 0):
        raise Exception("Trajectory must be read before calling AnimateMolecule")
    for atom in atom_list:
        if (atom.name == ''):
            raise Exception("No name found for an atom. Was PlotMolecule or PlotAtoms called properly?")
        if atom.name not in bpy.data.objects.keys():
            raise Exception("Atom " + atom.name + " not found. Was PlotMolecule or PlotAtoms called properly?")

        # deselect everything for good measure
        for item in bpy.context.selectable_objects:
            item.select = False

        atom_obj = bpy.data.objects[atom.name]

        for (iframe, position) in enumerate(atom.trajectory):
            atom_obj.location = position
            atom_obj.keyframe_insert(data_path='location', frame = iframe*kstride + 1)
    return

'''
import sys
sys.path.append("/Users/joshuaszekely/Desktop/Codes/Molecular-Blender")
import molecular_blender as mb
Atoms = mb.ImportXYZ("/Users/joshuaszekely/Desktop/Benzenes.xyz")
Base = mb.FormBaseSet(Atoms)
mb.MakeMaterials(Base)
Bonds=mb.ComputeBonds(Atoms)
mb.PlotMolecule(Atoms,name="Benzenes",bonds=Bonds)
'''

def PlotSingleBond(Atom1, Atom2):
    #Unselect everything first to be safe
    for item in bpy.context.selectable_objects:
        item.select = False
    bond_vector = tuple([x-y for x,y in zip(Atom2.position,Atom1.position)])
    bond_length = math.sqrt(math.pow(bond_vector[0],2)+math.pow(bond_vector[1],2)+math.pow(bond_vector[2],2))
    bond_center = tuple([0.5*(x+y) for x,y in zip(Atom2.position,Atom1.position)])
    bond_radius = 0.1*min([Atom1.el.vdw,Atom2.el.vdw])
    theta = math.acos(bond_vector[2]/bond_length)
    phi = math.atan2(bond_vector[1],bond_vector[0])
    bpy.ops.mesh.primitive_cylinder_add(radius=bond_radius,depth=bond_length,location=bond_center,rotation=(0,theta,phi))
    bond_name = '-'.join([Atom1.name, Atom2.name])
    bpy.context.object.name = bond_name
    bpy.context.object.data.name = bond_name

def PlotBonds(atoms, bond_list):
    for (x,y) in bond_list:
        PlotSingleBond(atoms[x],atoms[y])

def MoleculeFromFile(filename):
    atoms = ImportXYZ(filename)
    MakeMaterials(atoms)
    bonds = ComputeBonds(atoms)
    name = filename.rsplit('.', 1)[0].rsplit('/')[-1]
    PlotMolecule(atoms, name=name, bonds=bonds)

def AnimateFromFile(filename):
    atoms = ImportXYZTrajectory(filename)
    MakeMaterials(atoms)
    bonds = ComputeBonds(atoms)
    name = filename.rsplit('.', 1)[0].rsplit('/')[-1]
    PlotMolecule(atoms, name=name, bonds=bonds)
    AnimateMolecule(atoms, kstride=1)
