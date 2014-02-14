import periodictable as pt

import sys

class Atom():
    el = pt.element()
    position = (0.0, 0.0, 0.0)
    name = "Name"	#I think this will be useful in determining which atom corredsponds to which object for animations, it will be set in the PlotAtoms subroutine

    def __init__(self, symbol, position):
        self.el = pt.elements[symbol]
        self.position = position

#Given a set of strings containing all elements in the molecule, creates required materials
def MakeMaterials(atom_base_set):
	for atom in atom_base_set: 
		bpy.data.materials.new(atom)				#Creates new material
		bpy.data.materials[atom].diffuse_color = Color(pt.elements[atom].color)	#Sets color from atom dictionary
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
			while atom_name in bpy.data.objects.keys()
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
		else: 	#Create the base atom from which all other of same element will be copied
			atom.name = base_atom
			bpy.ops.mesh.primitive_uv_sphere_add(location=atom.position, size=atom.el.vdw)
			bpy.context.object.name = base_atom
			bpy.context.object.data.name = base_atom
			bpy.context.object.data.materials.append(bpy.data.materials[atom.el.symbol])
			bpy.ops.object.shade_smooth()
	return
