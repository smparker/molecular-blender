#
#  Molecular Blender
#  Filename: __init__.py
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

bl_info = {
    "name" : "Molecular Blender",
    "author" : "Shane Parker and Josh Szekely",
    "version" : (0,1,0),
    "location": "File > Import",
    "description" : "Import XYZ files",
    "category" : "Import-Export"
}

import bpy

from bpy_extras.io_utils import ImportHelper
from bpy.props import StringProperty, BoolProperty, EnumProperty, IntProperty, FloatProperty
from bpy.types import Operator

from .molecule_importer import BlendMolecule

class MolecularBlender(Operator, ImportHelper):
    """A class designed to conveniently import molecules to blender"""
    bl_idname = "import.molecular_blender"
    bl_label  = "Molecular Blender"

    # in the template, we'll figure it out later
    filename_ext = ".xyz"
    filter_glob  = StringProperty(
            default = "*.xyz",
            options = {'HIDDEN'},
          )

    type_of_object = EnumProperty(
                         name        = "Object type",
                         description = "Object type to use to draw atoms",
                         items       = (('meta', "Metaballs", "Metaballs"),
                                  ('nurbs', "NURBS", "NURBS"),
                                  ('mesh', "Mesh", "Mesh"),
                                  ('wireframe', "Wireframe", "Wireframe")),
                         default     = 'mesh'
                        )

    type_of_plot = EnumProperty(
                    name        = "Type",
                    description = "What to draw",
                    items       = (('frame', "Single Frame", "Plot a single frame from an XYZ"),
                             ('animate', "Animation", "Animate from an XYZ")),
                    default     = 'frame',
                   )

    style_of_plot = EnumProperty(
                     name        = "Style",
                     description = "Plot style",
                     items       = (('vdw', "Van der Waals", "Plot using VDW radii"),
                              ('bs', "Ball-and-stick", "Ball and stick with bond radii determined automaticaly"),
                              ('fixedbs', "Fixed ball-and-stick", "Ball and stick with a fixed bond radius"),
                              ('sticks', "Stick model", "Just sticks")),
                     default     = 'fixedbs'
                    )

    bond_thickness = FloatProperty(
                  name        = "Thickness of bonds",
                  description = "Determine overall thickness of bonds (0.0 turns bonds off)",
                  default     = 0.15,
                  min         = 0.0,
                  max         = 50.0,
                  step        = 0.05,
                  precision   = 4
                 )

    keystride = IntProperty(
                  name        = "Keystride",
                  description = "Striding between keyframes in animation",
                  default     = 2
                )

    animate_bonds = EnumProperty(
                        name = "Animate bonds",
                        description = "Determine how the bonds are handled for an animation",
                        items = (   ('staticfirst', "Static: first frame", "Use bonds determined from only the first frame"),
                                    ('staticall', "Static: all frames", "Draw all bonds that are formed during any frame"),
                                    ('dynamic', "Dynamically draw frames", "Dynamically form and break bonds during animation")
                                ),
                        default = 'staticfirst'
                    )

    universal_bonds = BoolProperty(
                        name = "Universal bonds",
                        description = "Use one bond type for whole plot",
                        default = True)

    find_aromatic = BoolProperty(
                  name        ="Plot Aromatics",
                  description ="Find closed rings and if planar, fill in with object",
                  default     = False)

    plot_gradient = BoolProperty(
                    name = "Plot gradient",
                    description = "Draw arrows for data found in gradients column",
                    default = False)

    recycle_materials = BoolProperty(
                        name = "Recycle materials",
                        description = "Re-use materials generated for previously imported molecules",
                        default = False)

    color_palette = EnumProperty(
                        name = "Colors",
                        description = "Palette of colors to use",
                        items = ( ('default', 'default', 'Colors defined by molecular blender'),
                                  ('vmd', 'vmd', 'Same colors defined by VMD Element')
                                ),
                        default = 'default')

    def execute(self, context):
        BlendMolecule(context, self.filepath,
                      bonds          = (self.bond_thickness!=0.0),
                      bond_thickness = self.bond_thickness,
                      plot_style     = self.style_of_plot,
                      plot_type      = self.type_of_plot,
                      object_type    = self.type_of_object,
                      keystride      = self.keystride,
                      animate_bonds  = self.animate_bonds,
                      universal_bonds  = self.universal_bonds,
                      plot_gradient  = self.plot_gradient,
                      find_aromatic  = self.find_aromatic,
                      recycle_materials = self.recycle_materials,
                      colors = self.color_palette)

        return {'FINISHED'}

def menu_func_import(self, context):
    self.layout.operator(MolecularBlender.bl_idname, text="Molecular Blender (Import XYZ)")

def register():
    bpy.utils.register_class(MolecularBlender)
    bpy.types.INFO_MT_file_import.append(menu_func_import)

def unregister():
    bpy.utils.unregister_class(MolecularBlender)
    bpy.types.INFO_MT_file_import.remove(menu_func_import)
