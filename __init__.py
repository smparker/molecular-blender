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
    "category" : "Import-Export"
}
import bpy

from bpy_extras.io_utils import ImportHelper
from bpy.props import StringProperty, BoolProperty, EnumProperty, IntProperty
from bpy.types import Operator

from . import molecular_blender

class MolecularBlender(Operator, ImportHelper):
    """A class designed to conveniently import molecules to blender"""
    bl_idname = "import.molecular_blender"
    bl_label = "Molecular Blender"

    # in the template, we'll figure it out later
    filename_ext = ".xyz"
    filter_glob = StringProperty(
            default="*.xyz",
            options={'HIDDEN'},
          )

    type_of_object = EnumProperty(
                         name = "Object type",
                         description = "Object type to use to draw atoms",
                         items = (('meta', "Metaballs", "Metaballs"),
                                  ('nurbs', "NURBS", "NURBS"),
                                  ('mesh', "Mesh", "Mesh"),
                                  ('wireframe', "Wireframe", "Wireframe")),
                         default = 'mesh'
                        )

    type_of_plot = EnumProperty(
                    name = "Type",
                    description = "What to draw",
                    items = (('frame', "Single Frame", "Plot a single frame from an XYZ"),
                             ('animate', "Animation", "Animate from an XYZ")),
                    default = 'frame',
                   )

    draw_bonds = BoolProperty(
                  name = "Draw bonds",
                  description = "Toggles whether bonds are drawn",
                  default = True,
                 )

    keystride = IntProperty(
                  name = "Keystride",
                  description = "Striding between keyframes in animation",
                  default = 2,
                )

    def execute(self, context):
        molecular_blender.BlendMolecule(context, self.filepath,
                                              bonds=self.draw_bonds,
                                              plot_type=self.type_of_plot,
                                              object_type=self.type_of_object,
                                              keystride=self.keystride)

        return {'FINISHED'}

def menu_func_import(self, context):
    self.layout.operator(MolecularBlender.bl_idname, text="Molecular Blender (Import XYZ)")

def register():
    bpy.utils.register_class(MolecularBlender)
    bpy.types.INFO_MT_file_import.append(menu_func_import)

def unregister():
    bpy.utils.unregister_class(MolecularBlender)
    bpy.types.INFO_MT_file_import.remove(menu_func_import)
