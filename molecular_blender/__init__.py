# -*- coding: utf-8 -*-
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

"""Blender plug-in to read files commonly used in Quantum Chemistry"""

import bpy
from bpy.props import StringProperty, BoolProperty, EnumProperty, IntProperty, FloatProperty
import bpy_extras
from .plotter import BlendMolecule
from .marching_cube import CYTHON_ENABLED

bl_info = {  # pylint: disable=invalid-name
    "name": "Molecular Blender",
    "author": "Shane Parker and Josh Szekely",
    "version": (0, 4, 0),
    "location": "File > Import",
    "blender" : (2, 80, 0),
    "description": "Import XYZ files",
    "category": "Import-Export"
}

class MolecularBlender(bpy.types.Operator, bpy_extras.io_utils.ImportHelper):
    """A class designed to conveniently import molecules to blender"""
    bl_idname = "import.molecular_blender"
    bl_label = "Molecular Blender"
    bl_category = "Import"

    # in the template, we'll figure it out later
    filter_glob: StringProperty(
        default="*.xyz;*.molden;*.cube;*.cub;*.json",
        options={'HIDDEN'})

    ## Style Options
    object_type: EnumProperty(
        name="Object type",
        description="Object type to use to draw atoms",
        items=(('meta', "Metaballs", "Metaballs"),
               ('nurbs', "NURBS", "NURBS"),
               ('mesh', "Mesh", "Mesh"),
               ('wireframe', "Wireframe", "Wireframe")),
        default='mesh')

    plot_style: EnumProperty(
        name="Style",
        description="Plot style",
        items=(('vdw', "Van der Waals", "Plot using VDW radii"),
               ('bs', "Ball-and-stick",
                "Ball and stick with bond radii determined automaticaly"),
               ('fixedbs', "Fixed ball-and-stick",
                "Ball and stick with a fixed bond radius"),
               ('sticks', "Stick model", "Just sticks")),
        default='fixedbs')

    bond_thickness: FloatProperty(
        name="Thickness of bonds",
        description="Determine overall thickness of bonds (0.0 turns bonds off)",
        default=0.15,
        min=0.0,
        max=50.0,
        step=0.05,
        precision=4)

    find_aromatic: BoolProperty(
        name="Plot Aromatics",
        description="Find closed rings and if planar, fill in with object",
        default=False)

    ignore_hydrogen: BoolProperty(
        name="Ignore Hydrogens",
        description="Ignores Hydrogen atoms for cleaner images",
        default=False)

    colors: EnumProperty(
        name="Colors",
        description="Palette of colors to use",
        items=(('default', 'default', 'Colors defined by molecular blender'),
               ('vmd', 'vmd', 'Same colors defined by VMD Element')
              ),
        default='default')

    ## Animation
    animate_bonds: EnumProperty(
        name="Bonds",
        description="Determine how the bonds are handled for an animation",
        items=(('staticfirst', "Static: first frame", "Use bonds from only the first frame"),
               ('staticall', "Static: all frames",
                "Draw all bonds that are formed during any frame"),
               ('dynamic', "Dynamically draw frames",
                "Dynamically form and break bonds during animation")
              ),
        default='staticfirst')

    keystride: IntProperty(
        name="Keystride",
        description="Striding between keyframes in animation",
        default=1)

    ## Isosurfaces
    isovalues: StringProperty(
        name="Isovalues to plot",
        description="List of isovalues to plot densities or orbitals",
        default="0.25,0.50")

    cumulative: BoolProperty(
        name="Cumulative isovalues",
        description="Set isovalues so they contain given proportion of the total density",
        default=True)

    volume: EnumProperty(
        name="Volume",
        description="Type of volume to plot",
        items=(('density', 'density', "Density"),
               ('orbital', 'orbital', "Plot orbitals (positive and negative isovalues)")),
        default='orbital')

    orbital: StringProperty(
        name="Orbital",
        description="Plot this orbital from molden file. Accepts comma separated list, homo-n, lumo+m. (ignored for other inputs)",
        default="homo")

    resolution: FloatProperty(
        name="Resolution",
        description="Desired spacing between points in isosurface (Angstrom)",
        default=0.05,
        min=0.001,
        max=0.5,
        step=0.01,
        precision=3)

    remesh: IntProperty(
        name="Remesh Octree Depth",
        description="Octree depth for remesh modifier",
        default=6)

    ## Properties
    gradient: BoolProperty(
        name="Plot gradient",
        description="Draw arrows for data found in gradients column",
        default=False)

    charges: EnumProperty(
        name="Plot charges",
        description="Style in which to plot atomic charges",
        items=(('none', 'none', 'No charges plotted'),
               ('scale', 'scale',
                'Charge magnitude encoded by scale of sphere')
              ),
        default='none')

    charge_offset: FloatProperty(
        name="chgoff",
        description="Use chgfac*(charge + chgoff) to control visibility of charges",
        default=0.9,
        min=-5.0,
        max=5.0,
        step=0.05,
        precision=4
    )

    charge_factor: FloatProperty(
        name="chgfac",
        description="Use chgfac*(charge + chgoff) to control visibility of charges",
        default=1.0,
        min=-100.0,
        max=100.0,
        step=0.05,
        precision=2
    )

    ## Import details
    hook_atoms: EnumProperty(
        name="Hook bonds",
        description="Hook bonds to atoms to make bonds follow manual manipulation of the atoms",
        items=(('auto', "Automatic", "Off for single frames, on for animations"),
               ('on', "On", "Hook bonds (slow for molecules with hundreds or more atoms)"),
               ('off', "Off", "Do not hook bonds to atoms (faster, but manipulating atoms in" +
                " blender will not move the bonds alongside)")),
        default='auto')

    universal_bonds: BoolProperty(
        name="Universal bonds",
        description="Use one bond type for whole plot",
        default=True)

    recycle_materials: BoolProperty(
        name="Recycle materials",
        description="Re-use materials generated for previously imported molecules",
        default=True)

    use_cython: BoolProperty(
        name="Cython Acceleration",
        description="Accelerate isosurface generation and orbital computation using cython",
        default=CYTHON_ENABLED)

    def draw(self, context):
        layout = self.layout

        box = layout.box()
        box.label(text="Style Options")
        box.row().prop(self, "object_type", expand=True)
        box.row().prop(self, "plot_style", expand=True)
        box.row().prop(self, "bond_thickness")
        box.row().prop(self, "find_aromatic")
        box.row().prop(self, "ignore_hydrogen")
        box.row().prop(self, "colors")

        box = layout.box()
        box.label(text="Animation")
        box.row().prop(self, "animate_bonds")
        box.row().prop(self, "keystride")

        box = layout.box()
        box.label(text="Isosurfacing")
        box.row().prop(self, "isovalues")
        box.row().prop(self, "cumulative")
        box.row().prop(self, "volume")
        box.row().prop(self, "orbital")
        box.row().prop(self, "resolution")
        box.row().prop(self, "remesh")

        box = layout.box()
        box.label(text="Properties")
        box.row().prop(self, "gradient")
        box.row().prop(self, "charges")
        box.row().prop(self, "charge_offset")
        box.row().prop(self, "charge_factor")

        box = layout.box()
        box.label(text="Import Details")
        box.row().prop(self, "hook_atoms")
        box.row().prop(self, "universal_bonds")
        box.row().prop(self, "recycle_materials")
        row = box.row()
        row.prop(self, "use_cython")
        row.enabled = False

    def execute(self, context):
        """Function called to plot molecule"""
        BlendMolecule(context, self.filepath,
                      bonds=(self.bond_thickness != 0.0),
                      bond_thickness=self.bond_thickness,
                      hook_atoms=self.hook_atoms,
                      plot_style=self.plot_style,
                      object_type=self.object_type,
                      keystride=self.keystride,
                      animate_bonds=self.animate_bonds,
                      universal_bonds=self.universal_bonds,
                      ignore_hydrogen=self.ignore_hydrogen,
                      gradient=self.gradient,
                      charges=self.charges,
                      charge_offset=self.charge_offset,
                      charge_factor=self.charge_factor,
                      find_aromatic=self.find_aromatic,
                      recycle_materials=self.recycle_materials,
                      isovalues=self.isovalues,
                      cumulative=self.cumulative,
                      volume=self.volume,
                      orbital=self.orbital,
                      resolution=self.resolution,
                      remesh=self.remesh,
                      colors=self.colors)

        return {'FINISHED'}


def menu_func_import(self, _context):
    """Define layout of menu"""
    self.layout.operator(MolecularBlender.bl_idname,
                         text="Molecular Blender (Import XYZ)")

def register():
    from bpy.utils import register_class
    register_class(MolecularBlender)
    bpy.types.TOPBAR_MT_file_import.append(menu_func_import)

def unregister():
    from bpy.utils import unregister_class
    unregister_class(MolecularBlender)
    bpy.types.TOPBAR_MT_file_import.remove(menu_func_import)

if __name__ == "__main__":
    register()
