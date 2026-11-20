# -*- coding: utf-8 -*-
#
#  Molecular Blender
#  Filename: __init__.py
#  Copyright (C) 2014-2025 Shane Parker, Joshua Szekely
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

from .marching_cube import CYTHON_ENABLED
from .compute_ao import CYTHON_ENABLED as AO_CYTHON_ENABLED

from . import importers
from . import orbitals
from . import isosurfaces

try:
    import bpy
except ModuleNotFoundError as exc:
    if exc.name != 'bpy':
        raise
    _HAS_BPY = False
else:
    _HAS_BPY = True

if _HAS_BPY:
    from bpy.props import StringProperty, BoolProperty, EnumProperty, IntProperty, FloatProperty
    import bpy_extras
    from . import import_dialog
    from .plotter import BlendMolecule, BlendMoleculeSetup, compute_surface_vertex_sets, create_surface_meshes
    from .surface_operator import MoldenSurfaceProperties, SurfaceJobOperator
    from .surface_jobs import cancel_all_jobs
    from bpy.app.handlers import persistent
    cancel_all_jobs = persistent(cancel_all_jobs)

bl_info = {  # pylint: disable=invalid-name
    "name": "Molecular Blender",
    "author": "Shane Parker and Josh Szekely",
    "version": (0, 6, 3),
    "location": "File > Import",
    "blender" : (3, 6, 7),
    "description": "Import XYZ, SDF/MOL, Molden, cube, and JSON molecular data",
    "category": "Import-Export"
}

if _HAS_BPY:
    class MolecularBlender(MoldenSurfaceProperties, SurfaceJobOperator, bpy.types.Operator, bpy_extras.io_utils.ImportHelper):
        """A class designed to conveniently import molecules to blender"""
        bl_idname = "import.molecular_blender"
        bl_label = "Molecular Blender"
        bl_category = "Import"

        # in the template, we'll figure it out later
        filter_glob: StringProperty(
            default=";".join("*" + suffix for suffixes in importers.FORMAT_SUFFIXES.values() for suffix in suffixes),
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

        draw_bond_orders: BoolProperty(
            name="Draw bond orders",
            description="Draw supplied bond orders, or perceive them from coordinates (requires all hydrogens)",
            default=False)

        molecular_charge: IntProperty(
            name="Molecular charge",
            description="Total integer charge for bond-order perception (not partial atomic charges)",
            default=0)

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

        trajectory_interpolation: EnumProperty(name="Playback", items=(
            ('blender', 'Blender Default', 'Use the keyframe interpolation from Blender preferences'),
            ('hold', 'Hold', 'Show each sampled geometry until the next snapshot'),
            ('linear', 'Linear', 'Interpolate atom positions and charge sizes linearly between snapshots')))

        bond_break_margin: FloatProperty(name="Break margin (%)", default=5.0, min=0, soft_max=25,
            description="Extra distance above the formation cutoff before an inferred bond disappears")

        frame_start: IntProperty(name="First input frame", default=0, min=0,
            description="Zero-based position in the XYZ trajectory or JSON entry list")
        frame_end: IntProperty(name="Last input frame", default=-1, min=-1,
            description="Inclusive input frame, or -1 to read through the end")
        frame_step: IntProperty(name="Read every Nth frame", default=1, min=1,
            description="Sample starting at the first selected input frame and preserve source timing")

        cube_field_type: EnumProperty(
            name="Cube field type",
            description="How cumulative cube surfaces integrate values",
            items=(('auto', 'Auto / scalar magnitude', 'Use file metadata; untyped fields integrate absolute magnitude'),
                   ('orbital', 'Orbital amplitude', 'Integrate squared amplitude'),
                   ('density', 'Density', 'Integrate absolute density'),
                   ('signed', 'Signed scalar', 'Integrate absolute field magnitude')),
            default='auto')

        cube_dataset: IntProperty(
            name="Cube dataset index", description="Zero-based dataset index within a cube file",
            default=0, min=0)

        dataset_choice: EnumProperty(
            name="Dataset", items=import_dialog.dataset_items, get=import_dialog.get_dataset_choice,
            set=lambda op, value: setattr(op, 'cube_dataset', value), options={'SKIP_SAVE'})

        cube_both_signs: BoolProperty(
            name="Include both signs", default=True,
            description="Generate both signs for each cube level. Turn off to use exactly the signs entered")

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
            items=(('auto', "Automatic", "Hook bonds for animated trajectories and Molden vibrational modes"),
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
            if importers.file_format(import_dialog.browser_path(self, context)) in ('xyz', 'json'):
                box.prop(self, "trajectory_interpolation")
                if self.animate_bonds == 'dynamic':
                    box.prop(self, "bond_break_margin")
                row = box.row(align=True)
                row.prop(self, "frame_start")
                row.prop(self, "frame_end")
                box.prop(self, "frame_step")

            import_dialog.draw_surfaces(layout, self, context)

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
            box.row().prop(self, "draw_bond_orders")
            if self.draw_bond_orders:
                box.row().prop(self, "molecular_charge")
            box.row().prop(self, "recycle_materials")
            box.label(text="Orbital evaluator: " + ("Cython" if AO_CYTHON_ENABLED else "Python"))
            box.label(text="Surface mesher: " + ("Cython" if CYTHON_ENABLED else "Python"))

        def _build_options(self, source=None, from_browser=False):
            options = dict(
                bonds=(self.bond_thickness != 0.0),
                bond_thickness=self.bond_thickness,
                draw_bond_orders=self.draw_bond_orders,
                perceive_bond_orders=self.draw_bond_orders,
                molecular_charge=self.molecular_charge,
                hook_atoms=self.hook_atoms,
                plot_style=self.plot_style,
                object_type=self.object_type,
                keystride=self.keystride,
                trajectory_interpolation=self.trajectory_interpolation,
                bond_break_margin=self.bond_break_margin,
                frame_start=self.frame_start, frame_end=self.frame_end, frame_step=self.frame_step,
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
                cube_dataset=self.cube_dataset,
                cube_field_type=self.cube_field_type,
                orbital=self.orbital,
                orbital_spin=self.orbital_spin,
                resolution=self.resolution,
                remesh=self.remesh,
                orca_fix=self.orca_fix,
                colors=self.colors,
            )

            if source is not None:
                kind, labels = source
                if kind in ('xyz', 'sdf'):
                    options['isovalues'] = ''
                elif kind == 'cube':
                    if from_browser:
                        options['cube_dataset'] = import_dialog.dataset_index(self.cube_dataset, labels)
                    if from_browser or self.properties.is_property_set('cube_both_signs'):
                        options['cube_both_signs'] = self.cube_both_signs
            return options

        def execute(self, context):
            """Load molecule, plot atoms/bonds, then compute isosurfaces in background if needed."""
            try:
                source = import_dialog.surface_source(self.filepath, self.frame_start)
                settings = self._build_options(source, self.options.is_invoke)
                molecule, styler, options = BlendMoleculeSetup(context, self.filepath, **settings)
            except Exception as e:
                self.report({'ERROR'}, str(e))
                return {'CANCELLED'}

            if not options.get("isosurfaces"):
                return {'FINISHED'}

            return self.start_surface_job(context, molecule, styler, options, select_surfaces=True)


    def menu_func_import(self, _context):
        """Define layout of menu"""
        self.layout.operator(MolecularBlender.bl_idname,
                             text="Molecular Blender (XYZ, SDF/MOL, Molden, Cube, JSON)")

def register():
    if not _HAS_BPY:
        raise ImportError("register() requires bpy (must be run inside Blender)")
    from bpy.utils import register_class
    register_class(MolecularBlender)
    bpy.types.TOPBAR_MT_file_import.append(menu_func_import)
    bpy.app.handlers.load_pre.append(cancel_all_jobs)
    from . import scene_state, scene_panel, scene_edit, bond_edit
    scene_state.register()
    scene_edit.register()
    bond_edit.register()
    scene_panel.register()

def unregister():
    if not _HAS_BPY:
        raise ImportError("unregister() requires bpy (must be run inside Blender)")
    from bpy.utils import unregister_class
    cancel_all_jobs()
    from . import scene_state, scene_panel, scene_edit, bond_edit
    scene_panel.unregister()
    bond_edit.unregister()
    scene_edit.unregister()
    scene_state.unregister()
    if cancel_all_jobs in bpy.app.handlers.load_pre:
        bpy.app.handlers.load_pre.remove(cancel_all_jobs)
    unregister_class(MolecularBlender)
    bpy.types.TOPBAR_MT_file_import.remove(menu_func_import)

if __name__ == "__main__":
    register()
