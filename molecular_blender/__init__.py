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

from .diagnostics import report_exception
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
    from .options import IMPORT_PROPERTIES
    from .ui_options import option_property
    from .plotter import BlendMolecule, BlendMoleculeSetup, compute_surface_vertex_sets, create_surface_meshes
    from .surface_operator import MOLDEN_FIELDS, MoldenSurfaceProperties, SurfaceJobOperator
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
        object_type: option_property('object_type')

        plot_style: option_property('plot_style')

        bond_thickness: option_property('bond_thickness')

        find_aromatic: option_property('find_aromatic')

        draw_bond_orders: option_property('draw_bond_orders')

        molecular_charge: option_property('molecular_charge')

        ignore_hydrogen: option_property('ignore_hydrogen')

        colors: option_property('colors')

        ## Animation
        animate_bonds: option_property('animate_bonds')

        keystride: option_property('keystride')
        trajectory_timing: option_property('trajectory_timing')
        time_step: option_property('time_step')
        time_unit: option_property('time_unit')
        playback_ps_per_second: option_property('playback_ps_per_second')
        playback_fps: option_property('playback_fps')
        trajectory_meshes: option_property('trajectory_meshes')

        trajectory_interpolation: option_property('trajectory_interpolation')

        bond_break_margin: option_property('bond_break_margin')

        frame_start: option_property('frame_start')
        frame_end: option_property('frame_end')
        frame_step: option_property('frame_step')

        def _set_response_field(self, value):
            self['show_response'] = value == len(MOLDEN_FIELDS)
            if not self['show_response']:
                self.volume = MOLDEN_FIELDS[value][0]

        field_with_response: EnumProperty(name='Field',
            items=MOLDEN_FIELDS + (('response', 'Response density', 'Signed Turbomole polarization density'),),
            get=lambda op: len(MOLDEN_FIELDS) if op.get('show_response') else
                           next(i for i, field in enumerate(MOLDEN_FIELDS) if field[0] == op.volume),
            set=_set_response_field)
        response_choice: EnumProperty(name='Eigenpair', items=import_dialog.response_items,
            get=import_dialog.get_response_choice, set=import_dialog.set_response_choice,
            options={'SKIP_SAVE'})

        dataset_choice: EnumProperty(
            name="Dataset", items=import_dialog.dataset_items, get=import_dialog.get_dataset_choice,
            set=lambda op, value: setattr(op, 'cube_dataset', value), options={'SKIP_SAVE'})

        ## Properties
        gradient: option_property('gradient')

        charges: option_property('charges')

        charge_offset: option_property('charge_offset')

        charge_factor: option_property('charge_factor')

        ## Import details
        hook_atoms: option_property('hook_atoms')

        universal_bonds: option_property('universal_bonds')

        recycle_materials: option_property('recycle_materials')

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

            input_format = importers.file_format(import_dialog.browser_path(self, context))
            if input_format in ('xyz', 'json'):
                box.prop(self, 'trajectory_timing')
                if self.trajectory_timing == 'physical':
                    row = box.row(align=True)
                    row.prop(self, 'time_step')
                    row.prop(self, 'time_unit', text='')
                    box.prop(self, 'playback_ps_per_second')
                    box.prop(self, 'playback_fps')
                else:
                    box.prop(self, 'keystride')
                if input_format == 'json':
                    box.prop(self, 'trajectory_meshes')
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
            options = {key: getattr(self, key) for key in IMPORT_PROPERTIES}
            options.update(bonds=self.bond_thickness != 0.0, perceive_bond_orders=self.draw_bond_orders)

            if source is not None:
                kind, labels = source
                if (kind == 'molden' and self.field_with_response == 'response'
                        and import_dialog.response_path(self.filepath, self.frame_start)):
                    options['volume'] = 'response'
                    if from_browser:
                        items = import_dialog.response_items(self, None)
                        if str(self.response_eigenpair) not in {item[0] for item in items}:
                            options['response_eigenpair'] = int(items[0][0])
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
                report_exception(self, e)
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
    from . import scene_state, scene_panel, scene_edit, bond_edit, render_cli, trajectory_meshes
    scene_state.register()
    scene_edit.register()
    bond_edit.register()
    scene_panel.register()
    trajectory_meshes.register()
    render_cli.register()

def unregister():
    if not _HAS_BPY:
        raise ImportError("unregister() requires bpy (must be run inside Blender)")
    from bpy.utils import unregister_class
    cancel_all_jobs()
    from . import scene_state, scene_panel, scene_edit, bond_edit, render_cli, trajectory_meshes
    trajectory_meshes.unregister()
    render_cli.unregister()
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
