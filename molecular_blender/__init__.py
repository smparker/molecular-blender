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
    from .plotter import BlendMolecule, BlendMoleculeSetup, compute_surface_vertex_sets, create_surface_meshes
    from .surface_jobs import SurfaceJob, cancel_all_jobs
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
    class MolecularBlender(bpy.types.Operator, bpy_extras.io_utils.ImportHelper):
        """A class designed to conveniently import molecules to blender"""
        bl_idname = "import.molecular_blender"
        bl_label = "Molecular Blender"
        bl_category = "Import"

        # in the template, we'll figure it out later
        filter_glob: StringProperty(
            default="*.xyz;*.sdf;*.sd;*.mol;*.molden;*molden.input;*.cube;*.cub;*.json",
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

        ## Isosurfaces
        isovalues: StringProperty(
            name="Isovalues to plot",
            description="List of isovalues to plot densities or orbitals",
            default="0.25,0.50")

        cumulative: BoolProperty(
            name="Cumulative isovalues",
            description="Set isovalues so they contain given proportion of the total density",
            default=True)

        orca_fix: BoolProperty(
                name="Orca molden fix",
                description="Fix for orca molden files with incorrect AO norms",
                default=False)

        volume: EnumProperty(
            name="Volume",
            description="Type of volume to plot",
            items=(('density', 'Total density', "Total electron density"),
                   ('alpha', 'Alpha density', "Alpha electron density"),
                   ('beta', 'Beta density', "Beta electron density"),
                   ('spin', 'Spin density', "Alpha minus beta density"),
                   ('orbital', 'orbital', "Plot orbitals (positive and negative isovalues)")),
            default='orbital')

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

        orbital_spin: EnumProperty(
            name="Orbital spin", description="Channel for frontier-orbital names; explicit numbers always use file indices",
            items=(('all', 'All', 'All orbitals'), ('alpha', 'Alpha', 'Alpha frontier orbitals'),
                   ('beta', 'Beta', 'Beta frontier orbitals')), default='all')

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
            default=8)

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
            if self.volume in {"alpha", "beta", "spin"}:
                box.label(text="Spin populations use labels and occupation-based defaults")
            box.row().prop(self, "cube_dataset")
            box.row().prop(self, "cube_field_type")
            box.row().prop(self, "orbital")
            box.row().prop(self, "orbital_spin")
            box.row().prop(self, "resolution")
            box.row().prop(self, "remesh")
            box.row().prop(self, "orca_fix")

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

        def _build_options(self):
            return dict(
                bonds=(self.bond_thickness != 0.0),
                bond_thickness=self.bond_thickness,
                draw_bond_orders=self.draw_bond_orders,
                perceive_bond_orders=self.draw_bond_orders,
                molecular_charge=self.molecular_charge,
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
                cube_dataset=self.cube_dataset,
                cube_field_type=self.cube_field_type,
                orbital=self.orbital,
                orbital_spin=self.orbital_spin,
                resolution=self.resolution,
                remesh=self.remesh,
                orca_fix=self.orca_fix,
                colors=self.colors,
            )

        def execute(self, context):
            """Load molecule, plot atoms/bonds, then compute isosurfaces in background if needed."""
            try:
                molecule, styler, options = BlendMoleculeSetup(context, self.filepath, **self._build_options())
            except Exception as e:
                self.report({'ERROR'}, str(e))
                return {'CANCELLED'}

            if not options.get("isosurfaces"):
                return {'FINISHED'}

            self._molecule = molecule
            self._styler = styler
            self._options = options
            self._timer = None
            self._job = None
            self._cancelled = False
            self._owner = molecule.scene_owner
            self._collection = molecule.scene_collection
            try:
                self._job = SurfaceJob(bpy.app.binary_path, molecule.surface_input, options)
                wm = context.window_manager
                self._timer = wm.event_timer_add(0.1, window=context.window)
                wm.modal_handler_add(self)
                wm.progress_begin(0, 100)
                context.workspace.status_text_set("Molecular Blender: computing isosurfaces... (Esc to cancel)")
            except Exception as exc:
                self.cancel(context)
                self.report({'ERROR'}, str(exc))
                return {'CANCELLED'}
            return {'RUNNING_MODAL'}

        def _owner_exists(self):
            try:
                return (self._owner is not None and self._collection is not None
                        and bpy.data.objects.get(self._owner.name) == self._owner
                        and bpy.data.collections.get(self._collection.name) == self._collection)
            except ReferenceError:
                return False

        def _cleanup(self, context):
            if self._job is not None:
                self._job.close()
                self._job = None
            if self._timer is not None:
                context.window_manager.event_timer_remove(self._timer)
                self._timer = None
            context.window_manager.progress_end()
            if context.workspace is not None:
                context.workspace.status_text_set(None)

        def modal(self, context, event):
            if event.type == 'ESC' or self._cancelled or not self._owner_exists():
                self.cancel(context)
                return {'CANCELLED'}
            if event.type != 'TIMER':
                return {'PASS_THROUGH'}
            try:
                progress = self._job.progress()
                context.window_manager.progress_update(
                    int(100 * progress.get('completed', 0) / max(progress.get('total', 1), 1)))
                if self._job.running:
                    return {'RUNNING_MODAL'}
                vertex_sets, trajectory_map = self._job.result()
                # A renamed owner remains valid; scene generation uses these handles.
                self._molecule.scene_owner = self._owner
                self._molecule.scene_collection = self._collection
                create_surface_meshes(self._molecule, self._styler, context,
                                      self._options, vertex_sets, trajectory_map)
            except Exception as exc:
                self._cleanup(context)
                self.report({'ERROR'}, str(exc))
                return {'CANCELLED'}
            finally:
                if self._job is not None and not self._job.running:
                    self._cleanup(context)
            return {'FINISHED'}

        def cancel(self, context):
            self._cancelled = True
            if self._job is not None:
                self._job.cancel()
            self._cleanup(context)


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
    from . import scene_state, scene_panel
    scene_state.register()
    scene_panel.register()

def unregister():
    if not _HAS_BPY:
        raise ImportError("unregister() requires bpy (must be run inside Blender)")
    from bpy.utils import unregister_class
    cancel_all_jobs()
    from . import scene_state, scene_panel
    scene_panel.unregister()
    scene_state.unregister()
    if cancel_all_jobs in bpy.app.handlers.load_pre:
        bpy.app.handlers.load_pre.remove(cancel_all_jobs)
    unregister_class(MolecularBlender)
    bpy.types.TOPBAR_MT_file_import.remove(menu_func_import)

if __name__ == "__main__":
    register()
