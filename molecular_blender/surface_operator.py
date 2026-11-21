"""Shared Molden controls and modal surface computation for Blender operators."""
import bpy
from bpy.props import BoolProperty, EnumProperty, FloatProperty, IntProperty, StringProperty

from . import import_dialog
from .plotter import SurfaceMeshStage, create_surface_meshes
from .surface_jobs import SurfaceJob


DENSITIES = (('density', 'Total density', 'Total electron density'),
             ('alpha', 'Alpha density', 'Alpha electron density'),
             ('beta', 'Beta density', 'Beta electron density'),
             ('spin', 'Spin density', 'Alpha minus beta density'))
ORBITAL_SPINS = (('all', 'All', 'All orbitals'), ('alpha', 'Alpha', 'Alpha frontier orbitals'),
                 ('beta', 'Beta', 'Beta frontier orbitals'))
MOLDEN_FIELDS = DENSITIES + (('orbital', 'Molecular orbital', 'Plot positive and negative orbital lobes'),)


class MoldenSurfaceProperties:
    cube_both_signs: BoolProperty(name="Include both signs", default=True,
        description="Generate both signs for each cube level. Turn off to use exactly the signs entered")

    isovalues: StringProperty(
        name="Isovalues to plot",
        description="Comma-separated surface levels. Leave empty to import just the structure",
        default="0.25,0.50")

    cumulative: BoolProperty(
        name="Cumulative isovalues",
        description="Set isovalues so they contain given proportion of the total density",
        default=True)

    molden_convention: EnumProperty(
        name="Molden convention", default='auto',
        items=(('auto', 'Auto', 'Use ORCA for an orca_2mkl producer title, otherwise standard Molden'),
               ('standard', 'Standard Molden', 'Use standard Molden normalization and phases'),
               ('orca', 'ORCA', 'Use orca_2mkl normalization and phases')))

    # Existing Blender scripts and presets can still explicitly select either convention.
    orca_fix: BoolProperty(options={'HIDDEN'},
        get=lambda op: op.molden_convention == 'orca',
        set=lambda op, value: setattr(op, 'molden_convention', 'orca' if value else 'standard'))

    volume: EnumProperty(
        name="Field",
        description="Field to calculate from the Molden wavefunction",
        items=MOLDEN_FIELDS,
        default='orbital')

    orbital_spin: EnumProperty(
        name="Orbital Spin", description="Channel for frontier-orbital names; explicit numbers always use file indices",
        items=ORBITAL_SPINS, default='all')

    align_orbital_phase: BoolProperty(name="Align orbital phase across frames", default=False,
        description="Use the overlap with the preceding displayed orbital to keep lobe colors consistent")

    orbital: StringProperty(
        name="Orbital",
        description="Plot this orbital from molden file. Accepts comma separated list, homo-n, lumo+m. (ignored for other inputs)",
        default="homo")

    # These UI selectors edit the existing RNA properties used by scripts and presets.
    orbital_choice: EnumProperty(
        name="Orbital", items=(('homo', 'HOMO', '', 0), ('lumo', 'LUMO', '', 1),
                              ('custom', 'Custom…', 'Enter indices, frontier expressions, or a list', 2)),
        get=import_dialog.get_orbital_choice, set=import_dialog.set_orbital_choice,
        options={'SKIP_SAVE'})

    level_mode: EnumProperty(
        name="Level mode", items=(('fraction', 'Fraction enclosed', '', 1),
                                 ('absolute', 'Absolute isovalue', '', 0)),
        get=lambda op: int(op.cumulative), set=lambda op, value: setattr(op, 'cumulative', bool(value)),
        options={'SKIP_SAVE'})

    surface_advanced: BoolProperty(name="Advanced", default=False, options={'SKIP_SAVE'})

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


class SurfaceJobOperator:
    def start_surface_job(self, context, molecule, styler, options, *, select_surfaces):
        self._molecule = molecule
        self._styler = styler
        self._options = options
        self._timer = None
        self._job = None
        self._cancelled = False
        self._select_surfaces = select_surfaces
        self._owner = molecule.scene_owner
        self._collection = molecule.scene_collection
        try:
            trajectory = bool(molecule.volume_trajectory or molecule.orbitals_trajectory)
            self._stage = SurfaceMeshStage(molecule, styler, options) if trajectory else None
            self._job = SurfaceJob(bpy.app.binary_path, molecule.surface_input, options,
                                   stream_frames=trajectory,
                                   on_close=self._stage.discard if trajectory else None)
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
        if event.type == 'ESC' or self._cancelled or not self._owner_exists() or self._job.cancelled:
            self.cancel(context)
            return {'CANCELLED'}
        if event.type != 'TIMER':
            return {'PASS_THROUGH'}
        try:
            progress = self._job.progress()
            context.window_manager.progress_update(
                int(100 * progress.get('completed', 0) / max(progress.get('total', 1), 1)))
            finished = not self._job.running
            if self._stage is not None:
                batch = self._job.take_frame()
                if batch is not None:
                    self._stage.add(*batch)
                    return {'RUNNING_MODAL'}
            if not finished:
                return {'RUNNING_MODAL'}
            if self._stage is not None:
                self._job.check_finished()
                self._stage.finish(context)
            else:
                vertex_sets, trajectory_map = self._job.result()
                if not self._select_surfaces:
                    vertex_sets = [v for v in vertex_sets if len(v["faces"] if "faces" in v else v["triangles"])]
                    if not vertex_sets:
                        self._cleanup(context)
                        self.report({'INFO'}, 'No surface at the requested levels')
                        return {'CANCELLED'}
                self.apply_surface_result(context, vertex_sets, trajectory_map)
        except Exception as exc:
            self._cleanup(context)
            self.report({'ERROR'}, str(exc))
            return {'CANCELLED'}
        self._cleanup(context)
        return {'FINISHED'}

    def apply_surface_result(self, context, vertex_sets, trajectory_map):
        create_surface_meshes(self._molecule, self._styler, context, self._options,
                              vertex_sets, trajectory_map, select_surfaces=self._select_surfaces)

    def cancel(self, context):
        self._cancelled = True
        if self._job is not None:
            self._job.cancel()
        self._cleanup(context)
