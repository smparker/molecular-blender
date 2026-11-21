"""Compute and apply a Molden potential to the selected density surface."""
import bpy
from bpy.props import EnumProperty, FloatProperty, IntProperty, StringProperty
import numpy as np

from .diagnostics import report_exception
from .esp_jobs import ESPRequest
from .field_palettes import KEEP_PALETTE, PALETTES
from .importers import file_format
from .scene_esp import BASE_MESH, _display_mesh, _map_values, density_group, field_map_settings, mesh_points, validate_limits
from .scene_state import read_record, resolve_owner
from .scene_surfaces import surface_groups
from .surface_jobs import SurfaceJob
from .wavefunction_esp import HARTREE_IN_VOLTS


class MB_OT_compute_esp(bpy.types.Operator):
    bl_idname = 'molecular_blender.compute_esp'
    bl_label = 'Compute Wavefunction ESP'
    bl_description = 'Color this density surface using the total all-electron Molden wavefunction'
    bl_options = {'UNDO'}

    group: StringProperty(options={'HIDDEN', 'SKIP_SAVE'})
    charge: IntProperty(name='Net molecular charge', description='Charge of the entire all-electron molecule')
    quality: EnumProperty(name='Sampling', items=(
        ('final', 'Final', 'Evaluate every distinct displayed vertex'),
        ('preview', 'Coarse Preview', 'Evaluate one surface vertex per 0.45 Å cell')))
    def _change_units(self, context):
        previous = getattr(self, '_range_units', 'Hartree/e')
        if self.units != previous:
            factor = HARTREE_IN_VOLTS if self.units == 'V' else 1/HARTREE_IN_VOLTS
            self.minimum *= factor
            self.maximum *= factor
        self._range_units = self.units

    units: EnumProperty(name='Potential units', update=_change_units, items=(
        ('Hartree/e', 'Hartree/e', 'Atomic units'), ('V', 'Volts', 'Convert the calculated potential to volts')))
    palette: EnumProperty(name='Color map', items=(KEEP_PALETTE, *PALETTES))
    minimum: FloatProperty(name='Minimum', default=-.05, precision=4)
    maximum: FloatProperty(name='Maximum', default=.05, precision=4)

    def _selection(self, context):
        owner = getattr(self, '_owner', None) or resolve_owner(context.object)
        if owner is None:
            raise ValueError('Select a Molden density surface')
        record = read_record(owner)
        group = surface_groups(owner, record).get(self.group)
        if (file_format(record['source']) != 'molden' or not record['orbitals']
                or not group or not density_group(group)):
            raise ValueError('Select a single-frame Molden density surface')
        return record, group

    def invoke(self, context, event):
        try:
            self._owner = resolve_owner(context.object)
            record, group = self._selection(context)
            mat = group['objects'][0].active_material
            self.charge = mat.get('mb_esp_charge', record['settings']['molecular_charge'])
            settings = field_map_settings(mat)
            if settings and settings['kind'] == 'esp':
                self.units = settings['units']
                self.minimum, self.maximum = settings['limits']
            return context.window_manager.invoke_props_dialog(self, width=410)
        except (KeyError, ValueError, ReferenceError) as exc:
            report_exception(self, exc)
            return {'CANCELLED'}

    def draw(self, context):
        self.layout.label(text='Requires a complete all-electron occupied space.')
        for key in ('charge', 'quality', 'units', 'palette'):
            self.layout.prop(self, key)
        row = self.layout.row(align=True)
        row.prop(self, 'minimum')
        row.prop(self, 'maximum')

    def execute(self, context):
        self._job = self._timer = None
        try:
            self._owner = getattr(self, '_owner', None) or resolve_owner(context.object)
            record, group = self._selection(context)
            if not self.properties.is_property_set('charge'):
                self.charge = record['settings']['molecular_charge']
            validate_limits('esp', self.minimum, self.maximum)
            self._members = [(obj, obj.data, obj.active_material) for obj in group['objects']]
            self._positions = []
            for obj, _, _ in self._members:
                mesh = _display_mesh(context, obj)
                try:
                    self._positions.append(mesh_points(mesh))
                finally:
                    bpy.data.meshes.remove(mesh)
            points = np.concatenate(self._positions)
            if not len(points):
                raise ValueError('The selected density surface has no vertices')
            convention = (group.get('sampling') or {}).get('molden_convention',
                            record.get('molden_convention') or record['settings']['molden_convention'])
            self._request = ESPRequest.from_record(record, points, convention, self.charge, self.quality, self.units)
            if bpy.app.background:
                self._apply(context, self._request.calculate({}))
                return {'FINISHED'}
            self._job = SurfaceJob(bpy.app.binary_path, self._request, {})
            wm = context.window_manager
            self._timer = wm.event_timer_add(.1, window=context.window)
            wm.modal_handler_add(self)
            wm.progress_begin(0, 100)
            context.workspace.status_text_set('Computing wavefunction ESP…  Esc to cancel')
            return {'RUNNING_MODAL'}
        except Exception as exc:
            self.cancel(context)
            report_exception(self, exc)
            return {'CANCELLED'}

    def _apply(self, context, values):
        record, group = self._selection(context)
        current = [(obj, obj.data, obj.active_material) for obj in group['objects']]
        if current != self._members or self._request.source.path != record['source']:
            raise ValueError('The density surface changed during the ESP calculation. Compute it again.')
        offset = 0
        by_object = {}
        for (obj, _, _), points in zip(current, self._positions):
            by_object[obj] = (points, values[offset:offset+len(points)])
            offset += len(points)

        def sample(obj, mesh):
            points, batch = by_object[obj]
            if not np.array_equal(points, mesh_points(mesh)):
                raise ValueError('Displayed vertices changed during the ESP calculation. Compute it again.')
            return batch

        mat = current[0][2]
        keep = self.palette == 'current' and BASE_MESH in current[0][0] and mat.get('mb_field_kind', 'esp') != 'unsigned'
        quality = self.quality
        metadata = dict(mb_wavefunction_esp=True, mb_esp_charge=self.charge,
                        mb_esp_quality=quality, mb_esp_convention=self._request.convention)
        label = 'Wavefunction ESP' + (' (preview)' if quality == 'preview' else '')
        _map_values(context, group, sample, record['source'], 0, self.units, self.minimum, self.maximum,
                    'esp', label, 'default' if self.palette == 'current' else self.palette,
                    metadata, mat if keep else None)

    def modal(self, context, event):
        if event.type == 'ESC' or self._job.cancelled:
            self.cancel(context)
            return {'CANCELLED'}
        if event.type != 'TIMER':
            return {'PASS_THROUGH'}
        try:
            if self._job.running:
                progress = self._job.progress()
                context.window_manager.progress_update(100*progress.get('completed', 0)/max(1, progress.get('total', 1)))
                return {'RUNNING_MODAL'}
            self._apply(context, self._job.result())
        except Exception as exc:
            self.cancel(context)
            report_exception(self, exc)
            return {'CANCELLED'}
        self.cancel(context)
        return {'FINISHED'}

    def cancel(self, context):
        if self._job is not None:
            self._job.close()
            self._job = None
        if self._timer is not None:
            context.window_manager.event_timer_remove(self._timer)
            self._timer = None
            context.window_manager.progress_end()
            context.workspace.status_text_set(None)
