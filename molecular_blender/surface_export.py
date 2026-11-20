"""Render one transparent PNG per saved surface using the scene's camera."""
from contextlib import contextmanager
from pathlib import Path

import bpy
from bpy.props import BoolProperty, CollectionProperty, StringProperty

from .field_legend import write_field_legend
from .field_scale import scale_values
from .scene_esp import field_map_settings
from .scene_state import read_record, resolve_owner
from .scene_surfaces import restore_visibility, set_surface_visibility, surface_groups, visibility_state


def set_transparent_png(render):
    render.film_transparent = True
    render.image_settings.file_format = 'PNG'
    render.image_settings.color_mode = 'RGBA'


class MB_OT_transparent_png(bpy.types.Operator):
    bl_idname = 'molecular_blender.transparent_png'
    bl_label = 'Use Transparent PNG'
    bl_description = 'Set the scene output to RGBA PNG with a transparent film background'
    bl_options = {'UNDO'}

    def execute(self, context):
        set_transparent_png(context.scene.render)
        return {'FINISHED'}


def draw_png_preset(self, context):
    self.layout.operator(MB_OT_transparent_png.bl_idname, icon='IMAGE_DATA')


@contextmanager
def png_render_settings(render):
    image = render.image_settings
    previous = (render.filepath, render.film_transparent, render.use_multiview, render.use_sequencer,
                image.file_format, image.color_depth, image.color_mode)
    try:
        set_transparent_png(render)
        render.use_multiview = False
        render.use_sequencer = False
        yield
    finally:
        (render.filepath, render.film_transparent, render.use_multiview, render.use_sequencer,
         image.file_format, image.color_depth, image.color_mode) = previous


def reserve_outputs(directory, stem, legend):
    number = 1
    while True:
        name = stem if number == 1 else f'{stem}-{number}'
        png = directory / (name + '.png')
        svg = directory / (name + '.svg') if legend else None
        try:
            png.touch(exist_ok=False)
            try:
                if svg:
                    svg.touch(exist_ok=False)
            except BaseException:
                png.unlink()
                raise
            return png, svg
        except FileExistsError:
            number += 1


def _legend_material(group):
    materials = {obj.active_material for obj in group['objects']}
    signatures = set()
    for mat in materials:
        settings = field_map_settings(mat)
        signatures.add((settings['label'], settings['units'], settings['kind'] == 'unsigned', scale_values(mat))
                       if settings else None)
    if len(signatures) != 1:
        raise ValueError(f"{group['label']} uses several color scales. Choose a common scale to export one legend.")
    return next(iter(materials)) if None not in signatures else None


def export_surface_pngs(context, owner, groups, keys, directory, *, legends=False):
    scene = context.scene
    if scene.camera is None:
        raise ValueError('Set an active scene camera before exporting surfaces')
    if not directory.is_dir():
        raise ValueError('Choose an existing output folder')
    if not keys:
        raise ValueError('Choose at least one surface to export')
    if any(key not in groups for key in keys):
        raise ValueError('A selected surface was removed. Open the export dialog again.')
    mapped = {key: _legend_material(groups[key]) for key in keys} if legends else {}
    visibility = {obj: visibility_state(obj) for group in groups.values() for obj in group['objects']}
    layer = context.view_layer
    selected = [obj for obj in layer.objects if obj.select_get(view_layer=layer)]
    active = layer.objects.active
    wm = context.window_manager
    wm.progress_begin(0, len(keys))
    try:
        with png_render_settings(scene.render):
            for i, key in enumerate(keys):
                for group_key, group in groups.items():
                    for obj in group['objects']:
                        if obj.animation_data:
                            # Collection gating preserves the current frame's visibility keys.
                            set_surface_visibility([obj], group_key == key)
                        else:
                            obj.hide_render = group_key != key
                stem = bpy.path.clean_name(owner.name + '-' + groups[key]['label'], replace='-')
                path, svg = reserve_outputs(directory, stem, mapped.get(key) is not None)
                try:
                    if svg:
                        write_field_legend(mapped[key], svg)
                    scene.render.filepath = str(path)
                    status = bpy.ops.render.render(write_still=True, scene=scene.name)
                    if 'FINISHED' not in status or path.stat().st_size == 0:
                        raise RuntimeError('Render did not finish for ' + groups[key]['label'])
                except BaseException:
                    path.unlink(missing_ok=True)
                    if svg:
                        svg.unlink(missing_ok=True)
                    raise
                wm.progress_update(i + 1)
    finally:
        for obj, state in visibility.items():
            restore_visibility(obj, state)
        layer.update()
        for obj in selected:
            if not obj.select_get(view_layer=layer):
                obj.select_set(True, view_layer=layer)
        layer.objects.active = active
        wm.progress_end()


class MB_ExportSurface(bpy.types.PropertyGroup):
    include: BoolProperty(name='Export', default=True)
    label: StringProperty()


class MB_OT_export_surface_series(bpy.types.Operator):
    bl_idname = 'molecular_blender.export_surface_series'
    bl_label = 'Export Surface PNGs'
    bl_description = 'Render the selected surface groups as separate transparent PNGs with the active scene camera'

    directory: StringProperty(name='Output folder', subtype='DIR_PATH', options={'SKIP_SAVE'})
    filter_folder: BoolProperty(default=True, options={'HIDDEN'})
    surfaces: CollectionProperty(type=MB_ExportSurface, options={'SKIP_SAVE'})
    legends: BoolProperty(name='Include field legends', description='Write a matching SVG for each mapped panel')

    def invoke(self, context, event):
        self._owner = resolve_owner(context.object)
        if self._owner is None:
            return {'CANCELLED'}
        if context.scene.camera is None:
            self.report({'WARNING'}, 'Set an active scene camera before exporting surfaces')
            return {'CANCELLED'}
        self.surfaces.clear()
        for key, group in surface_groups(self._owner, read_record(self._owner)).items():
            item = self.surfaces.add()
            item.name = key
            levels = ', '.join(format(level, '.4g') for level in sorted(group['levels']))
            item.label = '{} ({})'.format(group['label'], levels)
        context.window_manager.fileselect_add(self)
        return {'RUNNING_MODAL'}

    def draw(self, context):
        self.layout.prop(self, 'legends')
        self.layout.label(text='Surfaces to render')
        for item in self.surfaces:
            self.layout.prop(item, 'include', text=item.label)

    def execute(self, context):
        try:
            owner = self._owner if hasattr(self, '_owner') else resolve_owner(context.object)
            if owner is None or not self.directory:
                raise ValueError('Choose a molecule and an output folder')
            groups = surface_groups(owner, read_record(owner))
            keys = [item.name for item in self.surfaces if item.include] if self.surfaces else list(groups)
            directory = Path(bpy.path.abspath(self.directory)).resolve()
            export_surface_pngs(context, owner, groups, keys, directory, legends=self.legends)
        except (OSError, ValueError, RuntimeError, ReferenceError) as exc:
            self.report({'ERROR'}, str(exc))
            return {'CANCELLED'}
        self.report({'INFO'}, 'Exported {} PNG(s) to {}'.format(len(keys), directory))
        return {'FINISHED'}
