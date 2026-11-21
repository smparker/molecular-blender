"""Render one transparent PNG per saved surface using the scene's camera."""
from contextlib import contextmanager
from pathlib import Path

import bpy
from bpy.props import BoolProperty, CollectionProperty, StringProperty
from mathutils import Quaternion, Vector

from .field_legend import write_field_legend
from .field_scale import scale_values
from .scene_esp import field_map_settings
from .scene_state import owned_objects, read_record, resolve_owner
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


def _view_rotation(context):
    view = context.region_data or getattr(context.space_data, 'region_3d', None)
    if view:
        return (view.view_matrix.to_quaternion().inverted() if view.view_perspective == 'CAMERA'
                else view.view_rotation.copy())
    camera = context.scene.camera
    return camera.matrix_world.to_quaternion() if camera else Quaternion()


@contextmanager
def series_camera(context, owner, groups, keys, rotation):
    """Fit one orthographic render frame around the selected groups in the current pose."""
    if rotation is None:
        yield
        return
    if not keys or any(key not in groups for key in keys):
        raise ValueError('Choose at least one available surface group')
    surfaces = {obj for key in keys for obj in groups[key]['objects']}
    states = {obj: visibility_state(obj) for obj in surfaces}
    scene = context.scene
    low, high = [float('inf')] * 3, [float('-inf')] * 3
    try:
        set_surface_visibility(surfaces, True)
        context.view_layer.update()
        graph = context.evaluated_depsgraph_get()
        to_view = rotation.inverted().to_matrix().to_4x4()
        def rendered_collections(layer):
            if not layer.exclude and not layer.collection.hide_render:
                yield layer.collection
                for child in layer.children:
                    yield from rendered_collections(child)
        collections = set(rendered_collections(context.view_layer.layer_collection))
        for obj in [owner, *owned_objects(owner)]:
            if obj.type not in {'MESH', 'CURVE', 'SURFACE', 'META', 'FONT'} or obj.hide_render:
                continue
            if obj.get('mb_role') == 'surface' and obj not in surfaces:
                continue
            if not any(coll in collections for coll in obj.users_collection):
                continue
            evaluated = obj.evaluated_get(graph)
            bounds = evaluated.bound_box
            if all(tuple(corner) == (-1., -1., -1.) for corner in bounds):
                continue
            transform = to_view @ evaluated.matrix_world
            for corner in bounds:
                point = transform @ Vector(corner)
                for axis in range(3):
                    low[axis] = min(low[axis], point[axis])
                    high[axis] = max(high[axis], point[axis])
    finally:
        for obj, state in states.items():
            restore_visibility(obj, state)
        context.view_layer.update()
    if low[0] > high[0]:
        raise ValueError('The selected molecule and surfaces have no geometry on this frame')

    # Keep a reusable figure camera so the original camera and any animation survive.
    old_camera = scene.camera
    camera = scene.get('mb_series_camera')
    created = not camera or camera.name not in scene.objects or camera.type != 'CAMERA'
    if created:
        camera = bpy.data.objects.new('Molecular surface series', bpy.data.cameras.new('Surface series'))
        scene.collection.objects.link(camera)
        scene['mb_series_camera'] = camera
    data = camera.data
    matrix = camera.matrix_world.copy()
    mode = camera.rotation_mode
    attributes = ('type', 'shift_x', 'shift_y', 'ortho_scale', 'clip_end')
    previous = [getattr(data, name) for name in attributes]
    try:
        data.type = 'ORTHO'
        data.shift_x = data.shift_y = 0
        # view_frame accounts for output dimensions, pixel aspect and sensor fit.
        data.ortho_scale = 1.
        corners = data.view_frame(scene=scene)
        width = max(p.x for p in corners) - min(p.x for p in corners)
        height = max(p.y for p in corners) - min(p.y for p in corners)
        size = Vector(high) - Vector(low)
        render = scene.render
        left, right, bottom, top = ((render.border_min_x, render.border_max_x, render.border_min_y, render.border_max_y)
                                    if render.use_border else (0., 1., 0., 1.))
        if right <= left or top <= bottom:
            raise ValueError('The render border must enclose a nonempty area')
        data.ortho_scale = max(1.1 * max(size.x / (width * (right - left)), size.y / (height * (top - bottom))), .001)
        center = (Vector(high) + Vector(low)) * .5
        center.x -= (right + left - 1.) * .5 * data.ortho_scale * width
        center.y -= (top + bottom - 1.) * .5 * data.ortho_scale * height
        clearance = max(1., data.clip_start * 2, data.ortho_scale)
        center.z = high[2] + clearance
        data.clip_end = max(data.clip_end, clearance + size.z + 1.)
        camera.rotation_mode = 'QUATERNION'
        pose = rotation.to_matrix().to_4x4()
        pose.translation = rotation @ center
        camera.matrix_world = pose
        scene.camera = camera
        context.view_layer.update()
        yield
    except BaseException:
        scene.camera = old_camera
        if created:
            del scene['mb_series_camera']
            bpy.data.objects.remove(camera, do_unlink=True)
            bpy.data.cameras.remove(data)
        else:
            for name, value in zip(attributes, previous):
                setattr(data, name, value)
            camera.rotation_mode = mode
            camera.matrix_world = matrix
        raise


def export_surface_pngs(context, owner, groups, keys, directory, *, legends=False, fit_rotation=None):
    scene = context.scene
    if scene.camera is None and fit_rotation is None:
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
        with png_render_settings(scene.render), series_camera(context, owner, groups, keys, fit_rotation):
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
    fit_camera: BoolProperty(name='Fit series camera from view',
        description='Fit an orthographic render camera around the molecule and all checked surfaces')

    def invoke(self, context, event):
        self._owner = resolve_owner(context.object)
        if self._owner is None:
            return {'CANCELLED'}
        self._view_rotation = _view_rotation(context)
        if context.scene.camera is None:
            self.fit_camera = True
        self.surfaces.clear()
        for key, group in surface_groups(self._owner, read_record(self._owner)).items():
            item = self.surfaces.add()
            item.name = key
            levels = ', '.join(format(level, '.4g') for level in sorted(group['levels']))
            item.label = '{} ({})'.format(group['label'], levels)
        context.window_manager.fileselect_add(self)
        return {'RUNNING_MODAL'}

    def draw(self, context):
        self.layout.prop(self, 'fit_camera')
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
            rotation = getattr(self, '_view_rotation', None) or _view_rotation(context)
            export_surface_pngs(context, owner, groups, keys, directory, legends=self.legends,
                                fit_rotation=rotation if self.fit_camera else None)
        except (OSError, ValueError, RuntimeError, ReferenceError) as exc:
            self.report({'ERROR'}, str(exc))
            return {'CANCELLED'}
        self.report({'INFO'}, 'Exported {} PNG(s) to {}'.format(len(keys), directory))
        return {'FINISHED'}
