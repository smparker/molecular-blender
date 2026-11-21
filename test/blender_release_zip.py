"""Install and render the release ZIP with Blender 4.3 in an empty user configuration.

Run in a temporary working directory with BLENDER_USER_SCRIPTS pointing into it:
blender -b --factory-startup --python-exit-code 1 --python /path/to/this/script -- /path/to/addon.zip
"""
import importlib
from importlib.machinery import ExtensionFileLoader
import os
from pathlib import Path
import sys

import bpy
import numpy as np


def main(archive):
    if bpy.app.version[:2] != (4, 3):
        raise RuntimeError('Run the release ZIP check with Blender 4.3')
    scripts = Path(os.environ['BLENDER_USER_SCRIPTS']).resolve()
    installed = Path(bpy.utils.user_resource('SCRIPTS'), 'addons', 'molecular_blender').resolve()
    if installed != scripts / 'addons' / 'molecular_blender' or installed.exists():
        raise RuntimeError('The ZIP test requires an empty, isolated Blender add-ons directory')
    assert bpy.ops.preferences.addon_install(filepath=str(archive), overwrite=False) == {'FINISHED'}
    assert bpy.ops.preferences.addon_enable(module='molecular_blender') == {'FINISHED'}

    import molecular_blender as addon
    from molecular_blender.plotter import BlendMoleculeSetup, draw_surfaces
    from molecular_blender.scene_state import owned_objects

    assert Path(addon.__file__).resolve().parent == installed
    assert addon.AO_CYTHON_ENABLED and addon.CYTHON_ENABLED
    for name in ('compute_ao_cy', 'marching_cube_cy', 'esp_kernel'):
        module = importlib.import_module('molecular_blender.' + name)
        assert isinstance(module.__loader__, ExtensionFileLoader)
        assert Path(module.__file__).resolve().parent == installed

    source = Path(__file__).resolve().parents[1] / 'examples' / 'fake-atom.molden'
    model, styler, options = BlendMoleculeSetup(bpy.context, str(source), orbital='1',
        resolution=.2, remesh=0, isovalues='.05', cumulative=False)
    draw_surfaces(model, styler, bpy.context, options)
    surfaces = [o for o in owned_objects(model.scene_owner, 'surface') if len(o.data.polygons)]
    assert surfaces
    scene = bpy.context.scene
    for obj in scene.objects:
        if obj.type == 'MESH':
            obj.hide_render = obj not in surfaces
    scene.camera.location = (0, 0, 5)
    scene.camera.rotation_euler = (0, 0, 0)
    scene.camera.data.type = 'ORTHO'
    scene.camera.data.ortho_scale = 5
    scene.render.engine = 'CYCLES'
    scene.cycles.device = 'CPU'
    scene.cycles.samples = 4
    scene.cycles.use_denoising = False
    scene.render.film_transparent = True
    scene.render.resolution_x = scene.render.resolution_y = 64
    scene.render.resolution_percentage = 100
    scene.render.image_settings.file_format = 'PNG'
    scene.render.image_settings.color_mode = 'RGBA'
    scene.render.filepath = str(Path.cwd() / 'release-orbital.png')
    bpy.ops.render.render(write_still=True)
    image = bpy.data.images.load(scene.render.filepath)
    alpha = np.array(image.pixels[:]).reshape(64, 64, 4)[:, :, 3]
    assert np.count_nonzero(alpha > .05) > 40 and alpha[0, 0] == 0, 'Empty or unframed orbital render'
    print(f'Verified installed ZIP in Blender {bpy.app.version_string}: {installed}', flush=True)


if __name__ == '__main__':
    main(Path(sys.argv[sys.argv.index('--') + 1]).resolve())
