"""Compare resident and packed surface playback in separate fresh Blender processes.

blender -b --factory-startup --python-exit-code 1 --python benchmarks/trajectory_meshes.py -- packed
Repeat with resident to measure the same 24 signed cube snapshots and 48 lobes.
"""
import gc
import json
from pathlib import Path
import resource
from statistics import median
import sys
import tempfile
from time import perf_counter

import bpy
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from molecular_blender.plotter import BlendMoleculeSetup, draw_surfaces
from molecular_blender.scene_state import owned_objects


def write_inputs(folder):
    size = 46
    spacing = .13
    axis = (np.arange(size) - (size - 1) / 2) * spacing
    x, y, z = np.meshgrid(axis, axis, axis, indexing='ij')
    entries = []
    for step in range(24):
        path = folder / f'step-{step}.cube'
        shift = .4 * np.sin(step / 24 * 2 * np.pi)
        field = (x - shift) * np.exp(-((x - shift)**2 + y*y + z*z) / 2)
        with path.open('w') as stream:
            stream.write(f'Animated signed field\nGenerated benchmark\n1 {axis[0]} {axis[0]} {axis[0]}\n'
                         f'-{size} {spacing} 0 0\n-{size} 0 {spacing} 0\n-{size} 0 0 {spacing}\n2 2 0 0 0\n')
            np.savetxt(stream, field.reshape(-1, 8), fmt='%.7g')
        entries.append(dict(filename=path.name, frame=step))
    manifest = folder / 'series.json'
    manifest.write_text(json.dumps(dict(molecules=entries)))
    return manifest


def measure(folder, mode):
    path = write_inputs(folder)
    start = perf_counter()
    model, styler, options = BlendMoleculeSetup(
        bpy.context, str(path), trajectory_meshes=mode, remesh=0, cumulative=False,
        isovalues='.15', cube_both_signs=True)
    draw_surfaces(model, styler, bpy.context, options)
    wall = perf_counter() - start
    objects = owned_objects(model.scene_owner, 'surface')
    del model, styler, options
    gc.collect()
    assert len(objects) == (2 if mode == 'packed' else 48)

    scene = bpy.context.scene
    timings, counts = [], []
    for frame in range(1, 25):
        start = perf_counter()
        scene.frame_set(frame)
        bpy.context.view_layer.update()
        timings.append((perf_counter() - start) * 1000)
        counts.append(sum(len(obj.data.polygons) for obj in objects if not obj.hide_render))
    assert min(counts) > 0
    blend = folder / 'result.blend'
    bpy.ops.wm.save_as_mainfile(filepath=str(blend), compress=False)
    rss_units = 2**20 if sys.platform == 'darwin' else 1024
    peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / rss_units
    print(f'{mode}: {len(objects)} objects, {sum(counts)} faces across playback, '
          f'{wall:.3f}s import, {median(timings):.2f}ms median frame, {peak:.1f} MiB peak RSS, '
          f'{blend.stat().st_size / 2**20:.2f} MiB blend', flush=True)


if __name__ == '__main__':
    with tempfile.TemporaryDirectory() as tmp:
        measure(Path(tmp), sys.argv[sys.argv.index('--') + 1])
