"""Measure actual viewport redraws and Escape-to-worker-cleanup time in a fresh Blender window.

blender --factory-startup --enable-event-simulate --python-exit-code 1 --python benchmarks/surface_responsiveness.py
Append `-- --fallback` to enable the Python backends in both the UI and worker.
The benchmark closes its own Blender process after three trials.
"""
import json
import os
from pathlib import Path
import shutil
import sys
import tempfile
import time
import traceback

import addon_utils
import bpy
import numpy as np

if not {'--factory-startup', '--enable-event-simulate'}.issubset(sys.argv):
    raise RuntimeError('Launch a fresh Blender with the flags in this script’s docstring')
repo = Path(__file__).resolve().parents[1]
fallback = '--fallback' in sys.argv
temporary = tempfile.TemporaryDirectory()
package_root = repo
if fallback:
    package_root = Path(temporary.name)
    package = package_root / 'molecular_blender'
    package.mkdir()
    for module in (repo / package.name).glob('*.py'):
        shutil.copy(module, package)
sys.path.insert(0, str(package_root))
addon_utils.enable('molecular_blender')
from molecular_blender import AO_CYTHON_ENABLED, CYTHON_ENABLED, surface_jobs
from molecular_blender.plotter import BlendMoleculeSetup

assert AO_CYTHON_ENABLED == CYTHON_ENABLED == (not fallback)
source = repo / 'examples' / 'phenanthroline.molden'
draws, intervals, cancellations = [], [], []
stage = 'setup'
deadline = time.perf_counter() + 100


def on_draw():
    if stage == 'measuring':
        draws.append(time.perf_counter())


def finish(status):
    surface_jobs.cancel_all_jobs()
    temporary.cleanup()
    sys.stdout.flush()
    sys.stderr.flush()
    os._exit(status)


def tick():
    global stage, job, started, escape
    try:
        now = time.perf_counter()
        if now > deadline:
            raise TimeoutError('Benchmark timed out during ' + stage)
        window = bpy.context.window_manager.windows[0]
        area = next(a for a in window.screen.areas if a.type == 'VIEW_3D')
        region = next(r for r in area.regions if r.type == 'WINDOW')
        with bpy.context.temp_override(window=window, area=area, region=region):
            if stage == 'setup':
                model, _, _ = BlendMoleculeSetup(bpy.context, str(source), isovalues='')
                bpy.context.view_layer.objects.active = model.scene_owner
                bpy.types.SpaceView3D.draw_handler_add(on_draw, (), 'WINDOW', 'POST_PIXEL')
                stage = 'start'
            if stage == 'start':
                result = bpy.ops.molecular_blender.add_surface(orbital='homo-2,homo-1,homo,lumo,lumo+1,lumo+2',
                    resolution=.025, cumulative=False, isovalues='.05', remesh=0)
                if result != {'RUNNING_MODAL'}:
                    raise RuntimeError(str(result))
                job, = surface_jobs._jobs
                stage = 'wait_for_sampling'
            elif stage in ('wait_for_sampling', 'measuring'):
                if not job.running:
                    raise RuntimeError('Worker exited before the sampling window finished')
                if stage == 'wait_for_sampling' and (job.directory / 'progress.json').exists():
                    draws.clear()
                    stage, started = 'measuring', now
                elif stage == 'measuring' and now - started >= 1.2:
                    if len(draws) < 2:
                        raise RuntimeError('The viewport did not redraw twice')
                    intervals.extend(np.diff(draws) * 1000)
                    escape = time.perf_counter()
                    window.event_simulate(type='ESC', value='PRESS', x=area.x+area.width//2, y=area.y+area.height//2)
                    stage = 'wait_for_cleanup'
            elif stage == 'wait_for_cleanup' and job not in surface_jobs._jobs:
                assert job.cancelled and not job.running and not job.directory.exists()
                cancellations.append((now - escape) * 1000)
                if len(cancellations) == 3:
                    print(json.dumps(dict(blender=bpy.app.version_string, backend='Python' if fallback else 'Cython',
                        redraw_intervals=len(intervals), median_redraw_ms=float(np.median(intervals)),
                        p95_redraw_ms=float(np.percentile(intervals, 95)), cancellation_ms=cancellations)))
                    finish(0)
                stage = 'start'
            area.tag_redraw()
        return .01
    except Exception:
        traceback.print_exc()
        finish(1)


bpy.app.timers.register(tick, first_interval=.5)
