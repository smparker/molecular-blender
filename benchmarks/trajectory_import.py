"""Time trajectory import and its atom keyframing step in two scene sizes.

blender -b --factory-startup --python-exit-code 1 --python benchmarks/trajectory_import.py
"""
from pathlib import Path
from statistics import median
import sys
from time import perf_counter
from unittest.mock import patch

import bpy

root = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(root))
from molecular_blender import plotter

animate = plotter.AnimateMolecule
for extras in (0, 1000):
    imports, keyframes = [], []
    for trial in range(3):
        bpy.ops.wm.read_factory_settings(use_empty=True)
        scene = bpy.context.scene
        for i in range(extras):
            scene.collection.objects.link(bpy.data.objects.new(f'Figure object {i}', None))
        bpy.context.view_layer.update()

        def timed_animation(*args):
            start = perf_counter()
            animate(*args)
            keyframes.append(perf_counter() - start)

        with patch.object(plotter, 'AnimateMolecule', timed_animation):
            start = perf_counter()
            plotter.BlendMoleculeSetup(bpy.context, str(root / 'examples/benzopentyl_twist.xyz'),
                                       hook_atoms='off', ignore_hydrogen=False, isovalues='')
            imports.append(perf_counter() - start)
    print(f'{extras} extra objects: import {median(imports)*1000:.2f} ms, '
          f'keyframing {median(keyframes)*1000:.2f} ms (medians of 3)')
