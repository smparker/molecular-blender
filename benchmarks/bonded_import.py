"""Profile a 384-carbon ball-and-stick import with and without editable endpoint hooks.

blender -b --factory-startup --python-exit-code 1 --python benchmarks/bonded_import.py
"""
from collections import defaultdict
from math import sqrt
from pathlib import Path
from statistics import median
import sys
import tempfile
from time import perf_counter
from unittest.mock import patch

import bpy

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from molecular_blender import plotter


def write_sheet(path):
    # Four sites in a rectangular graphene cell, C–C distance in angstroms.
    d = 1.42
    h = sqrt(3) / 2
    sites = [(0, 0), (1, 0), (1.5, h), (2.5, h)]
    rows = [f'C {(3*i+x)*d:.6f} {(2*h*j+y)*d:.6f} 0'
            for i in range(12) for j in range(8) for x, y in sites]
    path.write_text(f'{len(rows)}\nGraphene patch\n' + '\n'.join(rows) + '\n')


def measure(path, hooks):
    bpy.ops.wm.read_factory_settings(use_empty=True)
    phases = defaultdict(float)

    class PhaseTimer:
        def __init__(self):
            self.last = perf_counter()

        def tick_print(self, label):
            now = perf_counter()
            phases[label] += now - self.last
            self.last = now

    def before_update(*_):
        nonlocal evaluation_start
        evaluation_start = perf_counter()

    def after_update(*_):
        phases['depsgraph (overlaps stages)'] += perf_counter() - evaluation_start

    evaluation_start = 0.
    handlers = bpy.app.handlers
    handlers.depsgraph_update_pre.append(before_update)
    handlers.depsgraph_update_post.append(after_update)
    try:
        with patch.object(plotter, 'Timer', PhaseTimer):
            start = perf_counter()
            model, _, _ = plotter.BlendMoleculeSetup(
                bpy.context, str(path), plot_style='bs', bonds=True, hook_atoms=hooks, isovalues='')
            bpy.context.view_layer.update()
            phases['total'] = perf_counter() - start
    finally:
        handlers.depsgraph_update_pre.remove(before_update)
        handlers.depsgraph_update_post.remove(after_update)
    bonds = [bpy.data.objects[bond.name[0]] for bond in model.bonds]
    count = sum(mod.type == 'HOOK' for obj in bonds for mod in obj.modifiers)
    assert len(model.atoms) == 384 and len(bonds) == 537, (len(model.atoms), len(bonds))
    assert count == (2 * len(bonds) if hooks == 'on' else 0)
    for atom in model.atoms:
        assert (bpy.data.objects[atom.name].matrix_world.translation - atom.position).length < 1.e-5
    for bond, curve in zip(model.bonds, bonds):
        for point, atom in zip(curve.data.splines[0].bezier_points, (bond.iatom, bond.jatom)):
            assert (curve.matrix_world @ point.co - atom.position).length < 1.e-5
    return phases


with tempfile.TemporaryDirectory() as folder:
    path = Path(folder) / 'sheet.xyz'
    write_sheet(path)
    samples = {'off': [], 'on': []}
    for run in range(4):
        for hooks in (('off', 'on') if run % 2 == 0 else ('on', 'off')):
            phases = measure(path, hooks)
            if run:  # Warm up each mode once.
                samples[hooks].append(phases)
    print('384 atoms, 537 bonds · median of three imports (ms)')
    print(f'{"Phase":32} {"Hooks off":>10} {"Hooks on":>10}')
    for label in ('atom creation', 'bond creation', 'link objects', 'update scene',
                  'set hooks', 'set parentage', 'depsgraph (overlaps stages)', 'total'):
        times = [median(s[label] for s in samples[hooks]) * 1000 for hooks in ('off', 'on')]
        print(f'{label:32} {times[0]:10.1f} {times[1]:10.1f}')
