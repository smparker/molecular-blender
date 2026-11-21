"""Profile a bonded carbon ball-and-stick import with and without editable endpoint hooks.

blender -b --factory-startup --python-exit-code 1 --python benchmarks/bonded_import.py
Append -- --details for call profiles of molecule tagging and an isolated material-slot timing.
"""
from collections import defaultdict
from contextlib import ExitStack
import cProfile
from math import sqrt
from pathlib import Path
import pstats
from statistics import median
import sys
import tempfile
from time import perf_counter
from unittest.mock import patch

import bpy

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from molecular_blender import plotter, scene_geometry, scene_state


def write_sheet(path, width, height):
    # Four sites in a rectangular graphene cell, C–C distance in angstroms.
    d = 1.42
    h = sqrt(3) / 2
    sites = [(0, 0), (1, 0), (1.5, h), (2.5, h)]
    rows = [f'C {(3*i+x)*d:.6f} {(2*h*j+y)*d:.6f} 0'
            for i in range(width) for j in range(height) for x, y in sites]
    path.write_text(f'{len(rows)}\nGraphene patch\n' + '\n'.join(rows) + '\n')
    return [tuple(map(float, row.split()[1:])) for row in rows]


def measure(path, hooks, atom_count):
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
        def timed(name, function):
            def call(*args, **kwargs):
                before = perf_counter()
                try:
                    return function(*args, **kwargs)
                finally:
                    phases[name] += perf_counter() - before
            return call
        class Proxy:
            def __init__(self, base, **attrs):
                self.base = base
                self.__dict__.update(attrs)
            def __getattr__(self, name):
                return getattr(self.base, name)
        mesh_ops = Proxy(bpy.ops.mesh, primitive_uv_sphere_add=timed('sphere primitive', bpy.ops.mesh.primitive_uv_sphere_add))
        ops = Proxy(bpy.ops, mesh=mesh_ops)
        with ExitStack() as stack:
            stack.enter_context(patch.object(scene_geometry, 'Timer', PhaseTimer))
            stack.enter_context(patch.object(scene_geometry, 'bpy', Proxy(bpy, ops=ops)))
            for module, name in [(plotter, 'make_materials'), (scene_state, 'persist_molecule'),
                                 (plotter, 'molecule_from_file'), (plotter.Molecule, 'determine_bonding')]:
                stack.enter_context(patch.object(module, name, timed(name, getattr(module, name))))
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
    assert len(model.atoms) == atom_count and len(bonds) > atom_count
    phases["bonds"] = len(bonds)
    assert count == (2 * len(bonds) if hooks == 'on' else 0)
    record = scene_state.read_record(model.scene_owner)
    assert len(record['atoms']) == atom_count and len(record['bonds']) == len(bonds)
    for atom, entry in zip(model.atoms, record['atoms']):
        obj = bpy.data.objects[atom.name]
        assert obj[scene_state.ID] == entry['id'] and obj['mb_atom_index'] == atom.index
        assert obj[scene_state.MOLECULE_ID] == record['id'] and obj[scene_state.ROLE] == 'atom'
        assert (obj.matrix_world.translation - atom.position).length < 1.e-5
    for bond, curve, entry in zip(model.bonds, bonds, record['bonds']):
        assert curve[scene_state.ID] == entry['id'] and list(curve['mb_endpoints']) == entry['atoms']
        assert curve[scene_state.MOLECULE_ID] == record['id'] and curve[scene_state.ROLE] == 'bond'
        for point, atom in zip(curve.data.splines[0].bezier_points, (bond.iatom, bond.jatom)):
            assert (curve.matrix_world @ point.co - atom.position).length < 1.e-5
    return phases


def atom_geometry(positions, nodes):
    # Isolate the representation cost using identical sphere geometry, material,
    # positions and collection. Neither prototype constructs bonds or scene records.
    bpy.ops.wm.read_factory_settings(use_empty=True)
    scene = bpy.context.scene
    bpy.ops.mesh.primitive_uv_sphere_add(radius=.425)
    sphere = bpy.context.object
    sphere.data.materials.append(bpy.data.materials.new('Carbon'))
    scene.collection.objects.unlink(sphere)
    source_mesh = sphere.data
    start = perf_counter()
    if nodes:
        mesh = bpy.data.meshes.new('Sites')
        mesh.from_pydata(positions, [], [])
        obj = bpy.data.objects.new('Atoms', mesh)
        scene.collection.objects.link(obj)
        group = bpy.data.node_groups.new('Atom instances', 'GeometryNodeTree')
        for direction in ('INPUT', 'OUTPUT'):
            group.interface.new_socket(name='Geometry', in_out=direction, socket_type='NodeSocketGeometry')
        source = group.nodes.new('GeometryNodeObjectInfo')
        source.inputs['Object'].default_value = sphere
        instance = group.nodes.new('GeometryNodeInstanceOnPoints')
        group.links.new(group.nodes.new('NodeGroupInput').outputs['Geometry'], instance.inputs['Points'])
        group.links.new(source.outputs['Geometry'], instance.inputs['Instance'])
        group.links.new(instance.outputs['Instances'], group.nodes.new('NodeGroupOutput').inputs['Geometry'])
        obj.modifiers.new('Atoms', 'NODES').node_group = group
    else:
        for point in positions:
            obj = bpy.data.objects.new('Atom', source_mesh)
            obj.location = point
            scene.collection.objects.link(obj)
    built = perf_counter()
    bpy.context.view_layer.update()
    evaluated = perf_counter()
    dg = bpy.context.evaluated_depsgraph_get()
    locations = [tuple(round(v, 4) for v in it.matrix_world.translation) for it in dg.object_instances
                 if not nodes or it.is_instance]
    assert sorted(locations) == sorted(tuple(round(v, 4) for v in xyz) for xyz in positions)
    return [(built-start)*1000, (evaluated-built)*1000, (evaluated-start)*1000]


def detail_profile(path, positions):
    for hooks in ('off', 'on'):
        profile = cProfile.Profile()
        persist = scene_state.persist_molecule
        with patch.object(scene_state, 'persist_molecule',
                          lambda *args: profile.runcall(persist, *args)):
            phases = measure(path, hooks, len(positions))
        print('Molecule record and tagging, hooks', hooks)
        pstats.Stats(profile).strip_dirs().sort_stats('tottime').print_stats(8)
    # Reproduce the slot assignments on fresh, unlinked datablocks, at the
    # same point in their lifecycle as in the geometry importer.
    bpy.ops.wm.read_factory_settings(use_empty=True)
    material = bpy.data.materials.new('Carbon').name
    times = []
    for trial in range(4):
        sphere = bpy.data.meshes.new('Atom template')
        curves = [bpy.data.curves.new('Bond', 'CURVE') for _ in range(phases['bonds'])]
        start = perf_counter()
        sphere.materials.append(bpy.data.materials[material])
        for curve in curves:
            curve.materials.append(bpy.data.materials[material])
        times.append((perf_counter() - start) * 1000)
        bpy.data.batch_remove([sphere, *curves])
    print(f"Material-slot assignment for one atom template and {phases['bonds']} curves: "
          f"{median(times[1:]):.2f} ms")


with tempfile.TemporaryDirectory() as folder:
    path = Path(folder) / 'sheet.xyz'
    if '--details' in sys.argv:
        detail_profile(path, write_sheet(path, 24, 16))
        sys.exit(0)
    for width, height in ((12, 8), (24, 16)):
        positions = write_sheet(path, width, height)
        samples = {'off': [], 'on': []}
        for run in range(4):
            for hooks in (('off', 'on') if run % 2 == 0 else ('on', 'off')):
                phases = measure(path, hooks, len(positions))
                if run:
                    samples[hooks].append(phases)
        print(f"{len(positions)} atoms, {phases['bonds']} bonds · median of three imports (ms)")
        print(f'{"Phase":32} {"Hooks off":>10} {"Hooks on":>10}')
        for label in ('molecule_from_file', 'determine_bonding', 'make_materials', 'sphere primitive',
                      'atom creation', 'bond creation', 'link objects', 'update scene', 'set hooks',
                      'set parentage', 'persist_molecule', 'depsgraph (overlaps stages)', 'total'):
            times = [median(s[label] for s in samples[hooks]) * 1000 for hooks in ('off', 'on')]
            print(f'{label:32} {times[0]:10.1f} {times[1]:10.1f}')
        for nodes in (False, True):
            trials = [atom_geometry(positions, nodes) for _ in range(4)]
            print('atom-only prototype', 'Geometry Nodes' if nodes else 'linked objects',
                  'build/update/total ms', *[round(median(column), 2) for column in zip(*trials[1:])], flush=True)
