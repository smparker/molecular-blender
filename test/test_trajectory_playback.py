import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import numpy as np

from molecular_blender.bond_perception import perceive_bond_orders, BondPerceptionError
from molecular_blender.importers import molecule_from_file
from molecular_blender.molecular_data import Molecule
from molecular_blender.options import get_frame_mapper, validated_options

try:
    import bpy
except ImportError:
    bpy = None


def write_reaction(path, invalid=False):
    frames = []
    for pair in (('H 0 0 4', 'H .72 0 4'), ('H 0 0 -1.09', 'H 2.21 0 0')):
        frames.append('\n'.join(['6', 'carbonyl hydrogenation', 'C 0 0 0', 'O 1.25 0 0',
                                 'H -.58 .9 0', 'H -.58 -.9 0', *pair]))
    if invalid:
        frames.append(frames[1].replace('H 2.21 0 0', 'H 20 0 0'))
    path.write_text('\n'.join(frames)+'\n')


class TestTrajectoryModel(unittest.TestCase):
    def test_physical_units_and_sparse_source_times(self):
        for unit, step in [('fs', 250.), ('ps', .25)]:
            opts = validated_options(dict(trajectory_timing='physical', time_unit=unit, time_step=step,
                playback_ps_per_second=2., playback_fps=30, keystride=11))
            self.assertEqual(get_frame_mapper(opts, [2, 5])(1), 19.75)
        with self.assertRaises(ValueError):
            validated_options(dict(playback_ps_per_second=0.))

    def test_fragment_charges_ambiguity_and_global_budget(self):
        stats = {}
        orders, charges = perceive_bond_orders(['cl','n','h','h','h','h'],
            [(1,2),(1,3),(1,4),(1,5)], diagnostics=stats)
        self.assertEqual(charges, [-1,1,0,0,0,0])
        self.assertEqual((stats['ambiguity'], stats['fragments'], stats['charge_cost']), ('unique',2,2))
        ring = [(i,(i+1)%6) for i in range(6)] + [(i,i+6) for i in range(6)]
        symbols = ['c']*6 + ['h']*6
        orders, _ = perceive_bond_orders(symbols, ring, diagnostics=stats)
        self.assertEqual(stats['ambiguity'], 'multiple')
        perceive_bond_orders(symbols, ring, fixed_orders={0:orders[0]}, diagnostics=stats)
        self.assertEqual(stats['ambiguity'], 'unique')
        perceive_bond_orders(symbols, ring, diagnostics=stats, max_iterations=2)
        self.assertEqual(stats['ambiguity'], 'undetermined')
        with self.assertRaises(BondPerceptionError):
            perceive_bond_orders(['h']*4, [(0,1),(2,3)], max_iterations=1)

    def test_solver_failure_names_frame_and_keeps_previous_graph(self):
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder)/'reaction.xyz'
            write_reaction(path, invalid=True)
            settings = validated_options(dict(animate_bonds='dynamic'))
            model = Molecule.from_dict('reaction', molecule_from_file(path,settings))
            model.determine_bonding(settings)
            original = model.bonds
            with self.assertRaisesRegex(BondPerceptionError, 'Input frame 2'):
                model.determine_bonding(dict(settings, perceive_bond_orders=True))
            self.assertIs(model.bonds, original)


@unittest.skipUnless(bpy, 'requires Blender')
class TestTrajectoryScene(unittest.TestCase):
    def test_order_changes_both_halves_and_survives_scrubbing(self):
        from molecular_blender.plotter import BlendMoleculeSetup
        from molecular_blender.scene_state import read_record, owned_objects
        from molecular_blender.scene_resources import remove_molecule
        scene = bpy.context.scene
        previous = (scene.render.fps, scene.render.fps_base, scene.frame_current, scene.frame_end)
        def restore():
            scene.render.fps, scene.render.fps_base, frame, scene.frame_end = previous
            scene.frame_set(frame)
        self.addCleanup(restore)
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder)/'reaction.xyz'
            write_reaction(path)
            model, _, _ = BlendMoleculeSetup(bpy.context,str(path), animate_bonds='dynamic',
                perceive_bond_orders=True, draw_bond_orders=True, ignore_hydrogen=False,
                plot_style='sticks', isovalues='', trajectory_timing='physical', time_step=.5,
                time_unit='ps', playback_ps_per_second=1., playback_fps=30)
            owner = model.scene_owner
            self.addCleanup(remove_molecule, owner)
            record = read_record(owner)
            self.assertEqual(record['bond_frames'], [1,16])
            self.assertEqual(scene.render.fps, 30)
            carbonyl = next(b for b in record['bonds'] if b['frame_orders'] == [2,1])
            curves = [o for o in owned_objects(owner,'bond') if o['mb_id'] == carbonyl['id']]
            self.assertEqual(len(curves), 4)
            for frame, order in ((1,2),(16,1),(8,2)):
                scene.frame_set(frame)
                visible = [c for c in curves if not c.hide_render and not c.hide_viewport]
                self.assertEqual([c['mb_order'] for c in visible], [order, order])
                self.assertTrue(all((c.data.bevel_object is None) == (order == 1) for c in visible))
            for i in (0,1):
                self.assertEqual(sum(a['frame_charges'][i] for a in record['atoms']), 0)
            self.assertEqual([r['fragments'] for r in record['perception']], [2,1])

    def test_packed_cube_frames_match_resident_meshes(self):
        from molecular_blender.plotter import BlendMoleculeSetup, draw_surfaces
        from molecular_blender.scene_resources import remove_molecule
        from molecular_blender.scene_state import owned_objects
        from molecular_blender.trajectory_meshes import ARCHIVE
        from test_esp import write_cube
        scene = bpy.context.scene
        self.addCleanup(scene.frame_set, scene.frame_current)
        self.addCleanup(setattr, scene.render, 'use_lock_interface', scene.render.use_lock_interface)
        with tempfile.TemporaryDirectory() as folder:
            folder = Path(folder)
            x,y,z = np.indices((9,9,9))*.5-2
            entries = []
            for i, t in enumerate([0,3,8]):
                path = folder/f'{i}.cube'
                write_cube(path, np.exp(-((x-.2*i)**2+y*y+z*z)) if i!=1 else np.zeros_like(x))
                entries.append(dict(filename=path.name,frame=t))
            manifest = folder/'frames.json'
            manifest.write_text(json.dumps(dict(molecules=entries)))
            modes = {}
            for mode in ('resident','packed'):
                model, styler, opts = BlendMoleculeSetup(bpy.context,str(manifest),trajectory_meshes=mode,
                    volume='density',cube_field_type='density', cumulative=False,isovalues='.3',remesh=0)
                self.addCleanup(remove_molecule, model.scene_owner)
                if mode == 'packed':
                    before = set(bpy.data.meshes)
                    prior_marker = scene.get(ARCHIVE)
                    with patch('molecular_blender.scene_meshes.create_surface_meshes', side_effect=RuntimeError('failed')):
                        with self.assertRaisesRegex(RuntimeError, 'failed'):
                            draw_surfaces(model, styler, bpy.context, opts)
                    self.assertEqual(set(bpy.data.meshes), before)
                    self.assertEqual(scene.get(ARCHIVE), prior_marker)
                draw_surfaces(model,styler,bpy.context,opts)
                modes[mode] = owned_objects(model.scene_owner,'surface')
            proxy, = modes['packed']
            self.assertIn(ARCHIVE,proxy.data)
            for t in (1,4,9,6,1):
                scene.frame_set(t)
                expected = [tuple(v.co) for o in modes['resident'] if not o.hide_render for v in o.data.vertices]
                self.assertEqual(expected,[tuple(v.co) for v in proxy.data.vertices])
                self.assertEqual(len(proxy.data.polygons),
                                 sum(len(o.data.polygons) for o in modes['resident'] if not o.hide_render))
