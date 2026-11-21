"""Process lifetime and Blender modal cancellation regressions."""
import os
import shutil
import subprocess
import time
from types import SimpleNamespace
import unittest
from unittest.mock import Mock, patch

from molecular_blender.surface_jobs import SurfaceJob, cancel_all_jobs

try:
    import bpy
except ImportError:
    bpy = None


class TestSurfaceJobLifetime(unittest.TestCase):
    def test_cancel_reaps_process_and_removes_files(self):
        popen = subprocess.Popen
        python = shutil.which('python3') or shutil.which('python')
        with patch('molecular_blender.surface_jobs.subprocess.Popen',
                   side_effect=lambda *args, **kw: popen(
                       [python, '-c', 'import time; time.sleep(60)'], **kw)):
            job = SurfaceJob('unused', {}, {})
        directory = job.directory
        self.assertTrue(job.running)
        start = time.monotonic()
        cancel_all_jobs()
        self.assertLess(time.monotonic() - start, 3)
        self.assertFalse(job.running)
        self.assertIsNotNone(job.process.returncode)
        self.assertFalse(directory.exists())
        with self.assertRaisesRegex(RuntimeError, 'cancelled'):
            job.result()
        job.cancel()  # repeated cleanup is harmless

    def test_launch_failure_cleans_temporary_directory(self):
        from tempfile import TemporaryDirectory
        temporary = TemporaryDirectory(prefix='surface-launch-test-')
        directory = temporary.name
        with patch('molecular_blender.surface_jobs.tempfile.TemporaryDirectory', return_value=temporary):
            with self.assertRaises(OSError):
                SurfaceJob('/nonexistent/molecular-blender', {}, {})
        self.assertFalse(os.path.exists(directory))


@unittest.skipUnless(bpy, 'requires Blender')
class TestSurfaceWorker(unittest.TestCase):
    def test_worker_matches_synchronous_result(self):
        import numpy as np
        from molecular_blender.importers import molecule_from_file
        from molecular_blender.molecule import Molecule
        from molecular_blender.plotter import default_options, process_options, compute_surface_vertex_sets
        data = molecule_from_file('examples/fake-atom.molden', {})
        options = default_options(dict(orbital='1', cumulative=False, isovalues='0.05', resolution=0.3))
        process_options(options)
        expected, mapping = compute_surface_vertex_sets(Molecule.from_dict('reference', data), options, indexed=True)
        job = SurfaceJob(bpy.app.binary_path, data, options)
        try:
            job.process.wait(timeout=60)
            actual, actual_mapping = job.result()
            self.assertEqual(actual_mapping, mapping)
            self.assertEqual(len(actual), len(expected))
            for a, b in zip(actual, expected):
                np.testing.assert_array_equal(a['faces'], b['faces'])
                np.testing.assert_array_equal(a['vertices'], b['vertices'])
            self.assertEqual(job.progress()['completed'], job.progress()['total'])
        finally:
            job.close()
        self.assertFalse(job.directory.exists())

    def test_molden_frames_keep_channels_levels_and_legacy_result_order(self):
        import numpy as np
        from molecular_blender.importers import molecule_from_file
        from molecular_blender.molecule import Molecule
        from molecular_blender.plotter import (compute_surface_vertex_sets, default_options,
                                              iter_surface_frames, process_options)
        model = Molecule.from_dict('reference', molecule_from_file('examples/water-sto-3g.molden', {}))
        options = default_options(dict(orbital='lumo,homo', cumulative=True, isovalues='.65',
                                       resolution=.35, keystride=2))
        process_options(options)
        reference, _ = compute_surface_vertex_sets(model, options, indexed=True)
        model.orbitals_trajectory = [model.orbitals]
        model.frame_numbers = [0, 4]
        batches = list(iter_surface_frames(model, options, indexed=True))
        self.assertEqual([list(mapping) for _, mapping in batches], [[1], [9]])
        for i, (fields, mapping) in enumerate(batches):
            self.assertEqual(len(fields), len(reference))
            for field, expected in zip(fields, reference):
                self.assertIn(field['name'], mapping[(1, 9)[i]])
                self.assertEqual(field.get('series', field['name']), expected['name'])
                self.assertEqual(field['isovalue'], expected['isovalue'])
                np.testing.assert_array_equal(field['vertices'], expected['vertices'])
                np.testing.assert_array_equal(field['faces'], expected['faces'])
        collected, _ = compute_surface_vertex_sets(model, options, indexed=True)
        count_lumo = sum(s['name'] == 'lumo' for s in reference) * 2
        count_homo = sum(s['name'] == 'homo' for s in reference) * 2
        self.assertEqual([s.get('series', s['name']) for s in collected],
                         ['lumo'] * count_lumo + ['homo'] * count_homo)

    def test_streamed_trajectory_advances_in_bounded_batches_and_rolls_back(self):
        import json
        from pathlib import Path
        from tempfile import TemporaryDirectory
        from molecular_blender.plotter import BlendMoleculeSetup, SurfaceMeshStage
        from molecular_blender.scene_state import owned_objects, persist_surfaces, read_record
        from molecular_blender.surface_operator import SurfaceJobOperator

        original = bpy.context.window.scene
        scene = bpy.data.scenes.new('Streamed surfaces')
        bpy.context.window.scene = scene
        def ids():
            return {item for kind in ('objects', 'meshes', 'curves', 'collections', 'materials', 'actions')
                    for item in getattr(bpy.data, kind)}
        original_ids = ids()
        def cleanup():
            bpy.context.window.scene = original
            bpy.data.batch_remove(ids() - original_ids)
            bpy.data.scenes.remove(scene)
        self.addCleanup(cleanup)
        with TemporaryDirectory() as directory:
            cube = Path(directory) / 'signed.cube'
            samples = [0] * 125
            samples[37], samples[87] = 1, -1
            cube.write_text('signed\nfield\n1 0 0 0\n-5 .2 0 0\n-5 0 .2 0\n-5 0 0 .2\n'
                            '6 0 0 0 0\n' + ' '.join(map(str, samples)))
            term = dict(filename=cube.name, coefficient=1)
            filled = dict(units='au', linear_combination=[term])
            empty = dict(units='au', linear_combination=[term, dict(term, coefficient=-1)])
            path = cube.with_suffix('.json')
            path.write_text(json.dumps(dict(molecules=[dict(item, frame=2*i) for i, item in
                                                        enumerate((filled, empty, filled, filled))])))
            model, styler, options = BlendMoleculeSetup(bpy.context, str(path), cumulative=False,
                                                       isovalues='.5,-.5', remesh=0)
        before = ids()
        obj = model.scene_owner
        record = read_record(obj)
        timer = SimpleNamespace(type='TIMER')
        for outcome in ('escape', 'disable', 'publish_failure', 'finish'):
            with self.subTest(outcome=outcome):
                stage = SurfaceMeshStage(model, styler, options)
                job = SurfaceJob(bpy.app.binary_path, model.surface_input, options,
                                 stream_frames=True, on_close=stage.discard)
                self.addCleanup(job.close)
                op = SurfaceJobOperator()
                op._job, op._stage = job, stage
                op._owner, op._collection = obj, model.scene_collection
                op._cancelled, op._timer = False, None
                op.report = Mock()
                deadline = time.monotonic() + 60
                while not (job.directory / 'frame-1.pkl').exists():
                    self.assertTrue(job.running)
                    self.assertLess(time.monotonic(), deadline)
                    time.sleep(.02)
                self.assertEqual(len(list(job.directory.glob('frame-*.pkl'))), 2)
                self.assertEqual(op.modal(bpy.context, timer), {'RUNNING_MODAL'})
                self.assertEqual(len(stage.definitions), 2)
                self.assertTrue(all(not d['mesh_data'].users for d in stage.definitions))
                self.assertEqual(read_record(obj), record)
                self.assertFalse(owned_objects(obj, 'surface'))
                if outcome == 'disable':
                    cancel_all_jobs()
                if outcome in ('escape', 'disable'):
                    event = SimpleNamespace(type='ESC') if outcome == 'escape' else timer
                    self.assertEqual(op.modal(bpy.context, event), {'CANCELLED'})
                else:
                    from contextlib import nullcontext
                    def fail_after_save(*args):
                        persist_surfaces(*args)
                        raise ValueError('Failed after saving surface definitions')
                    hook = (patch('molecular_blender.scene_state.persist_surfaces', fail_after_save)
                            if outcome == 'publish_failure' else nullcontext())
                    with hook:
                        result = {'RUNNING_MODAL'}
                        while result == {'RUNNING_MODAL'}:
                            self.assertLess(time.monotonic(), deadline)
                            result = op.modal(bpy.context, timer)
                            time.sleep(.01)
                    self.assertEqual(result, {'CANCELLED'} if outcome == 'publish_failure' else {'FINISHED'})
                self.assertFalse(job.directory.exists())
                if outcome != 'finish':
                    self.assertEqual(ids(), before)
                    self.assertEqual(read_record(obj), record)

        objects = owned_objects(obj, 'surface')
        self.assertEqual(len(objects), 6)
        self.assertEqual(len({s['group'] for s in read_record(obj)['surfaces']}), 1)
        for timepoint, count in ((1, 2), (3, 0), (4, 0), (5, 2), (7, 2)):
            scene.frame_set(timepoint)
            self.assertEqual(sum(not o.hide_viewport and not o.hide_render for o in objects), count)

    def test_worker_reports_errors(self):
        job = SurfaceJob(bpy.app.binary_path, {}, {})
        try:
            job.process.wait(timeout=60)
            with self.assertRaisesRegex(RuntimeError, 'atoms'):
                job.result()
        finally:
            job.close()

    def test_escape_and_deleted_owner_discard_results(self):
        from molecular_blender import MolecularBlender
        for event, owner_exists in [('ESC', True), ('TIMER', False)]:
            job = Mock()
            wm = Mock()
            context = SimpleNamespace(window_manager=wm, workspace=Mock())
            op = SimpleNamespace(_cancelled=False, _job=job, _timer='timer',
                                 _owner_exists=lambda: owner_exists)
            op._cleanup = lambda ctx: MolecularBlender._cleanup(op, ctx)
            op.cancel = lambda ctx: MolecularBlender.cancel(op, ctx)
            with patch('molecular_blender.surface_operator.create_surface_meshes') as create:
                self.assertEqual(MolecularBlender.modal(op, context, SimpleNamespace(type=event)), {'CANCELLED'})
                create.assert_not_called()
            job.cancel.assert_called_once()
            self.assertIsNone(op._job)
            self.assertIsNone(op._timer)
            wm.event_timer_remove.assert_called_once_with('timer')
