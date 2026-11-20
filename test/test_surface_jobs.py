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
