"""Regression tests for invoking the suite through shell and Blender arguments."""

import importlib.util
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('mb_test_runner', ROOT / 'test/run_tests.py')
runner = importlib.util.module_from_spec(spec)
spec.loader.exec_module(runner)


class TestRunnerArguments(unittest.TestCase):
    def test_blender_arguments(self):
        self.assertEqual(runner.blender_extra_args([
            '/Applications/Custom App', '--background', '--python', 'tests.py',
            '--', 'test_bond_orders.TestBondOrders']), ['test_bond_orders.TestBondOrders'])
        self.assertEqual(runner.blender_extra_args(['blender', '--python', 'tests.py']), [])
        self.assertEqual(runner.blender_extra_args(['run_tests.py', 'test_util']), ['test_util'])

    @unittest.skipUnless(shutil.which('bash'), 'requires bash')
    def test_shell_preserves_paths_and_arguments(self):
        with tempfile.TemporaryDirectory(prefix='molecular blender ') as temporary:
            root = Path(temporary)
            script = root / 'run tests.sh'
            shutil.copyfile(ROOT / 'run_test.sh', script)
            blender = root / 'Blender App'
            blender.write_text('#!/bin/sh\nprintf \'%s\\0\' "$@" > "$RUNNER_ARGS"\n'
                               'pwd > "$RUNNER_CWD"\n')
            blender.chmod(0o755)
            arguments = root / 'arguments'
            cwd = root / 'cwd'
            environment = dict(os.environ, RUNNER_ARGS=str(arguments), RUNNER_CWD=str(cwd))
            for mode, python_script in (('-t', 'test/run_tests.py'), ('-p', 'marching_profile.py')):
                with self.subTest(mode=mode):
                    subprocess.run(['bash', str(script), mode, '-b', str(blender),
                                    'name with spaces', '*.literal'], cwd='/',
                                   env=environment, check=True, capture_output=True)
                    actual = arguments.read_bytes().decode().split('\0')[:-1]
                    self.assertEqual(actual, ['--background', '--factory-startup',
                                             '--python-exit-code', '1', '--python',
                                             str(root / python_script), '--',
                                             'name with spaces', '*.literal'])
                    self.assertEqual(Path(cwd.read_text().strip()).resolve(), root.resolve())



@unittest.skipUnless(runner.bpy, 'requires Blender')
class TestSceneIsolation(unittest.TestCase):
    def test_cleanup_after_failed_setup_before_next_test(self):
        import io
        import bpy
        original = bpy.context.scene
        baseline = set(bpy.data.user_map())

        class FailingSetup(unittest.TestCase):
            def setUp(self):
                bpy.context.scene.frame_set(47)
                mesh = bpy.data.meshes.new('Leaked mesh')
                obj = bpy.data.objects.new('Leaked object', mesh)
                bpy.context.scene.collection.objects.link(obj)
                bpy.data.materials.new('Orphan material')
                raise RuntimeError('setup failed')

            def runTest(self):
                pass

        result = runner.GroupedTestResult(io.StringIO())
        FailingSetup().run(result)
        self.assertEqual(len(result.errors), 1)
        self.assertIn('setup failed', result.errors[0][1])
        self.assertEqual(bpy.context.scene, original)
        self.assertEqual(set(bpy.data.user_map()), baseline)


if __name__ == '__main__':
    unittest.main()
