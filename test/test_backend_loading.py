"""Missing optional modules may fall back; broken modules must remain visible."""

from pathlib import Path
import subprocess
import sys
import unittest

ROOT = Path(__file__).resolve().parents[1]


@unittest.skipIf('bpy' in sys.modules, 'requires a standalone Python interpreter')
class TestBackendLoading(unittest.TestCase):
    def load_with_failure(self, target, exception):
        code = f'''
import importlib.abc
import sys
sys.path.insert(0, {str(ROOT)!r})
class Failure(importlib.abc.MetaPathFinder):
    def find_spec(self, fullname, path=None, target=None):
        if fullname == {target!r}:
            raise {exception}
sys.meta_path.insert(0, Failure())
import molecular_blender
'''
        return subprocess.run([sys.executable, '-c', code], cwd=ROOT,
                              capture_output=True, text=True)

    def test_missing_extensions_fall_back(self):
        for target in ('molecular_blender.compute_ao_cy', 'molecular_blender.marching_cube_cy'):
            with self.subTest(target=target):
                result = self.load_with_failure(target, f'ModuleNotFoundError(name={target!r})')
                self.assertEqual(result.returncode, 0, result.stderr)

    def test_broken_extensions_do_not_fall_back(self):
        for target in ('molecular_blender.compute_ao_cy', 'molecular_blender.marching_cube_cy'):
            for exception in ("ImportError('broken binary')",
                              "ModuleNotFoundError('missing internal dependency', name='internal_dependency')"):
                with self.subTest(target=target, exception=exception):
                    result = self.load_with_failure(target, exception)
                    self.assertNotEqual(result.returncode, 0)
                    self.assertIn('Error', result.stderr)

    def test_missing_bpy_is_optional_but_broken_bpy_is_not(self):
        missing = self.load_with_failure('bpy', "ModuleNotFoundError(name='bpy')")
        self.assertEqual(missing.returncode, 0, missing.stderr)
        broken = self.load_with_failure('bpy', "ImportError('broken Blender dependency')")
        self.assertNotEqual(broken.returncode, 0)
        self.assertIn('broken Blender dependency', broken.stderr)
