#!/usr/bin/env python3
"""Standalone checks for an installed binary wheel, without requiring Blender.

Run with ``python /path/to/test/test_wheel.py``. Source-suite discovery skips
these checks; the separate fallback job exercises installations without binaries.
"""

import importlib
import importlib.machinery
from pathlib import Path
import sys
import unittest


@unittest.skipUnless(__name__ == '__main__', 'run standalone against an installed wheel')
class TestWheelInstallation(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        package = importlib.import_module('molecular_blender')
        cls.package_dir = Path(package.__file__).resolve().parent
        checkout_dir = Path(__file__).resolve().parents[1] / 'molecular_blender'
        if cls.package_dir == checkout_dir:
            raise AssertionError('Wheel checks must load the installed package, not the checkout')
        # Load all extensions before any test, independently of test ordering.
        cls.extensions = {
            name: importlib.import_module(f'molecular_blender.{name}')
            for name in ('compute_ao_cy', 'marching_cube_cy', 'esp_kernel')
        }
        print(f'Testing installed package: {cls.package_dir}', flush=True)

    def test_imports(self):
        """All compiled extensions must belong to the installed package."""
        for name, module in self.extensions.items():
            with self.subTest(module=name):
                self.assertIsInstance(module.__loader__, importlib.machinery.ExtensionFileLoader)
                self.assertEqual(Path(module.__file__).resolve().parent, self.package_dir)

    def test_basic_functionality(self):
        """A normalized Gaussian and a small closed surface give known results."""
        import numpy as np

        ao_cy = self.extensions['compute_ao_cy']
        center = np.zeros(3, dtype=np.float32)
        coeffs = np.ones(1, dtype=np.float32)
        shell = ao_cy.Shell(center, 0, coeffs, coeffs, 0)
        self.assertAlmostEqual(shell.value(0.0, 0.0, 0.0, coeffs),
                               (2.0 / np.pi) ** 0.75, places=6)

        mc_cy = self.extensions['marching_cube_cy']
        data = np.zeros((3, 3, 3), dtype=np.float32)
        data[1, 1, 1] = 1.0
        result = mc_cy.marching_cube_box(data, [0.5])
        self.assertEqual(len(result), 1)
        self.assertEqual(np.asarray(result[0]).shape, (8, 3, 3))
        self.assertTrue(np.isfinite(result[0]).all())

    def test_reference_esp(self):
        from installed_esp import check_water_esp
        check_water_esp()


def main():
    # CWD and PYTHONPATH must not let checkout code mask a broken installation.
    checkout = Path(__file__).resolve().parents[1]
    sys.path[:] = [entry for entry in sys.path if Path(entry or '.').resolve() != checkout]
    suite = unittest.defaultTestLoader.loadTestsFromTestCase(TestWheelInstallation)
    result = unittest.TextTestRunner(verbosity=2).run(suite)
    return 0 if result.wasSuccessful() else 1


if __name__ == '__main__':
    sys.exit(main())
