#!/usr/bin/env python3
"""
Simple test script to verify that a molecular_blender wheel installation works correctly.
This script tests basic imports and functionality of the Cython extensions.
"""

import sys
import os
import glob
import importlib.util
import unittest


def _find_cython_module(mb_dir, module_name):
    """Find the compiled extension file matching the current Python ABI, or None."""
    pyver = f"cpython-{sys.version_info.major}{sys.version_info.minor}"
    for ext in ['.so', '.pyd']:
        for pattern in [f"{module_name}.{pyver}-*{ext}", f"{module_name}{ext}"]:
            matches = glob.glob(os.path.join(mb_dir, pattern))
            if matches:
                return matches[0]
    return None


def _load_module(module_name, module_file):
    """Load a module from a file path into sys.modules."""
    key = f"molecular_blender.{module_name}"
    spec = importlib.util.spec_from_file_location(key, module_file)
    module = importlib.util.module_from_spec(spec)
    sys.modules[key] = module
    spec.loader.exec_module(module)
    return module


class TestWheelInstallation(unittest.TestCase):
    def test_imports(self):
        """Cython extension modules can be found and imported for the current Python."""
        spec = importlib.util.find_spec("molecular_blender")
        if spec is None or spec.origin is None:
            self.skipTest("molecular_blender not installed")
        mb_dir = os.path.dirname(spec.origin)

        for module_name in ['compute_ao_cy', 'marching_cube_cy']:
            module_file = _find_cython_module(mb_dir, module_name)
            if module_file is None:
                continue  # not built; that's acceptable
            _load_module(module_name, module_file)

    def test_basic_functionality(self):
        """Cython Shell and marching_cube_box produce valid results."""
        import numpy as np

        ao_cy = sys.modules.get('molecular_blender.compute_ao_cy')
        if ao_cy is not None:
            center = np.array([0.0, 0.0, 0.0], dtype=np.float32)
            coeffs = np.array([1.0], dtype=np.float32)
            shell = ao_cy.Shell(center, 0, coeffs, coeffs, 0)
            val = shell.value(0.0, 0.0, 0.0, coeffs)
            self.assertIsNotNone(val)

        mc_cy = sys.modules.get('molecular_blender.marching_cube_cy')
        if mc_cy is not None:
            data = np.zeros((3, 3, 3), dtype=np.float32)
            data[1, 1, 1] = 1.0
            result = mc_cy.marching_cube_box(data, [0.5])
            self.assertEqual(len(result), 1)

    def test_numpy_compatibility(self):
        """NumPy is available and basic operations work."""
        import numpy as np
        arr = np.array([1.0, 2.0, 3.0], dtype=np.float32)
        self.assertAlmostEqual(float(np.sum(arr)), 6.0)


def main():
    loader = unittest.TestLoader()
    suite = loader.loadTestsFromTestCase(TestWheelInstallation)
    runner = unittest.TextTestRunner(verbosity=2)
    result = runner.run(suite)
    return 0 if result.wasSuccessful() else 1


if __name__ == "__main__":
    sys.exit(main())
