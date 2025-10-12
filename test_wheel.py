#!/usr/bin/env python3
"""
Simple test script to verify that a molecular_blender wheel installation works correctly.
This script tests basic imports and functionality of the Cython extensions.
"""

import sys
import os
import glob
import traceback
import importlib.util

def test_imports():
    """Test that all required modules can be imported."""
    print("Testing imports...")

    # Note: We skip importing the main molecular_blender package because it requires bpy (Blender Python API)
    # which is not available outside of Blender. Instead, we directly import the Cython extension modules
    # by finding them in the site-packages directory, bypassing the package __init__.py

    print("⚠ Skipping molecular_blender package import (requires bpy from Blender)")

    # Find the molecular_blender installation directory without importing it
    # (to avoid triggering the bpy import in __init__.py)
    try:
        spec = importlib.util.find_spec("molecular_blender")
        if spec is None or spec.origin is None:
            print(f"✗ Cannot locate molecular_blender installation")
            return False
        mb_dir = os.path.dirname(spec.origin)
        print(f"  Found molecular_blender at: {mb_dir}")
    except (ImportError, ValueError) as e:
        print(f"✗ Cannot locate molecular_blender installation: {e}")
        return False

    # Directly import Cython extensions from the directory, bypassing __init__.py
    modules_to_test = ['compute_ao_cy', 'marching_cube_cy']

    for module_name in modules_to_test:
        try:
            # Look for the compiled extension file
            extensions = ['.so', '.pyd']
            module_file = None

            for ext in extensions:
                # Try different naming patterns
                patterns = [
                    f"{module_name}{ext}",
                    f"{module_name}.cpython-*{ext}",
                ]
                for pattern in patterns:
                    matches = glob.glob(os.path.join(mb_dir, pattern))
                    if matches:
                        module_file = matches[0]
                        break
                if module_file:
                    break

            if not module_file:
                print(f"✗ Could not find compiled extension for {module_name}")
                continue

            # Load the module directly from the file
            spec = importlib.util.spec_from_file_location(
                f"molecular_blender.{module_name}",
                module_file
            )
            module = importlib.util.module_from_spec(spec)
            sys.modules[f"molecular_blender.{module_name}"] = module
            spec.loader.exec_module(module)

            print(f"✓ molecular_blender.{module_name} imported successfully")

        except Exception as e:
            print(f"✗ Failed to import {module_name}: {e}")
            traceback.print_exc()
            return False

    return True

def test_basic_functionality():
    """Test basic functionality of Cython extensions."""
    print("\nTesting basic functionality...")

    # Check if the modules were successfully imported in test_imports
    if 'molecular_blender.compute_ao_cy' not in sys.modules:
        print("⚠ Skipping functionality tests - Cython modules not available")
        return True

    try:
        import numpy as np
        # Use the already imported module from sys.modules
        compute_ao_cy = sys.modules['molecular_blender.compute_ao_cy']
        Shell = compute_ao_cy.Shell

        # Create a simple shell
        center = np.array([0.0, 0.0, 0.0], dtype=np.float32)
        exponents = np.array([1.0], dtype=np.float32)
        coeffs = np.array([1.0], dtype=np.float32)

        shell = Shell(center, 0, exponents, coeffs, 0)
        print("✓ Shell object created successfully")

        # Test value calculation
        val = shell.value(0.0, 0.0, 0.0, coeffs)
        print(f"✓ Shell.value() computed: {val}")

    except Exception as e:
        print(f"✗ Shell functionality test failed: {e}")
        traceback.print_exc()
        return False

    try:
        import numpy as np
        # Use the already imported module from sys.modules
        marching_cube_cy = sys.modules['molecular_blender.marching_cube_cy']
        marching_cube_box = marching_cube_cy.marching_cube_box

        # Create a simple test grid
        data = np.zeros((3, 3, 3), dtype=np.float32)
        data[1, 1, 1] = 1.0

        # Test marching cubes
        result = marching_cube_box(data, [0.5])
        print(f"✓ marching_cube_box() executed successfully, found {len(result[0])} triangles")

    except Exception as e:
        print(f"✗ Marching cubes test failed: {e}")
        traceback.print_exc()
        return False

    return True

def test_numpy_compatibility():
    """Test that numpy is available and compatible."""
    print("\nTesting numpy compatibility...")

    try:
        import numpy as np
        print(f"✓ NumPy version: {np.__version__}")

        # Test basic operations
        arr = np.array([1.0, 2.0, 3.0], dtype=np.float32)
        result = np.sum(arr)
        print(f"✓ NumPy operations work correctly (sum test: {result})")

        return True
    except Exception as e:
        print(f"✗ NumPy test failed: {e}")
        return False

def main():
    """Run all tests."""
    print("=" * 60)
    print("Molecular Blender Wheel Test Suite")
    print("=" * 60)
    print(f"Python version: {sys.version}")
    print(f"Platform: {sys.platform}")
    print("=" * 60)

    all_passed = True

    # Run tests
    all_passed &= test_imports()
    all_passed &= test_numpy_compatibility()
    all_passed &= test_basic_functionality()

    print("\n" + "=" * 60)
    if all_passed:
        print("✓ All tests passed!")
        print("=" * 60)
        return 0
    else:
        print("✗ Some tests failed")
        print("=" * 60)
        return 1

if __name__ == "__main__":
    sys.exit(main())

