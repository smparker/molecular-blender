#!/usr/bin/env python3
"""
Simple test script to verify that a molecular_blender wheel installation works correctly.
This script tests basic imports and functionality of the Cython extensions.
"""

import sys
import traceback

def test_imports():
    """Test that all required modules can be imported."""
    print("Testing imports...")

    # Test basic package import
    try:
        import molecular_blender
        print("✓ molecular_blender package imported successfully")
    except ImportError as e:
        print(f"✗ Failed to import molecular_blender: {e}")
        return False

    # Test Cython extension imports
    try:
        import molecular_blender.compute_ao_cy
        print("✓ molecular_blender.compute_ao_cy imported successfully")
    except ImportError as e:
        print(f"✗ Failed to import compute_ao_cy: {e}")
        traceback.print_exc()
        return False

    try:
        import molecular_blender.marching_cube_cy
        print("✓ molecular_blender.marching_cube_cy imported successfully")
    except ImportError as e:
        print(f"✗ Failed to import marching_cube_cy: {e}")
        traceback.print_exc()
        return False

    return True

def test_basic_functionality():
    """Test basic functionality of Cython extensions."""
    print("\nTesting basic functionality...")

    try:
        import numpy as np
        from molecular_blender.compute_ao_cy import Shell

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
        from molecular_blender.marching_cube_cy import marching_cube_box

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

