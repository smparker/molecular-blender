"""Native extension definitions shared by Blender source builds and Python wheels."""
import ast
from pathlib import Path
import sys

import Cython
from Cython.Build import cythonize
import numpy as np
from setuptools import Extension, setup


def addon_version():
    # Read metadata without importing the package before its extensions exist.
    for node in ast.parse(Path("molecular_blender/__init__.py").read_text()).body:
        if isinstance(node, ast.Assign) and any(
                isinstance(target, ast.Name) and target.id == "bl_info" for target in node.targets):
            return ".".join(map(str, ast.literal_eval(node.value)["version"]))
    raise ValueError("Missing bl_info in molecular_blender/__init__.py")


# The kernels are serial C++. Keep binaries portable between CPUs of the same
# architecture, including ZIPs exported from a local Blender build.
flags = ["/O2", "/fp:fast"] if sys.platform == "win32" else ["-O3", "-ffast-math", "-std=c++11"]
extensions = [
    Extension(f"molecular_blender.{name}", [f"molecular_blender/{name}.pyx"],
              include_dirs=[np.get_include()],
              define_macros=[("NPY_NO_DEPRECATED_API", "NPY_1_7_API_VERSION"),
                             ("NPY_TARGET_VERSION", "NPY_1_19_API_VERSION")],
              extra_compile_args=flags, language="c++")
    for name in ("compute_ao_cy", "marching_cube_cy", "esp_kernel")
]

# A checkout can target several Blender/Python installations. Cython generates
# different C++ from the NumPy 1.x and 2.x declarations for the same .pyx.
generated = Path("build") / f"cython-{Cython.__version__}-{sys.implementation.cache_tag}-numpy-{np.__version__}"
setup(version=addon_version(), ext_modules=cythonize(extensions, language_level=3, build_dir=str(generated)))
