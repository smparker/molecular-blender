"""
Setup file for building wheels with cibuildwheel
Modified from the original setup.py for molecular_blender
This file is specifically for wheel building and CI/CD purposes.
For local development, continue using setup.py
"""

import os
import sys
import numpy as np
from setuptools import setup, Extension
from Cython.Build import cythonize

# Get version from __init__.py
def get_version():
    init_file = os.path.join(os.path.dirname(__file__), "molecular_blender", "__init__.py")
    with open(init_file, "r") as f:
        for line in f:
            if line.startswith("bl_info"):
                # Parse the bl_info dictionary to get version
                for i in range(20):  # Look ahead a few lines
                    line = f.readline()
                    if '"version"' in line or "'version'" in line:
                        # Extract version tuple like (1, 2, 3)
                        import re
                        match = re.search(r'\((\d+),\s*(\d+),\s*(\d+)\)', line)
                        if match:
                            return f"{match.group(1)}.{match.group(2)}.{match.group(3)}"
    return "0.0.0"

# Compiler flags - portable version (no -march=native for wheels)
extra_compile_args = ["-O3", "-ffast-math"]
if sys.platform != "win32":
    extra_compile_args.append("-std=c++11")

# Define Cython extensions
extensions = [
    Extension(
        "molecular_blender.compute_ao_cy",
        sources=["molecular_blender/compute_ao_cy.pyx"],
        include_dirs=[np.get_include()],
        define_macros=[("NPY_NO_DEPRECATED_API", "NPY_1_7_API_VERSION")],
        extra_compile_args=extra_compile_args,
        language="c++",
    ),
    Extension(
        "molecular_blender.marching_cube_cy",
        sources=["molecular_blender/marching_cube_cy.pyx"],
        include_dirs=[np.get_include()],
        define_macros=[("NPY_NO_DEPRECATED_API", "NPY_1_7_API_VERSION")],
        extra_compile_args=extra_compile_args,
        language="c++",
    ),
]

# Cythonize extensions
extensions = cythonize(
    extensions,
    compiler_directives={
        "language_level": "3",
        "cdivision": True,
        "profile": True,
    },
    include_path=[np.get_include()],
)

setup(
    name="molecular_blender",
    version=get_version(),
    description="Blender addon for molecular visualization",
    long_description=open("README.md").read() if os.path.exists("README.md") else "",
    long_description_content_type="text/markdown",
    author="Shane Parker",
    url="https://github.com/smparker/molecular-blender",
    packages=["molecular_blender"],
    ext_modules=extensions,
    python_requires=">=3.9",
    install_requires=[
        "numpy>=1.19.0",
    ],
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Python :: 3.12",
        "Programming Language :: Python :: 3.13",
        "Programming Language :: Cython",
        "Topic :: Scientific/Engineering :: Chemistry",
        "Topic :: Scientific/Engineering :: Visualization",
    ],
    zip_safe=False,
)

