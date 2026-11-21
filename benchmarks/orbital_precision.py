"""Compare controlled AO and marching-cubes builds for two ethene pi orbitals.

Run with python3 benchmarks/orbital_precision.py from a checkout with Cython,
setuptools, NumPy, SciPy and a C++ compiler accepting Clang/GCC math flags.
Only temporary copies of the native sources are built. PySCF is used offline
to generate the checked-in reference values and continuous contour crossings.
"""
import importlib.util
from functools import partial
from pathlib import Path
import re
import sys
import tempfile
from unittest.mock import patch

from Cython.Build import cythonize
import numpy as np
from scipy.spatial import cKDTree
from setuptools import Distribution, Extension

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT), str(ROOT / 'test')]

from ethene_reference import FIXTURES, LEVELS, outer_crossing_errors
from molecular_blender import isosurfaces, orbitals
from molecular_blender.constants import ang2bohr
from molecular_blender.importers import molecule_from_file


def build_kernels(folder, label, dtype, math_flag):
    """Use the same sources, headers, compiler and optimization level for every build."""
    package = ROOT / 'molecular_blender'
    extensions = []
    for name in ('compute_ao_cy', 'marching_cube_cy'):
        source = (package / (name + '.pyx')).read_text()
        if dtype == np.float64:
            # Widen the Cython C scalars, memoryviews, struct fields and NumPy
            # buffers together. Leave existing double-precision operations alone.
            source = re.sub(r'\bfloat\b', 'double', source).replace('np.float32', 'np.float64')
        path = folder / f'{name}_{label}.pyx'
        path.write_text(source)
        extensions.append(Extension('molecular_blender.' + path.stem, [str(path)],
                                    include_dirs=[np.get_include()], language='c++',
                                    extra_compile_args=['-O3', math_flag, '-std=c++11']))
    dist = Distribution({'ext_modules': cythonize(extensions, quiet=True, language_level=3)})
    build = dist.get_command_obj('build_ext')
    build.build_temp = str(folder / 'objects')
    build.build_lib = str(folder / 'lib')
    dist.run_command('build_ext')
    modules = []
    for ext in extensions:
        spec = importlib.util.spec_from_file_location(ext.name, build.get_ext_fullpath(ext.name))
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        modules.append(module)
    return modules


def orbital(geometry, ao, threshold=1.e-6):
    with patch.object(orbitals, 'Shell', partial(ao.Shell, thr=threshold)), \
         patch.object(orbitals, 'OrbitalCalculator', ao.OrbitalCalculater):
        model = orbitals.MOData.from_dict(geometry)
        return model.get_orbital(model.homo())


def contours(evaluator, marcher, dtype, spacing):
    # Lock the sampling box, tile boundaries, grid coordinates and imported
    # coefficients across builds. Blender output positions use float32 in all cases.
    lower = np.array([-5.13, -4.31, -4.37])
    upper = lower + [10.4, 8.8, 8.8]
    counts = np.rint((upper - lower) / spacing).astype(int)
    def grid(x, y, z):
        axes = [np.ascontiguousarray(axis * ang2bohr, dtype=dtype) for axis in (x, y, z)]
        return evaluator.box_values(*axes)
    with patch.object(isosurfaces, 'marching_cube_box', marcher.marching_cube_box):
        surfaces = isosurfaces.isosurface_tiled(lower, upper, counts, LEVELS, grid, np.eye(3), 'HOMO')
    shapes = [np.asarray(s['triangles']) for s in surfaces]
    for triangles in shapes:
        vertices = triangles.reshape(-1, 3)
        if np.min(vertices - lower) < spacing or np.min(upper - vertices) < spacing:
            raise ValueError('Comparison box has insufficient clearance around a lobe')
    return shapes


def compare_meshes(left, right):
    """Symmetric maximum nearest-vertex distance, in angstroms."""
    a = np.unique(left.reshape(-1, 3), axis=0)
    b = np.unique(right.reshape(-1, 3), axis=0)
    return max(cKDTree(a).query(b)[0].max(), cKDTree(b).query(a)[0].max())


def main():
    variants = [('fast32', np.float32, '-ffast-math'), ('strict32', np.float32, '-fno-fast-math'),
                ('strict64', np.float64, '-fno-fast-math')]
    with tempfile.TemporaryDirectory(prefix='mb-orbital-precision-') as tmp:
        kernels = {name: build_kernels(Path(tmp), name, dtype, flag) for name, dtype, flag in variants}
        for stem in ('planar', 'twisted'):
            geometry = molecule_from_file(FIXTURES / f'{stem}.molden', {})
            probes = np.loadtxt(FIXTURES / f'{stem}-samples.dat')
            roots = np.loadtxt(FIXTURES / f'{stem}-crossings.dat')
            evaluators = {}
            for name, _, _ in variants:
                ao, _ = kernels[name]
                evaluators[name] = orbital(geometry, ao)
                for threshold, field in [(1.e-6, evaluators[name]), (1.e-12, orbital(geometry, ao, 1.e-12))]:
                    values = np.array([field.value(*p) for p in probes[:, :3]])
                    print(f'{stem} {name} screening {threshold:g}: point max '
                          f'{max(abs(values - probes[:, 3])):.3g}', flush=True)
            for spacing in (.10, .05):
                shapes = {}
                for name, dtype, _ in variants:
                    shapes[name] = contours(evaluators[name], kernels[name][1], dtype, spacing)
                    errors = [max(outer_crossing_errors(triangles, roots, level))
                              for triangles, level in zip(shapes[name], LEVELS)]
                    print(f'{stem} {spacing:.2f} {name}: outer crossing max {max(errors):.6g} A', flush=True)
                for left, right in [('fast32', 'strict32'), ('strict32', 'strict64')]:
                    for level, a, b in zip(LEVELS, shapes[left], shapes[right]):
                        print(f'  {left}/{right} {level:+.3f}: {len(a)}/{len(b)} triangles, '
                              f'max vertex displacement {compare_meshes(a, b):.6g} A', flush=True)


if __name__ == '__main__':
    main()
