"""Profile a signed HOMO from field sampling through construction of two Blender meshes.

blender -b --factory-startup --python-exit-code 1 --python benchmarks/orbital_kernels.py -- [FILE ...]
Pass one file per fresh Blender process to compare peak resident memory.
Use -- --response to sample the first supplied phenanthroline DIPL response on a 100³ grid.
"""
from pathlib import Path
from statistics import median
import sys
from time import perf_counter
from types import SimpleNamespace
from unittest.mock import patch

import bpy
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from molecular_blender import compute_ao, isosurfaces, marching_cube
from molecular_blender.importers import molecule_from_file
from molecular_blender.orbitals import MOData
from molecular_blender.scene_meshes import create_surface_mesh_data
from molecular_blender.surface_geometry import index_field_triangles


def profile(path):
    model = MOData.from_dict(molecule_from_file(path, {}))
    field = model.scalar_field(model.homo())
    spacing = .06
    samples = []
    for run in range(4):
        times = dict(ao=0., marching=0., indexing=0., mesh=0.)

        def timed(stage, function):
            def measure(*args, **kwargs):
                start = perf_counter()
                value = function(*args, **kwargs)
                times[stage] += perf_counter() - start
                return value
            return measure

        evaluator = SimpleNamespace(box_values=timed('ao', field.evaluator.box_values),
                                    bounding_box=field.evaluator.bounding_box)
        with patch.object(field, 'evaluator', evaluator), \
             patch.object(isosurfaces, 'marching_cube_box', timed('marching', isosurfaces.marching_cube_box)):
            start = perf_counter()
            contours = isosurfaces.field_isosurface(field, [.05, -.05], spacing)
            sizes = []
            for surface in contours:
                vertices, faces = timed('indexing', index_field_triangles)(surface['triangles'], spacing)
                mesh = timed('mesh', create_surface_mesh_data)('HOMO', vertices, faces, None)
                sizes.append((len(mesh.vertices), len(mesh.polygons)))
                bpy.data.meshes.remove(mesh)
            times['total'] = perf_counter() - start
        if run:
            samples.append(times)
    print(f'{path.name} HOMO {model.homo()}, ±0.05 at {spacing} Å, (vertices, faces) {sizes}')
    for stage in times:
        print(f'  {stage:10} {1000 * median(s[stage] for s in samples):8.1f} ms')


def print_peak_memory():
    if sys.platform != 'win32':
        from resource import RUSAGE_SELF, getrusage
        scale = 1024 ** (2 if sys.platform == 'darwin' else 1)
        print(f'  Peak process RSS {getrusage(RUSAGE_SELF).ru_maxrss / scale:.1f} MiB')


def profile_response(path):
    model = MOData.from_dict(molecule_from_file(path, {}))
    axis = np.linspace(-30., 30., 100, dtype=np.float32)
    start = perf_counter()
    field = model.get_density(model.response_field(0))
    prepared = perf_counter()
    grid = field.box_values(axis, axis, axis)
    sampled = perf_counter()
    print(f'{path.name}, first DIPL response, 100³ points from −30 to +30 bohr')
    print(f'  Setup {1000 * (prepared - start):.1f} ms, sampling {1000 * (sampled - prepared):.1f} ms, '
          f'range [{grid.min():.6g}, {grid.max():.6g}]')


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--response', action='store_true', help='sample a DIPL response grid')
    parser.add_argument('files', nargs='*', type=Path)
    args = parser.parse_args(sys.argv[sys.argv.index('--') + 1:] if '--' in sys.argv else [])
    if not (compute_ao.CYTHON_ENABLED and marching_cube.CYTHON_ENABLED):
        parser.error('This benchmark requires both compiled backends')
    if args.response:
        for source in args.files or [ROOT / 'examples/phenanth_polar.json']:
            profile_response(source)
            print_peak_memory()
    else:
        for source in args.files or [ROOT / 'test/fixtures/pyscf_ethene/planar.molden',
                                     ROOT / 'examples/phenanthroline.molden']:
            profile(source)
            print_peak_memory()
