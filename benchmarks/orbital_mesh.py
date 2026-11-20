"""Time construction of both HOMO lobe meshes from the same precomputed triangles.

blender -b --factory-startup --python benchmarks/orbital_mesh.py -- molecule.molden
Unix peak RSS includes Blender startup and the input arrays. Run each version in a fresh process.
"""
from pathlib import Path
import sys
from statistics import median
from time import perf_counter

import bpy

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from molecular_blender.importers import molecule_from_molden
from molecular_blender.isosurfaces import molden_isosurface
from molecular_blender.orbitals import MOData
from molecular_blender.plotter import create_surface_mesh_data, index_field_triangles


def main(path):
    model = MOData.from_dict(molecule_from_molden(path, {}))
    orbital = model.get_orbital(model.homo())
    levels = orbital.isovalue_containing_proportion([.5, -.5])
    spacing = .035
    surfaces = molden_isosurface(orbital, levels, spacing)
    times = []
    for _ in range(3):
        meshes = []
        try:
            start = perf_counter()
            for surface in surfaces:
                vertices, faces = index_field_triangles(surface['triangles'], spacing)
                meshes.append(create_surface_mesh_data('HOMO', vertices, faces, weld_distance=None))
            times.append(perf_counter() - start)
        finally:
            for mesh in meshes:
                bpy.data.meshes.remove(mesh)
    print(f'Median mesh build for both lobes: {median(times):.3f} s')
    if sys.platform != 'win32':
        from resource import getrusage, RUSAGE_SELF
        scale = 1024 ** (2 if sys.platform == 'darwin' else 1)
        print(f'Peak process RSS including sampling: {getrusage(RUSAGE_SELF).ru_maxrss / scale:.1f} MiB')


if __name__ == '__main__':
    main(Path(sys.argv[sys.argv.index('--') + 1]))
