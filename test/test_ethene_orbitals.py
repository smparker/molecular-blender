"""Independent reference checks for an occupied pi orbital under a 30-degree twist."""
import unittest
from unittest.mock import patch

import numpy as np

from ethene_reference import FIXTURES, LEVELS, outer_crossing_errors
from molecular_blender import compute_ao_py, orbitals
from molecular_blender.constants import bohr2ang
from molecular_blender.importers import molecule_from_file

try:
    from molecular_blender import compute_ao_cy
except ImportError:
    compute_ao_cy = None

try:
    import bpy
except ImportError:
    bpy = None


BACKENDS = [compute_ao_py] + ([compute_ao_cy] if compute_ao_cy else [])


class TestEtheneOrbital(unittest.TestCase):
    def test_signed_point_and_grid_values_against_pyscf(self):
        for stem in ('planar', 'twisted'):
            probes = np.loadtxt(FIXTURES / f'{stem}-samples.dat')
            geometry = molecule_from_file(FIXTURES / f'{stem}.molden', {})
            for backend in BACKENDS:
                with self.subTest(stem=stem, backend=backend.__name__), \
                     patch.object(orbitals, 'Shell', backend.Shell), \
                     patch.object(orbitals, 'OrbitalCalculator', backend.OrbitalCalculater):
                    model = orbitals.MOData.from_dict(geometry)
                    self.assertEqual(model.homo(), 8)
                    field = model.scalar_field('homo')
                    np.testing.assert_allclose(field.batch(probes[:, :3] * bohr2ang), probes[:, 3],
                                               atol=2.e-6, rtol=0)
                    # The production mesher samples boxes. Exercise the sorted z
                    # axis along each off-axis line, including its signed node.
                    for x, y in np.unique(probes[:, :2], axis=0):
                        line = probes[(probes[:, 0] == x) & (probes[:, 1] == y)]
                        line = line[np.argsort(line[:, 2])]
                        values = field.grid(np.array([x]) * bohr2ang, np.array([y]) * bohr2ang,
                                            line[:, 2] * bohr2ang)
                        np.testing.assert_allclose(values.ravel(), line[:, 3], atol=2.e-6, rtol=0)

    @unittest.skipUnless(bpy, 'checks published Blender meshes')
    def test_published_signed_meshes_close_at_reference_crossings(self):
        from molecular_blender.plotter import BlendMoleculeSetup, draw_surfaces
        from molecular_blender.scene_resources import remove_molecule
        from molecular_blender.scene_state import ID, read_record, owned_objects

        spacing = .12
        for stem in ('planar', 'twisted'):
            roots = np.loadtxt(FIXTURES / f'{stem}-crossings.dat')
            for backend in BACKENDS:
                with self.subTest(stem=stem, backend=backend.__name__), \
                     patch.object(orbitals, 'Shell', backend.Shell), \
                     patch.object(orbitals, 'OrbitalCalculator', backend.OrbitalCalculater):
                    model, styler, opts = BlendMoleculeSetup(
                        bpy.context, str(FIXTURES / f'{stem}.molden'), orbital='homo',
                        cumulative=False, isovalues='.05,.005', resolution=spacing, remesh=0)
                    try:
                        draw_surfaces(model, styler, bpy.context, opts)
                        owner = model.scene_owner
                        levels = {s['id']: s['isovalue'] for s in read_record(owner)['surfaces']}
                        self.assertCountEqual(levels.values(), LEVELS)
                        bounds = model.orbitals.scalar_field('homo').bounds(.005 * .001)
                        for obj in owned_objects(owner, 'surface'):
                            level = levels[obj[ID]]
                            mesh = obj.data
                            vertices = np.array([v.co[:] for v in mesh.vertices])
                            faces = np.array([p.vertices[:] for p in mesh.polygons])
                            self.assertGreater(len(faces), 0)
                            edges = np.sort(faces[:, [0, 1, 1, 2, 2, 0]].reshape(-1, 2), axis=1)
                            unique, uses = np.unique(edges, axis=0, return_counts=True)
                            self.assertTrue(np.all(uses == 2), f'Open or nonmanifold edge at {level}')
                            neighbors = [set() for _ in vertices]
                            for i, j in unique:
                                neighbors[i].add(j)
                                neighbors[j].add(i)
                            visited, pending = {0}, [0]
                            while pending:
                                new = neighbors[pending.pop()] - visited
                                visited.update(new)
                                pending.extend(new)
                            self.assertEqual(len(visited), len(vertices), f'Disconnected lobe at {level}')
                            self.assertEqual(len(vertices) - len(unique) + len(faces), 2)
                            self.assertGreater(np.min(vertices - bounds[0]), spacing)
                            self.assertGreater(np.min(bounds[1] - vertices), spacing)
                            errors = outer_crossing_errors(vertices[faces], roots, level)
                            self.assertLess(max(errors), .02, f'Outer crossing at {level}')
                    finally:
                        remove_molecule(model.scene_owner)
