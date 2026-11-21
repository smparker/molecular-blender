from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import numpy as np

from molecular_blender import _HAS_BPY
from molecular_blender import wavefunction_esp as esp
from molecular_blender.esp_jobs import ESPRequest
from molecular_blender.importers import molecule_from_molden
from molecular_blender.molecular_data import Molecule
from molecular_blender.orbitals import orbital_metadata
from molecular_blender.surface_jobs import SavedSurfaceSource

ROOT = Path(__file__).resolve().parents[1]
REF = ROOT / 'test/fixtures/wavefunction_esp'
WATER = ROOT / 'examples/water_esp/water.molden'


class TestWavefunctionESP(unittest.TestCase):
    def test_independent_reference_points_and_isodensities(self):
        for name, path, charge in [('water', WATER, 0), ('cartesian', REF.parent/'pyscf_water/cartesian.molden', 0),
                                   ('orca', REF.parent/'orca_water_fg/water.molden', 0), ('heh', REF/'heh.molden', 1)]:
            with self.subTest(name=name):
                field = esp.WavefunctionESP(molecule_from_molden(path, {}), charge)
                samples = np.loadtxt(REF / (name+'.dat'))
                points = samples[:, :3]
                for kernel in dict.fromkeys((esp.coulomb_sum, None)):
                    with patch.object(esp, 'coulomb_sum', kernel):
                        values = field.potential_bohr(points)
                    nuclear = (field.charges / np.linalg.norm(points[:, None]-field.nuclei, axis=2)).sum(axis=1)
                    np.testing.assert_allclose(nuclear, samples[:, 3], atol=1e-6)
                    np.testing.assert_allclose(nuclear-values, samples[:, 4], atol=5e-7)
                    np.testing.assert_allclose(values, samples[:, 5], atol=5e-7)
                    self.assertAlmostEqual(values[2]*1e4, charge, delta=.0003)
                with self.assertRaisesRegex(ValueError, 'singular'):
                    field.potential_bohr(field.nuclei[:1])

    @unittest.skipUnless(esp.coulomb_sum, 'requires compiled ESP kernel')
    def test_native_workspace_rejects_uninitialized_dependencies_and_invalid_points(self):
        points = np.zeros((1, 3))
        exponents = np.ones(1)
        _, rules, counts = esp._indices(1)
        args = [points, 1, exponents, points, np.ones((1, 4)), rules.copy(), counts, esp._boys_table()]
        args[5][2, 1] = 1  # In bounds, but unavailable in the first recurrence step.
        with self.assertRaisesRegex(ValueError, 'indices'):
            esp.coulomb_sum(*args)
        args[5] = rules
        for point, exponent in [(np.inf, 1.), (np.nan, 1.), (0., -1.)]:
            with self.subTest(point=point, exponent=exponent):
                args[0] = np.full((1, 3), point)
                args[2] = np.array([exponent])
                with self.assertRaisesRegex(ValueError, 'finite'):
                    esp.coulomb_sum(*args)

    def test_overlapping_occupied_orbitals_in_one_spin_channel(self):
        data = molecule_from_molden(REF/'heh.molden', {})
        orbital = dict(data['mo'][0], occup='1')
        data['mo'] = [orbital, orbital.copy()]
        with self.assertRaisesRegex(ValueError, 'mutually orthogonal'):
            esp.WavefunctionESP(data, 1)
        data['mo'][1]['spin'] = 'Beta'
        esp.WavefunctionESP(data, 1)

    def test_input_contract_and_occupied_only_file(self):
        data = molecule_from_molden(REF/'heh.molden', {})
        self.assertEqual([a['nuclear_charge'] for a in data['atoms']], [2, 1])
        data['mo'] = [mo for mo in data['mo'] if float(mo['occup'])]
        esp.WavefunctionESP(data, 1)
        with self.assertRaisesRegex(ValueError, 'requires 3'):
            esp.WavefunctionESP(data, 0)
        with self.assertRaisesRegex(ValueError, 'core or ECP'):
            esp.WavefunctionESP(dict(data, has_ecp=True), 1)
        for z in [0, 1, None]:
            with self.subTest(z=z), self.assertRaisesRegex(ValueError, 'full nuclear charges'):
                esp.WavefunctionESP(dict(data, atoms=[dict(data['atoms'][0], nuclear_charge=z), data['atoms'][1]]), 1)
        data['mo'][0]['coeff'][0] *= .5
        with self.assertRaisesRegex(ValueError, 'Integrated density'):
            esp.WavefunctionESP(data, 1)

    def test_cached_points_units_and_source_invalidation(self):
        from molecular_blender.constants import bohr2ang
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp)/'ion.molden'
            path.write_bytes((REF/'heh.molden').read_bytes())
            data = molecule_from_molden(path, {})
            model = Molecule.from_dict('ion', data)
            saved = SavedSurfaceSource(str(path), [(a.index, a.element.symbol, tuple(a.position)) for a in model.atoms],
                                       orbital_metadata(model.orbitals))
            reference = np.loadtxt(REF/'heh.dat')
            req = ESPRequest(saved, 'auto', 1, 'final', 'V', reference[:, :3]*bohr2ang, Path(tmp))
            np.testing.assert_allclose(req.calculate({}), reference[:, -1]*esp.HARTREE_IN_VOLTS, atol=4e-6)
            with patch.object(esp.WavefunctionESP, 'batch', side_effect=AssertionError('cache miss')):
                req.units = 'Hartree/e'
                np.testing.assert_allclose(req.calculate({}), reference[:, -1], atol=5e-7)
            with patch.object(SavedSurfaceSource, 'read', side_effect=AssertionError('density cache miss')):
                req.quality = 'preview'
                req.calculate({})
            path.write_text(path.read_text()+'\n[Core]\n1 : 2\n')
            with self.assertRaisesRegex(ValueError, 'core or ECP'):
                req.calculate({})


@unittest.skipUnless(_HAS_BPY, 'requires Blender')
class TestWavefunctionMap(unittest.TestCase):
    def test_preview_to_final_remap_palette_and_legend(self):
        import bpy
        from molecular_blender import scene_panel
        from molecular_blender.field_legend import write_field_legend
        from molecular_blender.plotter import BlendMoleculeSetup, draw_surfaces
        from molecular_blender.scene_esp import BASE_MESH, clear_esp, mesh_points, remap_contour, scale_nodes
        from molecular_blender.scene_state import read_record
        from molecular_blender.scene_surfaces import surface_groups
        scene_panel.register()
        self.addCleanup(scene_panel.unregister)
        import xml.etree.ElementTree as ET
        ids = set(bpy.data.objects) | set(bpy.data.collections)
        self.addCleanup(lambda: bpy.data.batch_remove((set(bpy.data.objects) | set(bpy.data.collections)) - ids))
        model, styler, options = BlendMoleculeSetup(bpy.context, str(WATER), volume='density', remesh=4,
                                                  resolution=.23, cumulative=False, isovalues='.01')
        draw_surfaces(model, styler, bpy.context, options)
        owner = model.scene_owner
        record = read_record(owner)
        key, group = next(iter(surface_groups(owner, record).items()))
        obj, = group['objects']
        bpy.context.view_layer.objects.active = owner
        op = bpy.ops.molecular_blender.compute_esp
        self.assertEqual(op(group=key, quality='preview', palette='viridis'), {'FINISHED'})
        self.assertEqual(obj.active_material['mb_esp_quality'], 'preview')
        self.assertEqual(op(group=key), {'FINISHED'})
        self.assertEqual(obj.active_material['mb_esp_quality'], 'final')
        self.assertEqual(len(scale_nodes(obj.active_material)[0].color_ramp.elements), 32)
        field = esp.WavefunctionESP(molecule_from_molden(WATER, {}))
        points = mesh_points(obj.data)
        expected = field.batch(points[::19])
        np.testing.assert_allclose([v.value for v in obj.data.attributes['mb_esp'].data][::19], expected, atol=1e-6)
        with tempfile.TemporaryDirectory() as tmp:
            legend = Path(tmp)/'ramp.svg'
            for palette, start, end in [('viridis', '#440154', '#fde725'), ('magma', '#000004', '#fcfdbf')]:
                self.assertEqual(bpy.ops.molecular_blender.edit_field_scale(group=key, palette=palette), {'FINISHED'})
                write_field_legend(obj.active_material, legend)
                stops = ET.parse(legend).findall('.//{*}stop')
                self.assertEqual((stops[0].get('stop-color'), stops[-1].get('stop-color')), (start, end))
        remapped = remap_contour(bpy.context, record, obj, obj[BASE_MESH], {})
        try:
            self.assertIs(remapped.materials[0], obj.active_material)
            np.testing.assert_allclose([v.value for v in remapped.attributes['mb_esp'].data][::19], expected, atol=1e-6)
        finally:
            bpy.data.meshes.remove(remapped)
        from molecular_blender.esp_jobs import ContourESPRequests
        from molecular_blender.scene_surfaces import _contour_requests, replace_surface_geometry
        from molecular_blender.surface_jobs import SurfaceJob
        from molecular_blender.surface_calculation import compute_surface_vertex_sets
        options['isovalues'] = [.015]
        definitions, _ = compute_surface_vertex_sets(model, options, indexed=True)
        requests = _contour_requests(bpy.context, record, group, definitions, options['resolution'])
        old = obj.data
        cancelled = SurfaceJob(bpy.app.binary_path, ContourESPRequests(requests), {})
        cancelled.cancel()
        self.assertEqual(obj.data, old)
        job = SurfaceJob(bpy.app.binary_path, ContourESPRequests(requests), {})
        try:
            job.process.wait(timeout=60)
            results = job.result()
            ready = {level: (req, results[level]) for level, req in requests.items()}
            replace_surface_geometry(bpy.context, owner, record, key, group, definitions, options, styler, ready)
        finally:
            job.close()
        self.assertNotEqual(obj.data, old)
        expected = field.batch(mesh_points(obj.data)[::19])
        np.testing.assert_allclose([v.value for v in obj.data.attributes['mb_esp'].data][::19], expected, atol=1e-6)
        base = obj[BASE_MESH]
        clear_esp(group)
        self.assertEqual(obj.data, base)
