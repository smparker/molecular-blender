"""Contracts shared by numerical Python callers, Blender operators and scene cleanup."""
from pathlib import Path
import tempfile
import unittest

import numpy as np

from molecular_blender import _HAS_BPY
from molecular_blender.constants import ang2bohr
from molecular_blender.diagnostics import InputError
from molecular_blender.importers import molecule_from_file
from molecular_blender.molecular_data import Molecule, VolumeData
from molecular_blender.options import process_options, validated_options
from molecular_blender.surface_calculation import compute_surface_vertex_sets

EXAMPLES = Path(__file__).resolve().parents[1] / 'examples'


class TestScientificBoundaries(unittest.TestCase):
    def test_sparse_source_ids_in_ring_and_scene_free_contour(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'triangle.molden'
            path.write_text('[Molden Format]\n[Atoms] Angs\nC 9 6 0 0 0\nC 4 6 1 0 0\nC 11 6 0 1 0\n')
            options = validated_options({'find_aromatic': True})
            model = Molecule.from_dict('triangle', molecule_from_file(path, options))
            model.determine_bonding(options)
            self.assertEqual([a.index for a in model.atoms], [4, 9, 11])
            self.assertEqual([set(r) for r in model.rings], [{0, 1, 2}])

        options = validated_options(dict(isovalues='.1', cumulative=False, resolution=.25))
        model = Molecule.from_dict('water', molecule_from_file(EXAMPLES / 'water-sto-3g.molden', options))
        result, _ = compute_surface_vertex_sets(model, process_options(options), indexed=True)
        self.assertEqual(len(result), 2)
        for level in result:
            self.assertGreater(len(level['vertices']), 3)
            self.assertEqual(level['faces'].shape[1], 3)
            self.assertLess(level['faces'].max(), len(level['vertices']))

    def test_cube_and_wavefunction_field_contract(self):
        data = np.arange(8, dtype=np.float32).reshape(2, 2, 2)
        origin, axes = np.array([1., 2., 3.]), np.diag([-1., 2., 3.])
        vol = VolumeData(origin, axes, data.shape, data, datasets=np.stack([data, -data], axis=-1),
                         labels=[7, 12], units='V', field_type='signed')
        field = vol.scalar_field()
        vol.select_dataset(1)
        self.assertEqual(field.info.selection, 'dataset 7')
        self.assertEqual((field.info.signed, field.info.units), (True, 'V'))
        np.testing.assert_allclose(field.bounds(), [[0, 2, 3], [1, 4, 6]])
        np.testing.assert_array_equal(field.grid([1, 0], [2, 4], [3, 6]), data)
        midpoint = [.5, 3., 4.5]
        self.assertAlmostEqual(field.point(midpoint), 3.5)
        self.assertAlmostEqual(vol.scalar_field().point(midpoint), -3.5)

        path = EXAMPLES / 'water-sto-3g.molden'
        model = Molecule.from_dict('water', molecule_from_file(path, {})).orbitals
        for selection, signed in [('homo', True), ('density', False), ('spin', True)]:
            field = model.scalar_field(selection)
            self.assertEqual((field.info.source, field.info.signed), (str(path), signed))
            raw = model.get_density(selection) if model.is_density(selection) else model.get_orbital(model.homo())
            xyz = np.array([[.2, .3, -.4], [.4, .5, -.7]])
            expected = [raw.value(*p) for p in xyz * ang2bohr]
            np.testing.assert_allclose(field.batch(xyz), expected, rtol=1e-5, atol=1e-8)
            grid = field.grid(*xyz.T)
            np.testing.assert_allclose([grid[0, 0, 0], grid[1, 1, 1]], expected, atol=1e-6)
            np.testing.assert_allclose(np.asarray(field.bounds()) * ang2bohr, raw.bounding_box(), rtol=1e-6)
            self.assertEqual(field.grid([], [1], [0]).shape, (0, 1, 1))
            self.assertIsNone(field.gradient)

    def test_manifest_error_keeps_both_locations(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            path = root / 'input.json'
            path.write_text('{"molecules":[{"filename":"bad.molden"}]}')
            (root / 'bad.molden').write_text('[Molden Format]\n[Atoms] Angs\nH 42 1 0 xx 0\n')
            with self.assertRaises(InputError) as raised:
                molecule_from_file(path, {})
            self.assertEqual(raised.exception.context['entry'], 1)
            cause = raised.exception.__cause__
            self.assertEqual((cause.line, cause.context['atom']), (3, 42))
            self.assertIn('bad.molden:3: [Atoms] atom 42', str(raised.exception))
            path.write_text('{\n?}')
            with self.assertRaises(InputError) as raised:
                molecule_from_file(path, {})
            self.assertEqual((raised.exception.line, raised.exception.context['column']), (2, 1))


@unittest.skipUnless(_HAS_BPY, 'requires Blender')
class TestArchitectureInBlender(unittest.TestCase):
    def test_all_exposed_options_match_registered_rna(self):
        import bpy
        from molecular_blender import MolecularBlender
        from molecular_blender.options import OPTIONS, IMPORT_PROPERTIES
        if not MolecularBlender.is_registered:
            bpy.utils.register_class(MolecularBlender)
            self.addCleanup(bpy.utils.unregister_class, MolecularBlender)
        properties = getattr(bpy.ops, 'import').molecular_blender.get_rna_type().properties
        for name in IMPORT_PROPERTIES:
            spec, rna = OPTIONS[name], properties[name]
            with self.subTest(name=name):
                if isinstance(spec.default, float):
                    self.assertAlmostEqual(rna.default, spec.ui['default'])
                else:
                    self.assertEqual(rna.default, spec.ui['default'])
                if spec.choices:
                    self.assertEqual(set(rna.enum_items.keys()), {e[0] for e in spec.ui['items']})
        self.assertEqual(OPTIONS['plot_style'].default, 'sticks')
        self.assertEqual(properties['plot_style'].default, 'fixedbs')

    def test_cleanup_preserves_shared_data_and_user_decoration(self):
        import bpy
        from molecular_blender.plotter import BlendMoleculeSetup
        from molecular_blender.scene_resources import remove_molecule

        def ids():
            return {block for key in ('objects', 'meshes', 'curves', 'materials', 'actions', 'collections')
                    for block in getattr(bpy.data, key)}

        first_ids = ids()
        self.addCleanup(lambda: bpy.data.batch_remove(ids() - first_ids))
        args = dict(recycle_materials=True, draw_bond_orders=True, ignore_hydrogen=False,
                    charges='scale', isovalues='')
        path = str(EXAMPLES / 'bond_orders.sdf')
        first, _, _ = BlendMoleculeSetup(bpy.context, path, **args)
        second, _, _ = BlendMoleculeSetup(bpy.context, path, **args)
        self.assertEqual(first.materials['c'], second.materials['c'])
        bond = bpy.data.objects[first.bonds[0].name[0]]
        profile = bond.data.bevel_object
        shared = bpy.data.objects.new('External curve user', bond.data)
        bpy.context.scene.collection.objects.link(shared)
        decoration = bpy.data.objects.new('External decoration', None)
        bpy.context.scene.collection.objects.link(decoration)
        decoration.parent = first.scene_owner
        decoration.location = (1., 2., 3.)
        bpy.context.view_layer.update()
        placement = decoration.matrix_world.translation.copy()
        remove_molecule(first.scene_owner)
        self.assertEqual(shared.data.bevel_object, profile)
        self.assertTrue(profile.use_fake_user)
        self.assertIsNone(decoration.parent)
        np.testing.assert_allclose(placement, decoration.matrix_world.translation)
        self.assertIn(second.materials['c'], bpy.data.materials)
        self.assertTrue(all(bpy.data.objects[name].parent == second.scene_owner
                            for bond in second.bonds for name in bond.name))

        mat = bpy.data.materials[second.materials['c']]
        mat.node_tree.nodes.get('Principled BSDF').inputs['Roughness'].default_value = .91
        third, _, _ = BlendMoleculeSetup(bpy.context, path, **args)
        self.assertNotEqual(second.materials['c'], third.materials['c'])
        remove_molecule(third.scene_owner)
        remove_molecule(second.scene_owner)
        self.assertIn(profile.name, bpy.data.objects)
        # With no external users, removing an animation releases all its IDs.
        before = ids()
        model, _, _ = BlendMoleculeSetup(bpy.context, str(EXAMPLES / 'benzopentyl_twist.xyz'),
                                         isovalues='', charges='scale', gradient=True)
        self.assertTrue(set(bpy.data.actions) - before)
        remove_molecule(model.scene_owner)
        self.assertEqual(ids(), before)
