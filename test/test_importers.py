"""Tests for molecular_blender.importers and molecular_blender.turbomole_reader"""

import json
import os
from pathlib import Path
import tempfile
import unittest
import numpy as np

from molecular_blender.importers import (
    Reader, molecule_from_file, molecule_from_xyz, molecule_from_cube,
    molecule_from_molden, molecule_from_json,
    make_atom_dict, make_snap_dict,
)
from molecular_blender.turbomole_reader import (
    DiplParser, read_dipl_file, read_turbomole_polarizability,
)

EXAMPLES = os.path.join(os.path.dirname(__file__), "../examples")


# -- helper constructors --

class TestMakeAtomDict(unittest.TestCase):
    def test_basic(self):
        d = make_atom_dict("C", [1.0, 2.0, 3.0], 0, 0.5, [0.1, 0.2, 0.3], False)
        assert d["symbol"] == "c"
        assert d["position"] == [1.0, 2.0, 3.0]
        assert d["index"] == 0
        assert d["charge"] == 0.5
        assert d["gradient"] == [0.1, 0.2, 0.3]
        assert d["hidden"] is False
        assert d["trajectory"] == []

    def test_hidden(self):
        d = make_atom_dict("H", [0, 0, 0], 1, 0.0, [0, 0, 0], True)
        assert d["hidden"] is True


class TestMakeSnapDict(unittest.TestCase):
    def test_basic(self):
        d = make_snap_dict([1.0, 2.0, 3.0], 0.5, [0.1, 0.2, 0.3])
        assert d["position"] == [1.0, 2.0, 3.0]
        assert d["charge"] == 0.5
        assert d["gradient"] == [0.1, 0.2, 0.3]


# -- Reader class --

class TestReader(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp()

    def tearDown(self):
        import shutil
        shutil.rmtree(self.tmp, ignore_errors=True)

    def _write(self, name, content):
        path = os.path.join(self.tmp, name)
        with open(path, "w") as f:
            f.write(content)
        return path

    def test_context_manager(self):
        f = self._write("test.txt", "line1\nline2\n")
        with Reader(f) as r:
            assert r.readline().strip() == "line1"
            assert r.readline().strip() == "line2"

    def test_eof_raises(self):
        f = self._write("test.txt", "only\n")
        with Reader(f) as r:
            r.readline()
            with self.assertRaises(EOFError):
                r.readline()

    def test_marks(self):
        f = self._write("test.txt", "aaa\nbbb\nccc\n")
        with Reader(f) as r:
            line = r.readline()
            r.set_mark("first", line)
            r.readline()
            r.readline()
            assert r.is_marked("first")
            assert not r.is_marked("nonexistent")
            restored = r.restore_mark("first")
            assert restored.strip() == "aaa"


# -- molecule_from_file dispatch --

class TestMoleculeFromFileDispatch(unittest.TestCase):
    def test_xyz(self):
        geo = molecule_from_file(os.path.join(EXAMPLES, "tetracene_dimer.xyz"), {})
        assert len(geo["atoms"]) == 60

    def test_cube(self):
        geo = molecule_from_file(os.path.join(EXAMPLES, "water.cube"), {})
        assert len(geo["atoms"]) == 3
        assert "volume" in geo

    def test_molden(self):
        geo = molecule_from_file(os.path.join(EXAMPLES, "water-sto-3g.molden"), {})
        assert len(geo["atoms"]) > 0
        assert "basis" in geo

    def test_json(self):
        geo = molecule_from_file(os.path.join(EXAMPLES, "phenanth_polar.json"), {})
        assert len(geo["atoms"]) == 22
        assert "polar" in geo

    def test_unsupported_format(self):
        for filename in ("something.pdb", "something.mol2", "extensionless"):
            with self.subTest(filename=filename):
                with self.assertRaises(ValueError) as error:
                    molecule_from_file(filename, {})
                self.assertIn(filename, str(error.exception))
                self.assertIn("Supported formats: XYZ, SDF/MOL/SD", str(error.exception))


# -- XYZ reading --

class TestXYZRead(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp()

    def tearDown(self):
        import shutil
        shutil.rmtree(self.tmp, ignore_errors=True)

    def _write(self, name, content):
        path = os.path.join(self.tmp, name)
        with open(path, "w") as f:
            f.write(content)
        return path

    def test_basic_read(self):
        geo = molecule_from_xyz(os.path.join(EXAMPLES, "tetracene_dimer.xyz"), {})
        assert len(geo["atoms"]) == 60
        assert geo["atoms"][0]["symbol"] == "c"

    def test_read_with_charges(self):
        geo = molecule_from_xyz(
            os.path.join(EXAMPLES, "tio2_excited_state.xyz"),
            {"plot_type": "animate", "ignore_hydrogen": False}
        )
        atoms = geo["atoms"]
        assert len(atoms) > 0
        has_charge = any(a["charge"] != 0.0 for a in atoms)
        assert has_charge

    def test_animation(self):
        opts = {"plot_type": "animate"}
        geo = molecule_from_xyz(os.path.join(EXAMPLES, "benzopentyl_twist.xyz"), opts)
        atoms = geo["atoms"]
        assert len(atoms[0]["trajectory"]) > 1
        assert opts["plot_type"] == "animate"

    def test_ignore_hydrogen(self):
        geo = molecule_from_xyz(
            os.path.join(EXAMPLES, "tetracene_dimer.xyz"),
            {"ignore_hydrogen": True}
        )
        hidden_h = [a for a in geo["atoms"] if a["hidden"] and a["symbol"] == "h"]
        all_h = [a for a in geo["atoms"] if a["symbol"] == "h"]
        assert len(hidden_h) == len(all_h)

    def test_sample_source_positions_and_preserve_timing(self):
        path = self._write('frames.xyz', ''.join(f'1\nframe\nHe {i} 0 0\n' for i in range(7)))
        data = molecule_from_xyz(path, dict(plot_type='auto', frame_start=1, frame_end=5, frame_step=2))
        self.assertEqual(data['frame_numbers'], [1, 3, 5])
        atom, = data['atoms']
        self.assertEqual(atom['position'][0], 1)
        self.assertEqual([s['position'][0] for s in atom['trajectory']], [1, 3, 5])
        with self.assertRaisesRegex(ValueError, 'no XYZ'):
            molecule_from_xyz(path, {'frame_start': 7})

    def test_mismatched_atoms_raises(self):
        f = self._write("bad.xyz", "2\nframe1\nC 0 0 0\nH 1 0 0\n2\nframe2\nC 0 0 0\nC 1 0 0\n")
        with self.assertRaisesRegex(Exception, "order of the atoms"):
            molecule_from_xyz(f, {"plot_type": "animate"})

    def test_mismatched_natoms_raises(self):
        f = self._write("bad.xyz", "2\nframe1\nC 0 0 0\nH 1 0 0\n3\nframe2\nC 0 0 0\nH 1 0 0\nO 2 0 0\n")
        with self.assertRaisesRegex(Exception, "same number of atoms"):
            molecule_from_xyz(f, {"plot_type": "animate"})


# -- Cube reading --

class TestCubeRead(unittest.TestCase):
    def test_read_and_density(self):
        geo = molecule_from_cube(os.path.join(EXAMPLES, "water.cube"), {})
        assert len(geo["atoms"]) == 3
        vol = geo["volume"]
        dxyz = [np.linalg.norm(vol["axes"][:, i]) for i in range(3)]
        dvol = np.prod(dxyz)
        total = np.sum(vol["data"] ** 2) * dvol
        assert abs(total - 0.141058531145) < 1e-3

    def test_cube_has_volume_data(self):
        geo = molecule_from_cube(os.path.join(EXAMPLES, "water.cube"), {})
        vol = geo["volume"]
        assert "origin" in vol
        assert "axes" in vol
        assert "nres" in vol
        assert "data" in vol
        assert vol["data"].ndim == 3


# -- Molden reading --

class TestMoldenRead(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp()

    def tearDown(self):
        import shutil
        shutil.rmtree(self.tmp, ignore_errors=True)

    def _write(self, name, content):
        path = os.path.join(self.tmp, name)
        with open(path, "w") as f:
            f.write(content)
        return path

    def test_water_sto3g(self):
        geo = molecule_from_molden(os.path.join(EXAMPLES, "water-sto-3g.molden"), {})
        assert len(geo["atoms"]) == 3
        assert "basis" in geo
        assert "mo" in geo

    def test_fake_atom(self):
        geo = molecule_from_molden(os.path.join(EXAMPLES, "fake-atom.molden"), {})
        assert len(geo["atoms"]) == 1
        assert "basis" in geo
        assert "mo" in geo

    def test_phenanthroline(self):
        geo = molecule_from_molden(
            os.path.join(EXAMPLES, "phenanthroline_pbe0_ris.molden"), {}
        )
        assert len(geo["atoms"]) == 22
        elements = {a["symbol"] for a in geo["atoms"]}
        assert elements == {"c", "h", "n"}
        assert "basis" in geo
        assert "mo" in geo

    def test_missing_molden_format(self):
        f = self._write("bad.molden", "[Atoms]\nC 1 6 0.0 0.0 0.0\n")
        with self.assertRaisesRegex(Exception, "Molden Format"):
            molecule_from_molden(f, {})

    def test_missing_atoms(self):
        f = self._write("bad.molden", "[Molden Format]\n[GTO]\n")
        with self.assertRaisesRegex(Exception, "Atoms"):
            molecule_from_molden(f, {})

    def test_orca_fix(self):
        geo = molecule_from_molden(
            os.path.join(EXAMPLES, "water-sto-3g.molden"),
            {"orca_fix": True}
        )
        assert geo["denorm"] is False

    def test_sp_shells_atom_references_and_fortran_numbers(self):
        text = ('[Molden Format]\n[Atoms] Angs\n'
                'H 2 1 1D+00 0 0\n\nH 1 1 0 0 0\n[GTO]\n'
                '2 0\ns 1 1.0\n2d+00 1D+00\n'
                '1 0\nsp 1 1.0\n2D-01 0.5D+00 0.75d+00\n'
                '[MO]\nEne= -1D+00\nSpin= Alpha\nOccup= 1d+00\n'
                '1 1D+00\n2 0\n3 0\n4 0\n5 0')
        geo = molecule_from_molden(self._write('sp.molden', text), {})
        self.assertEqual([a['index'] for a in geo['atoms']], [1, 2])
        self.assertEqual([s['shell'] for s in geo['basis'][0]], ['s', 'p'])
        self.assertEqual(geo['basis'][0][0]['contractions'], [0.5])
        self.assertEqual(geo['basis'][0][1]['contractions'], [0.75])
        self.assertEqual(geo['basis'][1][0]['exponents'], [2.0])
        self.assertEqual(float(geo['mo'][0]['ene']), -1.0)
        from molecular_blender.orbitals import MOData
        orbital = MOData.from_dict(geo).get_orbital(1)
        self.assertGreater(orbital.value(0., 0., 0.), 0.)

    def test_standalone_7f_and_section_termination(self):
        text = ('[Molden Format]\n[Atoms] AU\nH 1 1 0 0 0\n'
                '[GTO]\n1 0\nf 1 1.0\n1 1\n[7F]\n[MO]\n'
                'Ene= -1\nSpin= Alpha\nOccup= 1\n'
                + ''.join(f'{i} {1 if i == 1 else 0}\n' for i in range(1, 8))
                + '[FREQ]\n123.4\n')
        geo = molecule_from_molden(self._write('7f.molden', text), {})
        self.assertEqual(len(geo['mo']), 1)
        self.assertEqual(len(geo['mo'][0]['coeff']), 10)

    def test_malformed_basis_and_coefficients(self):
        header = '[Molden Format]\n[Atoms] Angs\nH 1 1 0 0 0\n[GTO]\n'
        for basis in ('2 0\ns 1\n1 1\n', '1 0\nsp 1\n1 1\n',
                      '1 0\ns 0\n', '1 0\ns 2\n1 1\n',
                      '1 0\ns 1\n-1 1\n', '1 0\ns 1\nnan 1\n'):
            with self.subTest(basis=basis), self.assertRaises(ValueError):
                molecule_from_molden(self._write('bad.molden', header + basis), {})
        basis = '1 0\np 1\n1 1\n[MO]\nSpin= Alpha\nOccup= 1\n'
        for coefficients in ('1 1\n2 0\n', '1 1\n1 0\n3 0\n',
                             '1 1\n2 0\n3 0\n4 0\n', '1 nan\n2 0\n3 0\n'):
            with self.subTest(coefficients=coefficients), self.assertRaises(ValueError):
                molecule_from_molden(self._write('bad.molden', header + basis + coefficients), {})


# -- JSON reading --

class TestJsonRead(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp()

    def tearDown(self):
        import shutil
        shutil.rmtree(self.tmp, ignore_errors=True)

    def _write(self, name, content):
        path = os.path.join(self.tmp, name)
        with open(path, "w") as f:
            f.write(content)
        return path

    def test_phenanth_polar(self):
        geo = molecule_from_json(os.path.join(EXAMPLES, "phenanth_polar.json"), {})
        assert len(geo["atoms"]) == 22
        assert "polar" in geo
        assert "basis" in geo
        assert "mo" in geo
        assert geo["atoms"][0]["symbol"] == "c"

    def test_returns_none_for_missing_molecules(self):
        f = self._write("empty.json", '{"other": "stuff"}')
        result = molecule_from_json(f, {})
        assert result is None
        with self.assertRaisesRegex(ValueError, 'Missing molecules list'):
            molecule_from_file(f, {})

    def manifest(self, entries):
        return self._write('series.json', json.dumps({'molecules': entries}))

    def test_frame_indices_and_references_report_the_entry(self):
        self._write('a.xyz', '1\natom\nHe 0 0 0\n')
        first = {'filename': 'a.xyz', 'frame': 2}
        for invalid in (True, -1, 1.5, 1, 2):
            with self.subTest(invalid=invalid), self.assertRaisesRegex(ValueError, 'entry 2.*Frame indices'):
                molecule_from_json(self.manifest([first, dict(first, frame=invalid)]), {})
        for entry in ({'filename': 'missing.xyz'}, {'filename': 'a.xyz', 'polar': 'missing-polar'}, {}):
            with self.assertRaisesRegex(ValueError, 'entry 2'):
                molecule_from_json(self.manifest([first, entry]), {})
        data = molecule_from_json(self.manifest([first, {'filename':'a.xyz'}, dict(first, frame=8)]), {})
        self.assertEqual(data['frame_numbers'], [2, 3, 8])
        self.assertEqual(len(data['atoms'][0]['trajectory']), 3)

    def test_sdf_record_selection_and_explicit_topology_survive_merging(self):
        source = (Path(EXAMPLES) / 'bond_orders.sdf').read_text().split('$$$$')[0]
        shifted = source.replace('-3.4000', '-3.9000')
        self._write('a.sdf', source + '$$$$\n' + shifted + '$$$$\n')
        first = {'filename':'a.sdf', 'sdf_record':0}
        second = dict(first, sdf_record=1)
        options = {}
        merged = molecule_from_json(self.manifest([first, second]), options)
        self.assertEqual(options['plot_type'], 'animate')
        self.assertEqual(merged['bonds'][0]['order'], 2)
        self.assertEqual([s['position'][0] for s in merged['atoms'][0]['trajectory']], [-3.4, -3.9])
        self._write('a.sdf', source + '$$$$\n' + shifted.replace('  1  2  2', '  1  2  1') + '$$$$\n')
        with self.assertRaisesRegex(ValueError, 'entry 2.*bond orders'):
            molecule_from_json(self.manifest([first, second]), {})

    def test_molden_frames_keep_orbitals_and_reject_incompatible_basis(self):
        original = (Path(EXAMPLES) / 'water-sto-3g.molden').read_text()
        self._write('a.molden', original)
        self._write('b.molden', original.replace('-20.670655', '-21.0'))
        path = self.manifest([{'filename':x+'.molden'} for x in 'ab'])
        data = molecule_from_json(path, {'orca_fix':True})
        self.assertEqual(len(data['orbital_trajectory']), 1)
        self.assertFalse(data['orbital_trajectory'][0]['denorm'])
        self.assertNotEqual(data['mo'][0]['ene'], data['orbital_trajectory'][0]['mo'][0]['ene'])
        self._write('b.molden', original.replace('3.42525091', '3.52525091'))
        with self.assertRaisesRegex(ValueError, 'entry 2.*basis'):
            molecule_from_json(path, {})

    def test_cube_frames_keep_independent_grid_origins_and_check_dataset_identity(self):
        cube = 'field\ngrid\n-1 0 0 0\n2 1 0 0\n2 0 1 0\n2 0 0 1\n1 0 0 0 0\n1 15\n1 2 3 4 5 6 7 8\n'
        self._write('a.cube', cube)
        self._write('b.cube', cube.replace('-1 0 0 0', '-1 1 0 0'))
        path = self.manifest([{'filename': x+'.cube'} for x in 'ab'])
        data = molecule_from_json(path, {})
        self.assertGreater(data['volume_trajectory'][0]['origin'][0], data['volume']['origin'][0])
        self._write('b.cube', cube.replace('1 15', '1 17'))
        with self.assertRaisesRegex(ValueError, 'entry 2.*dataset IDs'):
            molecule_from_json(path, {})

    def test_sampling_skips_files_and_keeps_manifest_timeline(self):
        self._write('a.xyz', '1\natom\nHe 2 0 0\n')
        path = self.manifest([{'filename': name + '.xyz', 'frame': t} for name, t in
                              [('missing', 2), ('a', 5), ('unused', 7), ('a', 11), ('unused', 20)]])
        data = molecule_from_json(path, dict(frame_start=1, frame_end=3, frame_step=2))
        self.assertEqual(data['frame_numbers'], [5, 11])
        self.assertEqual(len(data['atoms'][0]['trajectory']), 2)
        with self.assertRaisesRegex(ValueError, 'no JSON'):
            molecule_from_json(path, {'frame_start': 5})

    def test_linear_combination_of_cube_datasets_and_units(self):
        import numpy as np
        cube = ('field\ngrid\n1 0 0 0 2\n2 1 0 0\n2 0 1 0\n2 0 0 1\n1 0 0 0 0\n'
                + ' '.join(f'{i} {7-i}' for i in range(8)) + '\n')
        self._write('a.cube', cube)
        self._write('b.cube', cube)
        terms = [dict(filename='a.cube', coefficient=.5, cube_dataset=0),
                 dict(filename='b.cube', coefficient=-2, cube_dataset=1),
                 dict(filename='a.cube', coefficient=.25, cube_dataset=0)]
        def manifest():
            return self.manifest([dict(linear_combination=terms, units='e/bohr^3')])
        volume = molecule_from_json(manifest(), {})['volume']
        np.testing.assert_allclose(volume['data'].ravel(), .75*np.arange(8)-2*np.arange(7, -1, -1))
        self.assertEqual(volume['datasets'].shape, (2, 2, 2, 1))
        self.assertEqual(volume['field_type'], 'signed')
        self.assertEqual(volume['linear_combination'], terms)
        for key, val, message in [('units', 'e/angstrom^3', 'units'), ('coefficient', None, 'coefficient'),
                                   ('cube_dataset', 2, 'cube_dataset')]:
            old = dict(terms[-1])
            terms[-1][key] = val
            with self.assertRaisesRegex(ValueError, 'entry 1.*' + message):
                molecule_from_json(manifest(), {})
            terms[-1] = old
        self._write('b.cube', cube.replace('1 0 0 0 2', '1 1 0 0 2'))
        with self.assertRaisesRegex(ValueError, 'Term 2.*origin'):
            molecule_from_json(manifest(), {})

    def test_resample_combination_on_first_grid_only_with_full_coverage(self):
        import numpy as np
        def cube(name, origin, axes, shape):
            indices = np.indices(shape).reshape(3, -1).T
            x, y, z = (origin + indices @ axes).T
            field = 2*x - y + x*y*z
            rows = ['reference', 'two datasets', '1 ' + ' '.join(map(str, origin)) + ' 2']
            rows += [' '.join(map(str, [-n, *axis])) for n, axis in zip(shape, axes)]
            rows += ['1 0 0 0 0', ' '.join(map(str, np.column_stack((np.ones_like(field), field)).ravel()))]
            self._write(name, '\n'.join(rows) + '\n')
            return field.reshape(shape)
        target = cube('target.cube', np.array([1., .25, .5]), np.diag([.25, .5, .75]), (5, 3, 3))
        axes = np.array([[0., .5, 0.], [1., 0., 0.], [0., 0., 1.]])
        cube('other.cube', np.zeros(3), axes, (4, 3, 3))
        terms = [dict(filename='target.cube', coefficient=1.5, cube_dataset=1),
                 dict(filename='other.cube', coefficient=-.5, cube_dataset=1)]
        def read(**options):
            return molecule_from_json(self.manifest([dict(units='1', linear_combination=terms, **options)]), {})
        with self.assertRaisesRegex(ValueError, 'different grid'):
            read()
        with self.assertRaisesRegex(ValueError, 'true or false'):
            read(resample='true')
        volume = read(resample=True)['volume']
        np.testing.assert_allclose(volume['data'], target, atol=1e-7)
        np.testing.assert_allclose(volume['origin'], [1., .25, .5])
        self.assertTrue(volume['resample'])
        cube('other.cube', np.zeros(3), axes, (3, 3, 3))
        with self.assertRaisesRegex(ValueError, 'Term 2.*cover'):
            read(resample=True)

    def test_json_referencing_json_raises(self):
        f = self._write("bad.json", '{"molecules": [{"filename": "inner.json"}]}')
        with self.assertRaisesRegex(Exception, "json"):
            molecule_from_json(f, {})


# -- Turbomole reader --

class TestPeriodicTable(unittest.TestCase):
    def test_element_lookup(self):
        from molecular_blender.periodictable import elements
        c = elements["c"]
        assert c.name == "Carbon"
        assert c.symbol == "c"
        assert c.mass > 0
        assert c.vdw > 0
        assert c.covalent > 0

    def test_element_str(self):
        from molecular_blender.periodictable import elements
        assert str(elements["h"]) == "h"

    def test_symbols_list(self):
        from molecular_blender.periodictable import symbols
        assert symbols[1] == "h"
        assert symbols[6] == "c"
        assert symbols[8] == "o"


class TestTurbomoleReader(unittest.TestCase):
    def test_dipl_parser(self):
        with open(os.path.join(EXAMPLES, "dipl_a"), "r") as f:
            content = f.read()
        parser = DiplParser(content)
        assert parser.get_section("symmetry") == "c1"
        assert parser.get_section("tensor space dimension") == 8883
        assert parser.get_section("current subspace dimension") == 6

    def test_eigenpairs(self):
        with open(os.path.join(EXAMPLES, "dipl_a"), "r") as f:
            content = f.read()
        parser = DiplParser(content)
        eigenpairs = parser.get_section("eigenpairs")
        assert eigenpairs is not None
        assert len(eigenpairs) > 0
        assert eigenpairs[0].index == 1
        assert eigenpairs[0].eigenvalue > 0

    def test_read_dipl_file(self):
        parser = read_dipl_file(os.path.join(EXAMPLES, "dipl_a"))
        assert parser.get_section("symmetry") == "c1"

    def test_read_polarizability(self):
        pols = read_turbomole_polarizability(os.path.join(EXAMPLES, "dipl_a"))
        assert len(pols) > 0
        assert "eigenvalue" in pols[0]
        assert "vector" in pols[0]

    def test_scientific_notation_parser(self):
        parser = DiplParser("$symmetry c1\n$end\n")
        val = parser._parse_scientific_notation("0.12345D-03")
        assert abs(val - 0.12345e-3) < 1e-10

    def test_fixed_width_numbers(self):
        parser = DiplParser("$symmetry c1\n$end\n")
        line = "0.14237842212942D-03-.22577892444500D-03"
        nums = parser._parse_fixed_width_numbers(line)
        assert len(nums) == 2
        assert abs(nums[0] - 0.14237842212942e-3) < 1e-15
