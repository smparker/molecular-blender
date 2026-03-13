"""Tests for molecular_blender.importers and molecular_blender.turbomole_reader"""

import os
import pytest
import numpy as np
import tempfile

from molecular_blender.importers import (
    Reader, molecule_from_file, molecule_from_xyz, molecule_from_cube,
    molecule_from_molden, molecule_from_json, molecule_from_babel,
    make_atom_dict, make_snap_dict,
)
from molecular_blender.turbomole_reader import (
    DiplParser, read_dipl_file, read_turbomole_polarizability,
)

EXAMPLES = os.path.join(os.path.dirname(__file__), "../examples")


# -- helper constructors --

class TestMakeAtomDict:
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


class TestMakeSnapDict:
    def test_basic(self):
        d = make_snap_dict([1.0, 2.0, 3.0], 0.5, [0.1, 0.2, 0.3])
        assert d["position"] == [1.0, 2.0, 3.0]
        assert d["charge"] == 0.5
        assert d["gradient"] == [0.1, 0.2, 0.3]


# -- Reader class --

class TestReader:
    def test_context_manager(self, tmp_path):
        f = tmp_path / "test.txt"
        f.write_text("line1\nline2\n")
        with Reader(str(f)) as r:
            assert r.readline().strip() == "line1"
            assert r.readline().strip() == "line2"

    def test_eof_raises(self, tmp_path):
        f = tmp_path / "test.txt"
        f.write_text("only\n")
        with Reader(str(f)) as r:
            r.readline()
            with pytest.raises(EOFError):
                r.readline()

    def test_marks(self, tmp_path):
        f = tmp_path / "test.txt"
        f.write_text("aaa\nbbb\nccc\n")
        with Reader(str(f)) as r:
            line = r.readline()
            r.set_mark("first", line)
            r.readline()
            r.readline()
            assert r.is_marked("first")
            assert not r.is_marked("nonexistent")
            restored = r.restore_mark("first")
            assert restored.strip() == "aaa"


# -- molecule_from_file dispatch --

class TestMoleculeFromFileDispatch:
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
        assert len(geo["atoms"]) > 0

    def test_babel_raises(self):
        with pytest.raises(Exception, match="babel"):
            molecule_from_file("something.pdb", {})


# -- XYZ reading --

class TestXYZRead:
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

    def test_mismatched_atoms_raises(self, tmp_path):
        f = tmp_path / "bad.xyz"
        f.write_text("2\nframe1\nC 0 0 0\nH 1 0 0\n2\nframe2\nC 0 0 0\nC 1 0 0\n")
        with pytest.raises(Exception, match="order of the atoms"):
            molecule_from_xyz(str(f), {"plot_type": "animate"})

    def test_mismatched_natoms_raises(self, tmp_path):
        f = tmp_path / "bad.xyz"
        f.write_text("2\nframe1\nC 0 0 0\nH 1 0 0\n3\nframe2\nC 0 0 0\nH 1 0 0\nO 2 0 0\n")
        with pytest.raises(Exception, match="same number of atoms"):
            molecule_from_xyz(str(f), {"plot_type": "animate"})


# -- Cube reading --

class TestCubeRead:
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

class TestMoldenRead:
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

    def test_missing_molden_format(self, tmp_path):
        f = tmp_path / "bad.molden"
        f.write_text("[Atoms]\nC 1 6 0.0 0.0 0.0\n")
        with pytest.raises(Exception, match="Molden Format"):
            molecule_from_molden(str(f), {})

    def test_missing_atoms(self, tmp_path):
        f = tmp_path / "bad.molden"
        f.write_text("[Molden Format]\n[GTO]\n")
        with pytest.raises(Exception, match="Atoms"):
            molecule_from_molden(str(f), {})

    def test_orca_fix(self):
        geo = molecule_from_molden(
            os.path.join(EXAMPLES, "water-sto-3g.molden"),
            {"orca_fix": True}
        )
        assert geo["denorm"] is False


# -- JSON reading --

class TestJsonRead:
    def test_phenanth_polar(self):
        geo = molecule_from_json(os.path.join(EXAMPLES, "phenanth_polar.json"), {})
        assert len(geo["atoms"]) > 0
        assert "polar" in geo

    def test_returns_none_for_missing_molecules(self, tmp_path):
        f = tmp_path / "empty.json"
        f.write_text('{"other": "stuff"}')
        result = molecule_from_json(str(f), {})
        assert result is None

    def test_json_referencing_json_raises(self, tmp_path):
        f = tmp_path / "bad.json"
        f.write_text('{"molecules": [{"filename": "inner.json"}]}')
        with pytest.raises(Exception, match="json"):
            molecule_from_json(str(f), {})


# -- Turbomole reader --

class TestPeriodicTable:
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


class TestTurbomoleReader:
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
