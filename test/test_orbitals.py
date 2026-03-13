"""Tests for molecular_blender.orbitals (the wrapper module using compute_ao)"""

import os
import pytest
import numpy as np

from molecular_blender.importers import molecule_from_file
from molecular_blender.orbitals import MOData

EXAMPLES = os.path.join(os.path.dirname(__file__), "../examples")


@pytest.fixture
def water_modata():
    geo = molecule_from_file(os.path.join(EXAMPLES, "water-sto-3g.molden"), {})
    return MOData.from_dict(geo)


@pytest.fixture
def fake_modata():
    geo = molecule_from_file(os.path.join(EXAMPLES, "fake-atom.molden"), {})
    return MOData.from_dict(geo)


@pytest.fixture
def phenanth_modata():
    geo = molecule_from_file(os.path.join(EXAMPLES, "phenanth_polar.json"), {})
    return MOData.from_dict(geo)


class TestMODataBasics:
    def test_homo(self, water_modata):
        h = water_modata.homo()
        assert h > 0
        assert h == water_modata.nocc

    def test_nvir(self, water_modata):
        assert water_modata.nvir >= 0
        assert water_modata.nocc + water_modata.nvir == water_modata.coeff.shape[0]

    def test_get_orbital_positive(self, water_modata):
        orb = water_modata.get_orbital(1)
        assert orb is not None

    def test_get_orbital_negative(self, water_modata):
        orb = water_modata.get_orbital(-1)
        assert orb is not None


class TestMODataDensity:
    def test_density(self, fake_modata):
        dens = fake_modata.get_density("density")
        assert dens is not None

    def test_is_density(self, fake_modata):
        assert fake_modata.is_density("density")
        assert not fake_modata.is_density("nonexistent")

    def test_density_norm(self, fake_modata):
        density = fake_modata.get_density()
        xyz = np.linspace(-30, 30, 100, dtype=np.float32)
        dvol = (xyz[1] - xyz[0]) ** 3
        integration = np.sum(density.box_values(xyz, xyz, xyz)) * dvol
        assert abs(integration - 4.0) < 0.5

    def test_orbital_norm(self, fake_modata):
        for i in range(1, 6):
            orb = fake_modata.get_orbital(i)
            xyz = np.linspace(-30, 30, 100, dtype=np.float32)
            dvol = (xyz[1] - xyz[0]) ** 3
            integration = np.sum(orb.box_values(xyz, xyz, xyz) ** 2) * dvol
            assert abs(integration - 1.0) < 0.2, f"Orbital {i} norm = {integration}"


class TestMODataPolarizability:
    def test_polar_density_available(self, phenanth_modata):
        assert phenanth_modata.is_density("polar 0")

    def test_polar_density(self, phenanth_modata):
        polar = phenanth_modata.get_density("polar 0")
        assert polar is not None
        xyz = np.linspace(-30, 30, 100, dtype=np.float32)
        dvol = (xyz[1] - xyz[0]) ** 3
        pxyz = polar.box_values(xyz, xyz, xyz)
        integrate = np.sum(pxyz) * dvol
        assert abs(integrate) < 0.5


class TestMODataSpinDensities:
    def test_spin_density(self, water_modata):
        spin = water_modata.get_density("spin")
        assert spin is not None

    def test_alpha_density(self, water_modata):
        alpha = water_modata.get_density("alpha")
        assert alpha is not None

    def test_beta_density(self, water_modata):
        beta = water_modata.get_density("beta")
        assert beta is not None

    def test_is_density_spin(self, water_modata):
        assert water_modata.is_density("spin")
        assert water_modata.is_density("alpha")
        assert water_modata.is_density("beta")
        assert water_modata.is_density("density")
        assert not water_modata.is_density("garbage")


class TestMODataFromDict:
    def test_denorm_flag(self):
        geo = molecule_from_file(os.path.join(EXAMPLES, "water-sto-3g.molden"),
                                 {"orca_fix": True})
        modata = MOData.from_dict(geo)
        assert modata is not None

    def test_coeff_shape(self, water_modata):
        nmo = water_modata.coeff.shape[0]
        nao = water_modata.coeff.shape[1]
        assert nmo > 0
        assert nao > 0
        assert nao == sum(sh.size for sh in water_modata.shells)


class TestOrbitalIsovalue:
    def test_cumulative_orbital(self, fake_modata):
        orb = fake_modata.get_orbital(1)
        isovals = orb.isovalue_containing_proportion(
            [0.5, 0.9], resolution=0.25, box=[[-3, -3, -3], [3, 3, 3]]
        )
        assert len(isovals) == 2
        assert isovals[0] > isovals[1]

    def test_cumulative_density(self, fake_modata):
        density = fake_modata.get_density()
        isovals = density.isovalue_containing_proportion(
            [0.5, 0.9], resolution=0.25, box=([-3, -3, -3], [3, 3, 3])
        )
        assert len(isovals) == 2
        assert isovals[0] > isovals[1]
