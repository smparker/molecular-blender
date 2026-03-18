"""Tests for molecular_blender.orbitals (the wrapper module using compute_ao)"""

import os
import unittest
import numpy as np

from molecular_blender.importers import molecule_from_file
from molecular_blender.orbitals import MOData

EXAMPLES = os.path.join(os.path.dirname(__file__), "../examples")


class TestMODataBasics(unittest.TestCase):
    def setUp(self):
        geo = molecule_from_file(os.path.join(EXAMPLES, "water-sto-3g.molden"), {})
        self.modata = MOData.from_dict(geo)

    def test_homo(self):
        h = self.modata.homo()
        assert h > 0
        assert h == self.modata.nocc

    def test_nvir(self):
        assert self.modata.nvir >= 0
        assert self.modata.nocc + self.modata.nvir == self.modata.coeff.shape[0]

    def test_different_orbitals_give_different_values(self):
        """get_orbital(1) and get_orbital(2) should produce distinct orbital values."""
        orb1 = self.modata.get_orbital(1)
        orb2 = self.modata.get_orbital(2)
        v1 = orb1.value(0.0, 0.0, 0.2)
        v2 = orb2.value(0.0, 0.0, 0.2)
        assert abs(v1 - v2) > 1e-6, "Different orbitals should have different values"

    def test_get_orbital_homo_indexing(self):
        """get_orbital(0) = HOMO = get_orbital(nocc)."""
        homo_orb = self.modata.get_orbital(0)
        homo_abs = self.modata.get_orbital(self.modata.nocc)
        v_rel = homo_orb.value(0.0, 0.1, 0.0)
        v_abs = homo_abs.value(0.0, 0.1, 0.0)
        assert abs(v_rel - v_abs) < 1e-6


class TestMODataDensity(unittest.TestCase):
    def setUp(self):
        geo = molecule_from_file(os.path.join(EXAMPLES, "fake-atom.molden"), {})
        self.modata = MOData.from_dict(geo)

    def test_density_nonnegative_on_grid(self):
        """Electron density should be non-negative everywhere."""
        dens = self.modata.get_density("density")
        xyz = np.linspace(-5, 5, 10, dtype=np.float32)
        box = dens.box_values(xyz, xyz, xyz)
        assert np.all(box >= -1e-6), f"Density has negative values: min={np.min(box)}"

    def test_density_peaks_at_nucleus(self):
        """For an s-type atom, density should peak near the nucleus."""
        dens = self.modata.get_density("density")
        xyz_near = np.array([0.0], dtype=np.float32)
        xyz_far = np.array([5.0], dtype=np.float32)
        val_origin = dens.box_values(xyz_near, xyz_near, xyz_near)[0, 0, 0]
        val_far = dens.box_values(xyz_far, xyz_far, xyz_far)[0, 0, 0]
        assert val_origin > val_far

    def test_is_density(self):
        assert self.modata.is_density("density")
        assert not self.modata.is_density("nonexistent")

    def test_density_norm(self):
        density = self.modata.get_density()
        xyz = np.linspace(-30, 30, 100, dtype=np.float32)
        dvol = (xyz[1] - xyz[0]) ** 3
        integration = np.sum(density.box_values(xyz, xyz, xyz)) * dvol
        assert abs(integration - 4.0) < 0.5

    def test_orbital_norm(self):
        for i in range(1, 6):
            orb = self.modata.get_orbital(i)
            xyz = np.linspace(-30, 30, 100, dtype=np.float32)
            dvol = (xyz[1] - xyz[0]) ** 3
            integration = np.sum(orb.box_values(xyz, xyz, xyz) ** 2) * dvol
            assert abs(integration - 1.0) < 0.2, f"Orbital {i} norm = {integration}"


class TestMODataPolarizability(unittest.TestCase):
    def setUp(self):
        geo = molecule_from_file(os.path.join(EXAMPLES, "phenanth_polar.json"), {})
        self.modata = MOData.from_dict(geo)

    def test_polar_density_available(self):
        assert self.modata.is_density("polar 0")

    def test_polar_density(self):
        polar = self.modata.get_density("polar 0")
        assert polar is not None
        xyz = np.linspace(-30, 30, 100, dtype=np.float32)
        dvol = (xyz[1] - xyz[0]) ** 3
        pxyz = polar.box_values(xyz, xyz, xyz)
        integrate = np.sum(pxyz) * dvol
        assert abs(integrate) < 0.5


class TestMODataSpinDensities(unittest.TestCase):
    def setUp(self):
        geo = molecule_from_file(os.path.join(EXAMPLES, "water-sto-3g.molden"), {})
        self.modata = MOData.from_dict(geo)

    def test_alpha_plus_beta_equals_total_on_grid(self):
        """Alpha + beta density values should equal total density at every grid point."""
        alpha = self.modata.get_density("alpha")
        beta = self.modata.get_density("beta")
        total = self.modata.get_density("density")
        xyz = np.linspace(-5, 5, 10, dtype=np.float32)
        alpha_box = alpha.box_values(xyz, xyz, xyz)
        beta_box = beta.box_values(xyz, xyz, xyz)
        total_box = total.box_values(xyz, xyz, xyz)
        np.testing.assert_allclose(
            alpha_box + beta_box, total_box, atol=1e-4
        )

    def test_spin_density_is_alpha_minus_beta_on_grid(self):
        """Spin density = alpha - beta at every grid point."""
        alpha = self.modata.get_density("alpha")
        beta = self.modata.get_density("beta")
        spin = self.modata.get_density("spin")
        xyz = np.linspace(-5, 5, 10, dtype=np.float32)
        alpha_box = alpha.box_values(xyz, xyz, xyz)
        beta_box = beta.box_values(xyz, xyz, xyz)
        spin_box = spin.box_values(xyz, xyz, xyz)
        np.testing.assert_allclose(
            alpha_box - beta_box, spin_box, atol=1e-4
        )

    def test_is_density_spin(self):
        assert self.modata.is_density("spin")
        assert self.modata.is_density("alpha")
        assert self.modata.is_density("beta")
        assert self.modata.is_density("density")
        assert not self.modata.is_density("garbage")


class TestMODataFromDict(unittest.TestCase):
    def setUp(self):
        geo = molecule_from_file(os.path.join(EXAMPLES, "water-sto-3g.molden"), {})
        self.modata = MOData.from_dict(geo)

    def test_orca_fix_changes_orbital_values(self):
        """orca_fix should change denormalization, producing different orbital values."""
        geo_normal = molecule_from_file(os.path.join(EXAMPLES, "water-sto-3g.molden"), {})
        geo_orca = molecule_from_file(os.path.join(EXAMPLES, "water-sto-3g.molden"),
                                      {"orca_fix": True})
        modata_normal = MOData.from_dict(geo_normal)
        modata_orca = MOData.from_dict(geo_orca)
        orb_n = modata_normal.get_orbital(1)
        orb_o = modata_orca.get_orbital(1)
        val_n = orb_n.value(0.0, 0.0, 0.0)
        val_o = orb_o.value(0.0, 0.0, 0.0)
        assert abs(val_n - val_o) > 1e-6, "orca_fix should change orbital values"

    def test_coeff_shape(self):
        nmo = self.modata.coeff.shape[0]
        nao = self.modata.coeff.shape[1]
        assert nmo > 0
        assert nao > 0
        assert nao == sum(sh.size for sh in self.modata.shells)


class TestOrbitalIsovalue(unittest.TestCase):
    def setUp(self):
        geo = molecule_from_file(os.path.join(EXAMPLES, "fake-atom.molden"), {})
        self.modata = MOData.from_dict(geo)

    def test_cumulative_orbital(self):
        orb = self.modata.get_orbital(1)
        isovals = orb.isovalue_containing_proportion(
            [0.5, 0.9], resolution=0.25, box=[[-3, -3, -3], [3, 3, 3]]
        )
        assert len(isovals) == 2
        assert isovals[0] > isovals[1]

    def test_cumulative_density(self):
        density = self.modata.get_density()
        isovals = density.isovalue_containing_proportion(
            [0.5, 0.9], resolution=0.25, box=([-3, -3, -3], [3, 3, 3])
        )
        assert len(isovals) == 2
        assert isovals[0] > isovals[1]
