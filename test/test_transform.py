"""Tests for molecular_blender.transform"""

import unittest
from pathlib import Path
import numpy as np

from molecular_blender.transform import get_sph_to_cart, transform_sph_to_cart


class TestGetSphToCart(unittest.TestCase):
    def _check_shape(self, shell, cart_size, sph_size):
        mat = get_sph_to_cart(shell)
        assert mat.shape == (cart_size, sph_size)

    def test_shape_s(self): self._check_shape("s", 1, 1)
    def test_shape_p(self): self._check_shape("p", 3, 3)
    def test_shape_d(self): self._check_shape("d", 6, 5)
    def test_shape_f(self): self._check_shape("f", 10, 7)
    def test_shape_g(self): self._check_shape("g", 15, 9)

    def test_s_is_identity(self):
        mat = get_sph_to_cart("s")
        np.testing.assert_array_equal(mat, np.eye(1))

    def test_p_is_identity(self):
        mat = get_sph_to_cart("p")
        np.testing.assert_array_equal(mat, np.eye(3))

    def test_d_normalization(self):
        """The five harmonics are orthonormal in the Cartesian overlap metric."""
        mat = get_sph_to_cart("d")
        metric = np.eye(6)
        metric[:3, :3] = 1 / 3
        np.fill_diagonal(metric, 1)
        np.testing.assert_allclose(mat.T @ metric @ mat, np.eye(5), atol=1e-14)

    def test_ethene_pi_orbital_against_orca_cube(self):
        from molecular_blender.importers import molecule_from_molden
        from molecular_blender.orbitals import MOData
        directory = Path(__file__).parent / 'fixtures' / 'orca_ethene'
        samples = np.loadtxt(directory / 'mo8.dat')
        orbitals = MOData.from_dict(molecule_from_molden(directory / 'ethene.molden', {'orca_fix': True}))
        self.assertEqual(orbitals.homo(), 8)
        pi = orbitals.get_orbital(8)
        actual = np.array([pi.value(*p) for p in samples[:, :3]])
        reference = samples[:, 3]
        sign = -1 if np.dot(reference, actual) < 0 else 1
        np.testing.assert_allclose(actual, sign * reference, atol=1e-7, rtol=5e-6)

    def _check_transform_rank(self, shell, ncart, nsph):
        """Transform matrix should have full column rank (nsph)."""
        mat = get_sph_to_cart(shell)
        rank = np.linalg.matrix_rank(mat, tol=1e-8)
        assert rank == nsph

    def test_transform_rank_d_6_5(self): self._check_transform_rank("d", 6, 5)
    def test_transform_rank_f_10_7(self): self._check_transform_rank("f", 10, 7)
    def test_transform_rank_g_15_9(self): self._check_transform_rank("g", 15, 9)


class TestTransformSphToCart(unittest.TestCase):
    def test_identity_for_non_spherical(self):
        """If no shells are in the spherical list, coefficients should be unchanged."""
        orbitals = [{"coeff": [1.0, 0.0, 0.0, 0.0]}]
        basis = [[{"shell": "s"}, {"shell": "p"}]]
        result = transform_sph_to_cart(orbitals, ["d"], basis)
        np.testing.assert_allclose(result[0]["coeff"], [1.0, 0.0, 0.0, 0.0])

    def test_d_transform_size(self):
        """Transforming d from spherical (5) to cartesian (6) should change coeff length."""
        nao_sph = 1 + 3 + 5  # s + p + d(sph)
        nao_cart = 1 + 3 + 6  # s + p + d(cart)
        orbitals = [{"coeff": list(np.zeros(nao_sph))}]
        basis = [[{"shell": "s"}, {"shell": "p"}, {"shell": "d"}]]
        result = transform_sph_to_cart(orbitals, ["d"], basis)
        assert len(result[0]["coeff"]) == nao_cart

    def test_d_transform_preserves_s_p_coefficients(self):
        """s and p coefficients should pass through the d-transform unchanged."""
        nao_sph = 1 + 3 + 5
        coeffs_sph = list(range(1, nao_sph + 1))
        orbitals = [{"coeff": list(coeffs_sph)}]
        basis = [[{"shell": "s"}, {"shell": "p"}, {"shell": "d"}]]
        result = transform_sph_to_cart(orbitals, ["d"], basis)
        assert result[0]["coeff"][:1] == coeffs_sph[:1]
        assert result[0]["coeff"][1:4] == coeffs_sph[1:4]

    def test_transform_with_molden_coefficient_count(self):
        """After reading a molden file, MO coefficients should have Cartesian AO count."""
        import os
        from molecular_blender.importers import molecule_from_file
        from molecular_blender.transform import cartdegen
        examples = os.path.join(os.path.dirname(__file__), "../examples")
        geo = molecule_from_file(os.path.join(examples, "water-sto-3g.molden"), {})
        expected_nao = sum(
            cartdegen[sh["shell"]]
            for atom_basis in geo["basis"]
            for sh in atom_basis
        )
        for mo in geo["mo"]:
            assert len(mo["coeff"]) == expected_nao
