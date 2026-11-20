"""Tests for molecular_blender.transform"""

import unittest
import math
import tempfile
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

    def test_harmonics_are_orthonormal_in_cartesian_metric(self):
        from molecular_blender.compute_ao_py import Shell
        for l, shell in enumerate('dfg', 2):
            with self.subTest(shell=shell):
                powers = np.array([[term.count(axis) for axis in 'xyz'] for term in Shell.shorder[l]])
                # Even Gaussian moments are proportional to (n-1)!!, odd moments vanish.
                n = powers[:, None, :] + powers[None, :, :]
                moments = np.array([0 if k % 2 else math.prod(range(1, k, 2)) for k in range(2*l+1)])
                gram = moments[n].prod(axis=2)
                metric = gram / np.sqrt(np.outer(gram.diagonal(), gram.diagonal()))
                matrix = get_sph_to_cart(shell)
                np.testing.assert_allclose(matrix.T @ metric @ matrix, np.eye(2*l+1), atol=1e-14)

    def test_standard_molden_f_g_order_phase_and_shape(self):
        from numpy.polynomial.legendre import Legendre
        from molecular_blender.compute_ao_py import polynomial_f, polynomial_g
        points = np.array([[1, 2, 3], [-3, 2, -1], [2, -3, 4], [1, 0, 0], [0, 0, 1]], dtype=float)
        x, y, z = (points / np.linalg.norm(points, axis=1)[:, None]).T
        phi = np.arctan2(y, x)
        for l, shell, polynomial in [(3, 'f', polynomial_f), (4, 'g', polynomial_g)]:
            with self.subTest(shell=shell):
                # Racah-normalized real harmonics on the unit sphere, in Molden
                # order. Legendre derivatives omit the Condon–Shortley phase.
                legendre = Legendre.basis(l)
                expected = [legendre(z)]
                for m in range(1, l+1):
                    factor = math.sqrt(2 * math.factorial(l-m) / math.factorial(l+m))
                    amplitude = factor * (1-z*z)**(m/2) * legendre.deriv(m)(z)
                    expected.extend([amplitude * np.cos(m*phi), amplitude * np.sin(m*phi)])
                mos = [dict(coeff=row) for row in np.eye(2*l+1)]
                transform_sph_to_cart(mos, [shell], [[dict(shell=shell)]])
                actual = np.array([mo['coeff'] for mo in mos]) @ polynomial(x, y, z)
                np.testing.assert_allclose(actual, expected, atol=1e-7)

    def test_orca_water_with_f_and_g_against_plotted_orbitals(self):
        from molecular_blender.importers import molecule_from_molden
        from molecular_blender.orbitals import MOData
        directory = Path(__file__).parent / 'fixtures' / 'orca_water_fg'
        samples = np.loadtxt(directory / 'mo-values.dat')
        model = MOData.from_dict(molecule_from_molden(directory / 'water.molden', {}))
        axis = np.linspace(-3., 3., 13, dtype=np.float32)
        indices = tuple(((samples[:, :3] + 3) * 2).astype(int).T)
        for number, reference in enumerate(samples[:, 3:].T, 1):
            with self.subTest(orbital=number):
                orbital = model.get_orbital(number)
                actual = [orbital.value(*p) for p in samples[:, :3]]
                grid = orbital.box_values(axis, axis, axis)
                # Same SCF solution and phase. Allow the reader's 1e-6 shell
                # screening threshold and ORCA's six significant print digits.
                np.testing.assert_allclose(actual, reference, rtol=6e-6, atol=1e-6)
                np.testing.assert_allclose(grid[indices], reference, rtol=6e-6, atol=1e-6)

    def test_pyscf_spherical_and_cartesian_orbitals_and_total_density(self):
        from molecular_blender.importers import molecule_from_molden
        from molecular_blender.orbitals import MOData
        root = Path(__file__).resolve().parents[1]
        fixture = root / 'test/fixtures/pyscf_water'
        for source, table, count in ((root / 'examples/water_esp/water.molden', 'samples.dat', 28),
                                     (fixture / 'cartesian.molden', 'cartesian.dat', 38)):
            samples = np.loadtxt(fixture / table)
            model = MOData.from_dict(molecule_from_molden(source, {}))
            axes = [np.unique(col).astype(np.float32) for col in samples[:, :3].T]
            self.assertEqual(samples.shape, (64, count + 4))
            self.assertEqual(len(model.occupations), count)
            fields = [model.get_density()] + [model.get_orbital(i) for i in range(1, count + 1)]
            for index, (field, reference) in enumerate(zip(fields, samples[:, 3:].T)):
                with self.subTest(source=source.name, field='density' if index == 0 else index):
                    values = [field.value(*p) for p in samples[:, :3]]
                    np.testing.assert_allclose(values, reference, rtol=2e-6, atol=2e-7)
                    np.testing.assert_allclose(field.box_values(*axes).ravel(), reference, rtol=2e-6, atol=2e-7)

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
    def test_orca_fix_respects_mixed_cartesian_f_and_spherical_g(self):
        from molecular_blender.importers import molecule_from_molden
        # In a [9G] file alone the ten f components are Cartesian.
        coefficients = np.arange(1., 20.)
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'mixed.molden'
            path.write_text('[Molden Format]\n[Atoms] AU\nHe 1 2 0 0 0\n[GTO]\n1 0\n'
                            'f 1 1.0\n0.75 1\ng 1 1.0\n0.55 1\n[9G]\n[MO]\n'
                            'Sym= A\nEne= 0\nSpin= Alpha\nOccup= 0\n'
                            + ''.join(f'{i} {v}\n' for i, v in enumerate(coefficients, 1)))
            geo = molecule_from_molden(path, {'orca_fix': True})
        f, g = geo['basis'][0]
        self.assertEqual(f['contractions'], [1.])
        np.testing.assert_allclose(g['contractions'], [1 / math.sqrt(35)])
        np.testing.assert_array_equal(geo['mo'][0]['coeff'][:10], coefficients[:10])

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
