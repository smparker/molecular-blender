from pathlib import Path
import unittest
from unittest.mock import patch

import numpy as np

from molecular_blender import wavefunction_esp as esp
from molecular_blender.importers import molecule_from_molden

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
