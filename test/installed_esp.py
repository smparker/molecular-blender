"""Shared scientific check for the installed wheel and Blender ZIP."""
from pathlib import Path

import numpy as np


def check_water_esp():
    from molecular_blender import esp_kernel, wavefunction_esp
    from molecular_blender.importers import molecule_from_molden

    assert wavefunction_esp.coulomb_sum is esp_kernel.coulomb_sum, 'ESP selected the NumPy fallback'
    root = Path(__file__).resolve().parents[1]
    sample = np.loadtxt(root / 'test/fixtures/wavefunction_esp/water.dat')[0]
    data = molecule_from_molden(root / 'examples/water_esp/water.molden', {})
    field = wavefunction_esp.WavefunctionESP(data)
    value, = field.potential_bohr(sample[None, :3])
    np.testing.assert_allclose(value, sample[-1], atol=5e-7, rtol=0)
    print(f'Native water ESP {value:.12f} Hartree/e (reference {sample[-1]:.12f})', flush=True)
