"""Generate the ethene pi HOMO fixtures with PySCF 2.11.0 and SciPy."""
from pathlib import Path

import numpy as np
from pyscf import gto, lib, scf
from pyscf.tools import molden
from scipy.optimize import brentq


HERE = Path(__file__).resolve().parent
LEVELS = (.05, .005)
# Off-axis probe lines above the two carbons and the middle of the bond, in bohr.
RAYS = ((-1.26, .43), (0., -.37), (1.26, .59))


def calculate(angle):
    """Rotate the right CH2 group around the C=C axis, holding the bond lengths fixed."""
    atoms = [('C', (-.667, 0., 0.)), ('C', (.667, 0., 0.))]
    for side, twist in [(-1, 0.), (1, np.deg2rad(angle))]:
        for sign in (-1, 1):
            atoms.append(('H', (side * 1.238, sign * .923 * np.cos(twist),
                                sign * .923 * np.sin(twist))))
    mol = gto.M(atom=atoms, basis='6-31+g(d)', unit='Angstrom', verbose=0)
    mf = scf.RHF(mol).run(conv_tol=1e-12)
    if not mf.converged:
        raise RuntimeError('Ethene RHF did not converge')
    return mf


def write_reference(stem, mf, index):
    coeff = mf.mo_coeff[:, index]
    def value(x, y, z):
        return (mf.mol.eval_gto('GTOval_sph', [[x, y, z]]) @ coeff).item()

    with (HERE / f'{stem}.molden').open('w') as stream:
        molden.header(mf.mol, stream)
        # Retain the occupied orbitals so HOMO selection uses the normal import path.
        occupied = slice(0, index + 1)
        molden.orbital_coeff(mf.mol, stream, mf.mo_coeff[:, occupied],
                             ene=mf.mo_energy[occupied], occ=mf.mo_occ[occupied])
    samples, crossings = [], []
    for x, y in RAYS:
        node = brentq(lambda z: value(x, y, z), -1., 1.)
        for z in [node - 1.e-3, node, node + 1.e-3, -8., -4., -1., 1., 4., 8.]:
            samples.append((x, y, z, value(x, y, z)))
        for level in LEVELS:
            for side in (-1, 1):
                signed = np.copysign(level, value(x, y, side * 2.))
                ends = sorted((side * 2., side * 12.))
                outer = brentq(lambda z: value(x, y, z) - signed, *ends)
                crossings.append((x, y, signed, outer))
    np.savetxt(HERE / f'{stem}-samples.dat', samples, fmt='%.14g',
               header='x y z (bohr) and signed pi orbital (bohr^-3/2)')
    np.savetxt(HERE / f'{stem}-crossings.dat', crossings, fmt='%.14g',
               header='x y (bohr), signed level (bohr^-3/2), outer crossing z (bohr)')


if __name__ == '__main__':
    lib.num_threads(1)
    planar, twisted = calculate(0), calculate(30)
    homo = np.flatnonzero(planar.mo_occ)[-1]
    anchor = planar.mol.eval_gto('GTOval_sph', [[-1.26, 0., 1.]]) @ planar.mo_coeff[:, homo]
    planar.mo_coeff[:, homo] *= np.sign(anchor.item())
    overlap = planar.mo_coeff[:, homo] @ gto.intor_cross('int1e_ovlp', planar.mol, twisted.mol) @ twisted.mo_coeff
    if np.argmax(np.abs(overlap)) != homo or abs(overlap[homo]) < .8 or twisted.mo_occ[homo] != 2:
        raise ValueError('Twisted HOMO no longer matches the planar pi orbital')
    twisted.mo_coeff[:, homo] *= np.sign(overlap[homo])
    print(f'RHF pi orbital {homo + 1}, aligned overlap {abs(overlap[homo]):.8f}')
    write_reference('planar', planar, homo)
    write_reference('twisted', twisted, homo)
