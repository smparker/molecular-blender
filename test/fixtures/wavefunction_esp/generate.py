"""Regenerate ESP reference values with PySCF 2.11, without importing the add-on."""
from pathlib import Path

import numpy as np
from scipy.optimize import brentq
from pyscf import gto, lib, scf
from pyscf.tools import molden

root = Path(__file__).resolve().parent
lib.num_threads(1)


def references(name, mol, dm):
    if dm.ndim == 3:
        dm = dm.sum(axis=0)
    center = mol.atom_coords().mean(axis=0)
    points = [[1.1, -1.5, .8], [3., -.2, 1.7], [0., 1.e4, 0.]]
    for ray in np.eye(3).tolist() + (-np.eye(3)).tolist():
        direction = np.asarray(ray)

        def residual(distance):
            ao = mol.eval_gto('GTOval', (center + distance * direction)[None, :])[0]
            return ao @ dm @ ao - .001

        points.append(center + brentq(residual, .5, 12.) * direction)
    rows = []
    for point in points:
        with mol.with_rinv_origin(point):
            electronic = np.einsum('ij,ij', dm, mol.intor('int1e_rinv'))
        nuclear = (mol.atom_charges() / np.linalg.norm(mol.atom_coords()-point, axis=1)).sum()
        rows.append([*point, nuclear, electronic, nuclear-electronic])
    np.savetxt(root / (name + '.dat'), rows, fmt='%.14g',
               header='x y z (bohr), nuclear, positive electronic integral, ESP (Hartree/e). Last six points lie on rho=0.001.')


def main():
    # Load the exact committed standard Molden files in PySCF.
    for name, source in [('water', root.parents[2] / 'examples/water_esp/water.molden'),
                         ('cartesian', root.parent / 'pyscf_water/cartesian.molden')]:
        mol, _, coeff, occup, _, _ = molden.load(str(source))
        references(name, mol, (coeff * occup) @ coeff.T)

    # Independent SCF with the same geometry and spherical basis as the ORCA fixture.
    oxygen = gto.basis.load('sto-3g', 'O') + [[3, [.75, 1.]], [4, [.55, 1.]]]
    mol = gto.M(atom='O 0 0 0; H .209 .818 .467; H -.754 -.402 .475',
                basis={'O': oxygen, 'H': 'sto-3g'}, verbose=0)
    references('orca', mol, scf.RHF(mol).run(conv_tol=1e-13).make_rdm1())

    mol = gto.M(atom='He 0 0 0; H 0 0 1.45', unit='Bohr', basis='sto-3g', charge=1, verbose=0)
    mf = scf.RHF(mol).run(conv_tol=1e-13)
    molden.from_scf(mf, str(root / 'heh.molden'))
    references('heh', mol, mf.make_rdm1())


if __name__ == '__main__':
    main()
