"""Regenerate the signed orbital references and the matching ESP figure inputs with PySCF 2.11."""
from pathlib import Path

import numpy as np
from pyscf import gto, lib, scf
from pyscf.dft import numint
from pyscf.tools import cubegen, molden

lib.num_threads(1)
fixture = Path(__file__).resolve().parent
example = fixture.parents[2] / 'examples' / 'water_esp'
example.mkdir(exist_ok=True)
oxygen = gto.basis.load('sto-3g', 'O')
oxygen += [[l, [1.2, 0.4], [0.45, 0.7]] for l in (2, 3, 4)]
mol = gto.M(atom='O 0 0 0; H 0.209 0.818 0.467; H -0.754 -0.402 0.475',
            unit='Angstrom', basis={'O': oxygen, 'H': 'sto-3g'}, cart=False)
mf = scf.RHF(mol).run(conv_tol=1e-12)
if not mf.converged:
    raise RuntimeError('RHF failed to converge')

molden.from_scf(mf, example / 'water.molden')
dm = mf.make_rdm1()
cubegen.mep(mol, example / 'potential.cube', dm, nx=28, ny=28, nz=28, margin=4.)

# A fixed off-axis grid samples all octants, including angular nodes and tails.
points = lib.cartesian_prod([[-2.1, -.8, .3, 1.7], [-1.8, -.2, .6, 2.3], [-2.4, -.5, .9, 1.6]])
ao = mol.eval_gto('GTOval_sph', points)
values = np.column_stack((points, numint.eval_rho(mol, ao, dm), ao @ mf.mo_coeff))
np.savetxt(fixture / 'samples.dat', values, fmt='%.12g',
           header='x y z (bohr), total density (e/bohr^3), then signed MOs in file order (bohr^-3/2)')
