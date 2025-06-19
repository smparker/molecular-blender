# -*- coding: utf-8 -*-
#
#  Molecular Blender
#  Filename: orbitals.py
#  Copyright (C) 2017-2025 Shane Parker, Joshua Szekely
#
#  This file is part of Molecular Blender.
#
#  Molecular Blender is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 3, or (at your option)
#  any later version.
#
#  Molecular Blender is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU Library General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with Molecular Blender; see COPYING.
#  If not, see <http://www.gnu.org/licenses/>.
#
import numpy as np

from .compute_ao import *

"""Classes and functions to compute orbital values in real space."""

class MOData:
    """Organizes data for computing orbitals in real space"""

    def __init__(self, shells, coeff, nocc, occupations, spins=None, densities=None):
        assert sum([sh.size for sh in shells]) == coeff.shape[1]

        self.shells = shells
        self.coeff = coeff
        self.nocc = nocc
        self.nvir = coeff.shape[0] - nocc
        self.occupations = occupations
        self.spins = spins
        # map of density names to AO densities for things like polarizabilities
        self.densities = {} if densities is None else densities

    def homo(self):
        return self.nocc

    def get_orbital(self, iorb):
        """Returns OrbitalCalculater for orbital"""
        i = (iorb if iorb > 0 else self.nocc + iorb) - 1
        return OrbitalCalculater(self.shells, self.coeff[i, :])

    def get_density(self, densityname="density"):
        """Returns DensityCalculater for the electronic density"""
        if densityname=="density":
            rho = np.einsum("pu,p,pv->uv", self.coeff, self.occupations, self.coeff).astype(np.float32)
            nelec = np.sum(self.occupations)
            return DensityCalculater(self.shells, rho, nelec)
        elif densityname == "spin":
            alphaoccs = np.array( [ self.occupations[i] if self.spins[i] == 'Alpha' else 0.0 for i in range(len(self.occupations)) ] )
            betaoccs = np.array( [ self.occupations[i] if self.spins[i] == 'Beta' else 0.0 for i in range(len(self.occupations)) ] )
            alpha = np.einsum("pu,p,pv->uv", self.coeff, alphaoccs, self.coeff)
            beta = np.einsum("pu,p,pv->uv", self.coeff, betaoccs, self.coeff)
            return DensityCalculater(self.shells, alpha - beta, np.sum(alphaoccs) + np.sum(betaoccs))
        elif densityname == "alpha":
            alphaoccs = self.occupations[self.spins == 'Alpha']
            alpha = np.einsum("pu,p,pv->uv", self.coeff, alphaoccs, self.coeff)
            return DensityCalculater(self.shells, alpha, np.sum(alphaoccs))
        elif densityname == "beta":
            betaoccs = self.occupations[self.spins == 'Beta']
            beta = np.einsum("pu,p,pv->uv", self.coeff, betaoccs, self.coeff)
            return DensityCalculater(self.shells, beta, np.sum(betaoccs))
        elif densityname in self.densities:
            rho = self.densities[densityname]
            nelec = np.sum(self.occupations)
            return DensityCalculater(self.shells, rho, nelec)

    def is_density(self, densityname):
        return densityname in ["density", "spin", "alpha", "beta"] or densityname in self.densities.keys()

    @classmethod
    def from_dict(cls, geo):
        """Builds MOData out of the dict object available in MolecularBlender"""
        shells = []
        ico = 0
        for atom, ash in zip(geo["atoms"], geo["basis"]):
            # make sure centers are given in bohr
            xyz = [x * ang2bohr for x in atom["position"]]
            for sh in ash:
                newshell = Shell(xyz, "spdfg".index(
                    sh["shell"]), sh["exponents"], sh["contractions"], start=ico)
                ico += newshell.size
                shells.append(newshell)

        occupations = np.array([ float(x["occup"]) for x in geo["mo"]])
        spins = [ x["spin"] for x in geo["mo"] ]
        nocc = np.sum(occupations > 0.0)
        nmo = len(occupations)
        nao = sum([sh.size for sh in shells])

        coeff = np.zeros([nmo, nao], dtype=np.float32)
        for i in range(nmo):
            coeff[i, :] = geo["mo"][i]["coeff"]

        densities = {}
        if "polar" in geo:
            # turbomole polarizabilities in occ-virt MO are stored in geo
            nvir = nmo - nocc
            for i, pair in enumerate(geo["polar"]):
                eigenvalue = pair["eigenvalue"]
                if False:
                    vec = pair["vector"].reshape(nvir, nocc)
                    dens = np.einsum("au,ai,iv->uv", coeff[nocc:, :], vec, coeff[:nocc, :])
                    densities[f"polar{i:2d}"] = 2.0 * (dens + dens.T)
                else:
                    vec = pair["vector"].reshape(nocc, nvir)
                    dens = np.einsum("iu,ia,av->uv", coeff[:nocc, :], vec, coeff[nocc:, :])
                    densities[f"polar{i:2d}"] = 2.0*(dens + dens.T)

        return cls(shells, coeff, nocc, occupations, spins, densities)
