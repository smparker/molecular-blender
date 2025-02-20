# -*- coding: utf-8 -*-
#
#  Molecular Blender
#  Filename: orbitals.py
#  Copyright (C) 2017-2024 Shane Parker, Joshua Szekely
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

    def __init__(self, shells, coeff, nocc, occupations):
        assert sum([sh.size for sh in shells]) == coeff.shape[1]

        self.shells = shells
        self.coeff = coeff
        self.nocc = nocc
        self.nvir = coeff.shape[0] - nocc
        self.occupations = occupations

    def homo(self):
        return self.nocc

    def get_orbital(self, iorb):
        """Returns OrbitalCalculater for orbital"""
        i = (iorb if iorb > 0 else self.nocc + iorb) - 1
        return OrbitalCalculater(self.shells, self.coeff[i, :])

    def get_density(self):
        """Returns DensityCalculater for the electronic density"""
        rho = np.einsum("pu,p,pv->uv", self.coeff, self.occupations, self.coeff).astype(np.float32)
        nelec = np.sum(self.occupations)
        return DensityCalculater(self.shells, rho, nelec)

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
        nocc = np.sum(occupations > 0.0)
        nmo = len(occupations)
        nao = sum([sh.size for sh in shells])

        coeff = np.zeros([nmo, nao])
        for i in range(nmo):
            coeff[i, :] = geo["mo"][i]["coeff"]

        return cls(shells, coeff, nocc, occupations)
