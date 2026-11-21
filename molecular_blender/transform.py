# -*- coding: utf-8 -*-
#
#  Molecular Blender
#  Filename: transform.py
#  Copyright (C) 2020 Shane Parker, Joshua Szekely
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

"""Classes and functions to transform MOs (partially) using spherical
AOs to using only cartesian AOs (CAOs)."""

from typing import Collection, Sequence

import numpy as np
from math import sqrt

cartdegen = {"s": 1, "sp": 4, "p": 3, "d": 6, "f": 10, "g": 15}
sphdegen = {"s": 1, "sp": 4, "p": 3, "d": 5, "f": 7, "g": 9}

def transform_sph_to_cart(orbitals: list[dict], spherical: Collection[str], basis: Sequence[list],
                          *, orca: bool = False) -> list[dict]:
    """Replace each (nao,) coefficient row in place with its Molden Cartesian expansion."""
    nao = len(orbitals[0]["coeff"])
    ncao = sum([sum([cartdegen[x["shell"]] for x in atom])
                       for atom in basis])
    nmo = len(orbitals)

    sph = np.zeros([nmo, nao], dtype=np.float32)
    cart = np.zeros([nmo, ncao], dtype=np.float32)

    for p, orb in enumerate(orbitals):
        sph[p,:] = orb["coeff"][:]

    iao = 0
    icao = 0
    for atom in basis:
        for shell in atom:
            sh = shell["shell"]

            ndeg = sphdegen[sh] if sh in spherical else cartdegen[sh]
            ncartdeg = cartdegen[sh]
            if sh in spherical:
                sph_to_cart = get_sph_to_cart(sh)
                if orca and sh in ("f", "g"):
                    # In orca_2mkl, both components of |m| = 3, 4 have
                    # the opposite phase to the standard real harmonics.
                    sph_to_cart[:, 5:] *= -1
                cart[:,icao:icao+ncartdeg] = np.einsum("cs,ps->pc", sph_to_cart, sph[:,iao:iao+ndeg])
            else:
                cart[:,icao:icao+ndeg] = sph[:,iao:iao+ndeg]

            iao += ndeg
            icao += ncartdeg

    for p, orb in enumerate(orbitals):
        orb["coeff"] = list(cart[p,:])

    return orbitals

def get_sph_to_cart(sh: str) -> np.ndarray:
    """Expand normalized real harmonics in normalized Cartesian Gaussians.

    Rows follow Molden Cartesian order. Columns follow m = 0, +1, -1, ...
    Regenerate with util/sph_to_cart.py --mo-sph-to-cart.
    """
    if sh=="s":
        return np.eye(1)
    elif sh=="p":
        return np.eye(3)
    elif sh=="d":
        # Expand normalized real harmonics in normalized Cartesian functions.
        # Columns follow Molden m = 0, +1, -1, +2, -2.
        return np.array(
                [[-0.5, 0, 0, sqrt(3)/2, 0],
                 [-0.5, 0, 0, -sqrt(3)/2, 0],
                 [1, 0, 0, 0, 0],
                 [0, 0, 0, 0, 1],
                 [0, 1, 0, 0, 0],
                 [0, 0, 1, 0, 0]]
                )
    elif sh=="f":
        return np.array([
            [0, -sqrt(6)/4, 0, 0, 0, sqrt(10)/4, 0],
            [0, 0, -sqrt(6)/4, 0, 0, 0, -sqrt(10)/4],
            [1, 0, 0, 0, 0, 0, 0],
            [0, -sqrt(30)/20, 0, 0, 0, -3*sqrt(2)/4, 0],
            [0, 0, -sqrt(30)/20, 0, 0, 0, 3*sqrt(2)/4],
            [-3*sqrt(5)/10, 0, 0, sqrt(3)/2, 0, 0, 0],
            [0, sqrt(30)/5, 0, 0, 0, 0, 0],
            [0, 0, sqrt(30)/5, 0, 0, 0, 0],
            [-3*sqrt(5)/10, 0, 0, -sqrt(3)/2, 0, 0, 0],
            [0, 0, 0, 0, 1, 0, 0],
        ])
    elif sh=="g":
        return np.array([
            [3/8, 0, 0, -sqrt(5)/4, 0, 0, 0, sqrt(35)/8, 0],
            [3/8, 0, 0, sqrt(5)/4, 0, 0, 0, sqrt(35)/8, 0],
            [1, 0, 0, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, -sqrt(35)/14, 0, 0, 0, sqrt(5)/2],
            [0, -3*sqrt(70)/28, 0, 0, 0, sqrt(10)/4, 0, 0, 0],
            [0, 0, 0, 0, -sqrt(35)/14, 0, 0, 0, -sqrt(5)/2],
            [0, 0, -3*sqrt(70)/28, 0, 0, 0, -sqrt(10)/4, 0, 0],
            [0, sqrt(70)/7, 0, 0, 0, 0, 0, 0, 0],
            [0, 0, sqrt(70)/7, 0, 0, 0, 0, 0, 0],
            [3*sqrt(105)/140, 0, 0, 0, 0, 0, 0, -3*sqrt(3)/4, 0],
            [-3*sqrt(105)/35, 0, 0, 3*sqrt(21)/14, 0, 0, 0, 0, 0],
            [-3*sqrt(105)/35, 0, 0, -3*sqrt(21)/14, 0, 0, 0, 0, 0],
            [0, 0, -3*sqrt(14)/28, 0, 0, 0, 3*sqrt(2)/4, 0, 0],
            [0, -3*sqrt(14)/28, 0, 0, 0, -3*sqrt(2)/4, 0, 0, 0],
            [0, 0, 0, 0, 3*sqrt(7)/7, 0, 0, 0, 0],
        ])
