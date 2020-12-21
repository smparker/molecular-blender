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

import numpy as np
from math import sqrt

cartdegen = {"s": 1, "sp": 4, "p": 3, "d": 6, "f": 10, "g": 15}
sphdegen = {"s": 1, "sp": 4, "p": 3, "d": 5, "f": 7, "g": 9}

def transform_sph_to_cart(orbitals, spherical, basis):
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
                # do transform
                sph_to_cart = get_sph_to_cart(sh)
                cart[:,icao:icao+ncartdeg] = np.einsum("cs,ps->pc", sph_to_cart, sph[:,iao:iao+ndeg])
            else:
                cart[:,icao:icao+ndeg] = sph[:,iao:iao+ndeg]

            iao += ndeg
            icao += ncartdeg

    for p, orb in enumerate(orbitals):
        orb["coeff"] = list(cart[p,:])

    return orbitals

def get_sph_to_cart(sh):
    """Returns real spherical to cartesian GTO transformation

    Data was generated with sph_to_cart.py script
    """
    if sh=="s":
        return np.eye(1)
    elif sh=="p":
        return np.eye(3)
    elif sh=="d":
        return np.array(
                [[-1/3, 0, 0, (1/3)*sqrt(3), 0],
                 [-1/3, 0, 0,  -1/3*sqrt(3), 0],
                 [2/3, 0, 0, 0, 0],
                 [0, 0, 0, 0, 1],
                 [0, 1, 0, 0, 0],
                 [0, 0, 1, 0, 0]]
                )
    elif sh=="f":
        return np.array(
                [[0, -1/10*sqrt(6), 0, 0, 0, (1/10)*sqrt(10), 0],
                 [0, 0, -1/10*sqrt(6), 0, 0, 0, -1/10*sqrt(10)],
                 [2/5, 0, 0, 0, 0, 0, 0],
                 [0, -1/30*sqrt(30), 0, 0, 0, -1/2*sqrt(2), 0],
                 [0, 0, -1/30*sqrt(30), 0, 0, 0, (1/2)*sqrt(2)],
                 [-1/5*sqrt(5), 0, 0, (1/3)*sqrt(3), 0, 0, 0],
                 [0, (2/15)*sqrt(30), 0, 0, 0, 0, 0],
                 [0, 0, (2/15)*sqrt(30), 0, 0, 0, 0],
                 [-1/5*sqrt(5), 0, 0, -1/3*sqrt(3), 0, 0, 0],
                 [0, 0, 0, 0, 1, 0, 0]]
                )
    elif sh=="g":
        return np.array(
                [[3/35, 0, 0, -2/35*sqrt(5), 0, 0, 0, (1/35)*sqrt(35), 0],
                 [3/35, 0, 0, (2/35)*sqrt(5), 0, 0, 0, (1/35)*sqrt(35), 0],
                 [8/35, 0, 0, 0, 0, 0, 0, 0, 0],
                 [0, 0, 0, 0, -1/35*sqrt(35), 0, 0, 0, (1/5)*sqrt(5)],
                 [0, -3/70*sqrt(70), 0, 0, 0, (1/10)*sqrt(10), 0, 0, 0],
                 [0, 0, 0, 0, -1/35*sqrt(35), 0, 0, 0, -1/5*sqrt(5)],
                 [0, 0, -3/70*sqrt(70), 0, 0, 0, -1/10*sqrt(10), 0, 0],
                 [0, (2/35)*sqrt(70), 0, 0, 0, 0, 0, 0, 0],
                 [0, 0, (2/35)*sqrt(70), 0, 0, 0, 0, 0, 0],
                 [(1/105)*sqrt(105), 0, 0, 0, 0, 0, 0, -1/3*sqrt(3), 0],
                 [-4/105*sqrt(105), 0, 0, (2/21)*sqrt(21), 0, 0, 0, 0, 0],
                 [-4/105*sqrt(105), 0, 0, -2/21*sqrt(21), 0, 0, 0, 0, 0],
                 [0, 0, -1/14*sqrt(14), 0, 0, 0, (1/2)*sqrt(2), 0, 0],
                 [0, -1/14*sqrt(14), 0, 0, 0, -1/2*sqrt(2), 0, 0, 0],
                 [0, 0, 0, 0, (2/7)*sqrt(7), 0, 0, 0, 0]]
                )
