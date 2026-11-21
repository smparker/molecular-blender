# -*- coding: utf-8 -*-
#
#  Molecular Blender
#  Filename: compute_ao.py
#  Copyright (C) 2014-2025 Shane Parker, Joshua Szekely
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

"""Evaluate real-space orbitals and densities in atomic units.

Coordinates passed to value(), plane_values() and box_values() are in Bohr.
Grid axes for the compiled API are 1-D float32 arrays. Strided axes are
packed before native grid evaluation. Results
are indexed [x, y] or [x, y, z] in the same order as the supplied axes.
"""

CYTHON_ENABLED = False
try:
    from .compute_ao_cy import *
    CYTHON_ENABLED = True
except ModuleNotFoundError as exc:
    if exc.name != __package__ + '.compute_ao_cy':
        raise
    print("Using the Python atomic-orbital evaluator.")
    from .compute_ao_py import *

# Both spellings refer to the same backend classes, including for old callers.
OrbitalCalculator = OrbitalCalculater
DensityCalculator = DensityCalculater
