# -*- coding: utf-8 -*-
#
#  Molecular Blender
#  Filename: containing_isovalues.py
#  Copyright (C) 2025 Shane Parker
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
import math
from .util import round_sigfigs

def isovalue_containing_proportion(values, volume_data_sh, dV, square=True, method='binary'):
    """
    Find the isovalue that contains a target percentage of the total integrated function value.
    """
    if method not in ['binary', 'sorted']:
        method = 'binary'
    if method == 'binary':
        return isovalue_containing_proportion_binary(values, volume_data_sh, dV, square)
    else: # method == 'sorted'
        return isovalue_containing_proportion_sorted(values, volume_data_sh, dV, square)

def isovalue_containing_proportion_sorted(values, volume_data_sh, dV, square=True):
    cumulative = 0.0
    volume_data = volume_data_sh.reshape(-1)
    data = np.zeros([len(volume_data), 2], dtype=np.float32)
    vkey = (lambda x: x*x) if square else (lambda x: abs(x))
    for i, v in enumerate(sorted(volume_data, key=vkey, reverse=True)):
        data[i,0] = abs(v)
        data[i,1] = cumulative
        cumulative += vkey(v) * dV

    # rescale containing values down according to integrated density
    # equivalent to uniformly scaling cumulative data so norm is 1
    values = [ cumulative * v for v in values ]
    isovals = [ None for x in values ]

    for i in range(1, data.shape[0]):
        x1, y1 = data[i-1,:]
        x2, y2 = data[i,:]

        for ival, val in enumerate(values):
            # does desired integration fall between y1 and y2
            if y2 >= abs(val) and y1 < abs(val):
                # linearly interpolate measured value
                print("x1: {}, y1: {}, x2: {}, y2: {}, val: {}".format(x1, y1, x2, y2, val))
                xval = x1 + (abs(val) - y1)/(y2 - y1)*(x2 - x1)
                isovals[ival] = math.copysign(xval, val)

    # round isovalues to 2 decimal places
    isovals = [ round_sigfigs(x, 2) for x in isovals ]

    return isovals

def isovalue_containing_proportion_binary(values, volume_data_sh, dV, square=True,
                                   tolerance=1e-3, max_iterations=100):
    """
    Find the isovalue that contains a target percentage of the total integrated function value.
    """
    # Calculate total integral
    total_integral = np.sum(volume_data_sh**2) * dV if square else np.sum(np.abs(volume_data_sh)) * dV

    out = []

    for target_percentage in values:
        target_integral = total_integral * abs(target_percentage)

        # Binary search for isovalue
        min_val = np.min(np.abs(volume_data_sh))
        max_val = np.max(np.abs(volume_data_sh))

        for _ in range(max_iterations):
            current_iso = (min_val + max_val) / 2
            mask = np.abs(volume_data_sh) >= current_iso
            current_integral = np.sum(volume_data_sh[mask]**2) * dV if square else np.sum(np.abs(volume_data_sh[mask])) * dV

            if abs(current_integral - target_integral) < tolerance * total_integral:
                break

            if current_integral > target_integral:
                min_val = current_iso
            else:
                max_val = current_iso

        # If we hit max iterations, return best estimate
        out.append(np.copysign(current_iso, target_percentage))

    isovals = [ round_sigfigs(x, 2) for x in out ]
    return isovals
