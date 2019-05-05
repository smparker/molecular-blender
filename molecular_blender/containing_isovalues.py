# -*- coding: utf-8 -*-

import numpy as np
import math
from .util import round_sigfigs

def isovalue_containing_proportion(values, volume_data, dV):
    cumulative = 0.0
    data = np.zeros([len(volume_data), 2], dtype=np.float32)
    for i, v in enumerate(sorted(volume_data, key=lambda x: x*x, reverse=True)):
        data[i,0] = abs(v)
        data[i,1] = cumulative
        cumulative += v*v*dV

    # rescale containing values according to integrated density
    # equivalent to uniformly scaling cumulative data so norm is 1
    values = [ cumulative * v for v in values ]
    isovals = [ None for x in values ]

    for i in range(1, data.shape[0]):
        x1, y1 = data[i-1,:]
        x2, y2 = data[i,:]

        for ival, val in enumerate(values):
            if y2 >= abs(val) and y1 < abs(val):
                xval = x1 + (abs(val) - y1)/(y2 - y1)*(x2 - x1)
                isovals[ival] = math.copysign(xval, val)

    # round isovalues to 2 decimal places
    isovals = [ round_sigfigs(x, 2) for x in isovals ]

    return isovals
