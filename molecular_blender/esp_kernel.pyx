# distutils: language = c++
# cython: language_level=3, boundscheck=False, wraparound=False, cdivision=True
"""Coulomb integrals of contracted Gaussian densities at surface points, s through g."""
from libc.math cimport sqrt
import numpy as np


def coulomb_sum(double[:, ::1] points, int degree, double[::1] p,
                double[:, ::1] centers, double[:, ::1] coefficients,
                int[:, ::1] rules, int[::1] counts, double[:, ::1] table):
    cdef Py_ssize_t k, pair
    cdef int n, h, j, grid
    cdef double delta[3]
    cdef double boys[9]
    cdef double aux[165]  # C(8+3, 3) derivatives for a pair of g shells
    cdef double t, dx, value, scale, contribution
    cdef double[::1] result = np.zeros(points.shape[0], dtype=np.float64)
    if not 0 <= degree <= 8 or points.shape[1] != 3 or centers.shape[1] != 3:
        raise ValueError('Invalid Gaussian integral dimensions')
    if (p.shape[0] != centers.shape[0] or p.shape[0] != coefficients.shape[0]
            or counts.shape[0] != degree+1 or rules.shape[1] != 4
            or counts[degree] != coefficients.shape[1] or rules.shape[0] != counts[degree]
            or rules.shape[0] > 165 or table.shape[0] < 641 or table.shape[1] < degree+7):
        raise ValueError('Inconsistent Gaussian integral workspace')
    for n in range(degree+1):
        if counts[n] != (n+1)*(n+2)*(n+3)//6:
            raise ValueError('Invalid Gaussian derivative count')
    for j in range(1, counts[degree]):
        if (rules[j, 0] < 0 or rules[j, 0] > 2 or rules[j, 1] < 0 or rules[j, 1] >= j
                or rules[j, 2] < 0 or rules[j, 2] >= j):
            raise ValueError('Invalid Gaussian derivative indices')
    for pair in range(p.shape[0]):
        for k in range(points.shape[0]):
            for j in range(3):
                delta[j] = centers[pair, j] - points[k, j]
            t = p[pair]*(delta[0]*delta[0] + delta[1]*delta[1] + delta[2]*delta[2])
            if t >= 40.:
                boys[0] = .886226925452758 / sqrt(t)
                for n in range(1, degree+1):
                    boys[n] = boys[n-1] * ((2*n-1)*.5) / t
            elif t >= 0.:
                grid = <int>(t*16+.5)
                dx = grid*.0625 - t
                for n in range(degree+1):
                    value = table[grid, n+6]
                    for h in range(6, 0, -1):
                        value = table[grid, n+h-1] + dx*value/h
                    boys[n] = value
            else:
                raise ValueError('Gaussian coordinates and exponents must be finite')
            scale = 1.
            for n in range(degree+1):
                boys[n] *= scale
                scale *= -2*p[pair]
            for n in range(degree, -1, -1):
                for j in range(counts[degree-n]-1, 0, -1):
                    aux[j] = delta[rules[j, 0]]*aux[rules[j, 1]] + rules[j, 3]*aux[rules[j, 2]]
                aux[0] = boys[n]
            contribution = 0.
            for j in range(counts[degree]):
                contribution += coefficients[pair, j]*aux[j]
            result[k] += contribution
    return np.asarray(result)
