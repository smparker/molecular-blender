#distutils: language = c++
#cython: cdivision = True
#cython: profile = True
# -*- coding: utf-8 -*-
#
#  Molecular Blender
#  Filename: compute_ao_cy.py
#  Copyright (C) 2017-2025 Shane Parker
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

"""Classes and functions to compute orbital values in real space."""

cimport cython

import math
import numpy as np
cimport numpy as np

from numpy cimport ndarray
from libc.math cimport sqrt, exp, fabs, ceil

DTYPE = np.float32
ctypedef np.float32_t DTYPE_t

DEF CYDEBUG = False
# The AO values and their density-matrix product coexist during contraction.
_DENSITY_WORKSPACE_BYTES = 64 * 1024**2

from .constants import ang2bohr
from .compute_ao_py import density_shell_weights, shell_screening_parameters
from .containing_isovalues import isovalue_containing_proportion

molden_shorder = [[""],
           ["x", "y", "z"],
           ["xx", "yy", "zz", "xy", "xz", "yz"],
           ["xxx", "yyy", "zzz", "xyy", "xxy",
            "xxz", "xzz", "yzz", "yyz", "xyz"],
           ["xxxx", "yyyy", "zzzz", "xxxy", "xxxz", "yyyx", "yyyz", "zzzx", "zzzy",
            "xxyy", "xxzz", "yyzz", "xxyz", "yyxz", "zzxy"]
          ]

cdef gamma2n(i):
    """Returns Gamma(n/2), only for half-integers"""
    return np.array([gamma2n_impl(int(round(x))) for x in i])

cdef gamma2n_impl(i):
    """Returns gamma2n_impl(n) where n is input as 2*n=i"""
    if i <= 0:
        raise Exception("gamma2n_impl function cannot accept negative numbers")
    elif i == 2:
        return 1
    elif i == 1:
        return sqrt(np.pi)
    else:
        return (0.5 * i - 1.) * gamma2n_impl(i - 2)

ctypedef struct Interval:
    float low, high

cdef inline Interval axis_bounds(const DTYPE_t* axis, Py_ssize_t count) noexcept nogil:
    cdef Interval bounds
    cdef Py_ssize_t i
    bounds.low = bounds.high = axis[0]
    for i in range(1, count):
        bounds.low = min(bounds.low, axis[i])
        bounds.high = max(bounds.high, axis[i])
    return bounds

cdef inline float distance_to_interval(float x, Interval bounds) noexcept nogil:
    return max(0., bounds.low - x, x - bounds.high)

cdef class Shell:
    """ Fast Cython container for shell data"""
    cdef readonly float X, Y, Z
    cdef readonly int l, start, size
    cdef readonly const DTYPE_t[::1] center
    cdef readonly object exponents
    cdef public object prim_coeffs
    cdef public object norms
    cdef readonly object denormed_coeffs
    cdef public float mxx
    cdef readonly float most_diffuse
    cdef double log_envelope
    cdef readonly float screen_log_weight
    cdef readonly float logthr
    cdef public float mxnorm

    cdef const DTYPE_t[::1] exponent_view, radial_view

    def __init__(self, center, int l,
            exponents, coeff,
            int start, float thr=1e-6, float mxx=1.0,
            bint denorm=True):
        if self.exponents is not None:
            raise ValueError("A native shell is already initialized")
        if not 0 <= l <= 4 or start < 0:
            raise ValueError("Expected an s through g shell and nonnegative AO offset")
        center = np.array(center, dtype=DTYPE)
        if center.shape != (3,) or not np.isfinite(center).all():
            raise ValueError("Expected three finite shell coordinates")
        center.flags.writeable = False
        self.center = center

        self.X = center[0]
        self.Y = center[1]
        self.Z = center[2]

        self.l = int(l)

        self.start = start
        self.size = (self.l + 1) * (self.l + 2) // 2

        self.exponents = np.array(exponents, dtype=DTYPE)
        self.prim_coeffs = np.array(coeff, dtype=DTYPE)

        if (self.exponents.ndim != 1 or not self.exponents.size or
                self.prim_coeffs.shape != self.exponents.shape or
                not np.isfinite(self.exponents).all() or np.any(self.exponents <= 0) or
                not np.isfinite(self.prim_coeffs).all() or not np.isfinite(thr) or thr <= 0):
            raise ValueError("Expected positive finite exponents and threshold, and one finite coefficient per primitive")

        self.norms = 1.0 / np.sqrt(self.__norms(), dtype=DTYPE)
        if not denorm:
            self.norms[:,:] = 1.0
        sh_norms = 1.0 / np.sqrt(self.__shellnorms(), dtype=DTYPE) if denorm else 1.0
        self.denormed_coeffs = np.array(coeff, dtype=DTYPE) * sh_norms

        # most diffuse exponent in shell
        cdef int imin = np.argmin(self.exponents)
        self.mxx = mxx  # save max x just for use later
        self.most_diffuse = np.min(self.exponents)
        self.mxnorm = np.max(self.norms[:, imin]) * abs(self.prim_coeffs[imin])
        self.log_envelope, self.screen_log_weight = shell_screening_parameters(l, self.most_diffuse, self.denormed_coeffs)
        self.logthr = self.screen_log_weight * (self.log_envelope - math.log(thr))

        # Pin primitive storage for point, plane and grid queries.
        self.exponents.flags.writeable = False
        self.denormed_coeffs.flags.writeable = False
        self.exponent_view = self.exponents
        self.radial_view = self.denormed_coeffs

    def __norms(self):
        """Returns [nL, nexp] array of norms"""
        out = np.zeros([self.size, len(self.exponents)], dtype=DTYPE)
        for ixyz, xyz in enumerate(molden_shorder[self.l]):
            cxyz = [xyz.count("x"), xyz.count("y"), xyz.count("z")]
            for ig, ex in enumerate(self.exponents):
                out[ixyz, ig] = gaussian_overlap(
                    self.center, ex, cxyz, self.center, ex, cxyz)

        return out

    def __shellnorms(self):
        """
        Returns [nexp] array of shell norms, i.e., norms of Gaussians with all angular
        in the x-direction. For example, for L=2, norms of exp(-zeta r^2) x x.
        """
        out = np.zeros_like(self.exponents, dtype=DTYPE)
        for ig, ex in enumerate(self.exponents):
            out[ig] = gaussian_overlap(self.center, ex, [ self.l, 0, 0 ], self.center, ex, [ self.l, 0, 0 ])
        return out

    def bounding_box_size(self, thr, logmxcoeff=0.0):
        """Returns half the edgelength of a box outside of which the shell
        is guaranteed to be below the given threshold"""
        rr = self.screen_log_weight * (self.log_envelope + logmxcoeff - math.log(thr))
        return math.sqrt(max(0.0, rr))

    @cython.boundscheck(CYDEBUG)
    @cython.wraparound(CYDEBUG)
    @cython.initializedcheck(CYDEBUG)
    @cython.nonecheck(CYDEBUG)
    cpdef float value(self, DTYPE_t x, DTYPE_t y, DTYPE_t z, ndarray[DTYPE_t, ndim=1] coeff,
            float logmxcoeff=0.):

        cdef float dx = x - self.X
        cdef float dy = y - self.Y
        cdef float dz = z - self.Z

        cdef float sqrt3 = sqrt(3.0)
        cdef float sqrt5 = sqrt(5.0)
        cdef float sqrt15 = sqrt(15.0)
        cdef float sqrt_35_5 = sqrt(35.0/5.0)
        cdef float sqrt_35_3 = sqrt(35.0/3.0)
        cdef float sqrt_35_1 = sqrt(35.0/1.0)

        cdef float logthresh = self.logthr + self.screen_log_weight * logmxcoeff

        cdef Py_ssize_t nexp = self.exponent_view.shape[0]
        cdef int l = self.l
        cdef int shellsize = self.size

        cdef const DTYPE_t[::1] exponents = self.exponent_view
        cdef const DTYPE_t[::1] denormed_coeffs = self.radial_view

        if coeff is None or coeff.shape[0] < shellsize:
            raise ValueError("Too few coefficients for this shell")
        cdef float rr = dx*dx + dy*dy + dz*dz
        if rr >= logthresh:
            return 0.0

        cdef float radial = 0.0

        for iexp in range(nexp):
            radial += denormed_coeffs[iexp] * exp(-exponents[iexp]*rr)

        cdef float out = 0.0

        if l == 0:
            out += radial * coeff[0]
        elif l == 1:
            out += radial * (
                      coeff[0] * dx
                    + coeff[1] * dy
                    + coeff[2] * dz)
        elif l == 2:
            out += radial * (
                      coeff[0] * dx*dx
                    + coeff[1] * dy*dy
                    + coeff[2] * dz*dz
                    + coeff[3] * dx*dy * sqrt3
                    + coeff[4] * dx*dz * sqrt3
                    + coeff[5] * dy*dz * sqrt3
                    )
        elif l == 3:
            out += radial * (
                      coeff[0] * dx*dx*dx
                    + coeff[1] * dy*dy*dy
                    + coeff[2] * dz*dz*dz
                    + coeff[3] * dx*dy*dy * sqrt5
                    + coeff[4] * dx*dx*dy * sqrt5
                    + coeff[5] * dx*dx*dz * sqrt5
                    + coeff[6] * dx*dz*dz * sqrt5
                    + coeff[7] * dy*dz*dz * sqrt5
                    + coeff[8] * dy*dy*dz * sqrt5
                    + coeff[9] * dx*dy*dz * sqrt15
                    )
        elif l == 4:
            out += radial * (
                      coeff[0] * dx*dx*dx*dx
                    + coeff[1] * dy*dy*dy*dy
                    + coeff[2] * dz*dz*dz*dz
                    + coeff[3] * dx*dx*dx*dy * sqrt_35_5
                    + coeff[4] * dx*dx*dx*dz * sqrt_35_5
                    + coeff[5] * dy*dy*dy*dx * sqrt_35_5
                    + coeff[6] * dy*dy*dy*dz * sqrt_35_5
                    + coeff[7] * dz*dz*dz*dx * sqrt_35_5
                    + coeff[8] * dz*dz*dz*dy * sqrt_35_5
                    + coeff[9] * dx*dx*dy*dy * sqrt_35_3
                    + coeff[10] * dx*dx*dz*dz * sqrt_35_3
                    + coeff[11] * dy*dy*dz*dz * sqrt_35_3
                    + coeff[12] * dx*dx*dy*dz * sqrt_35_1
                    + coeff[13] * dy*dy*dx*dz * sqrt_35_1
                    + coeff[14] * dz*dz*dx*dy * sqrt_35_1
                    )
        return out

    @cython.boundscheck(CYDEBUG)
    @cython.wraparound(CYDEBUG)
    @cython.initializedcheck(CYDEBUG)
    @cython.nonecheck(CYDEBUG)
    cpdef void add_plane_values(self, const DTYPE_t[::1] xx,
            const DTYPE_t[::1] yy, float z_a,
            const DTYPE_t[::1] coeff,
            DTYPE_t[:,::1] target, float logmxcoeff=0.0):

        cdef float X = self.X
        cdef float Y = self.Y
        cdef float Z = self.Z

        cdef float logthresh = self.logthr + self.screen_log_weight * logmxcoeff

        cdef Py_ssize_t nx = len(xx)
        cdef Py_ssize_t ny = len(yy)

        cdef Py_ssize_t nexp = self.exponent_view.shape[0]

        cdef int l = self.l
        cdef int shellsize = self.size

        if coeff.shape[0] < shellsize or target.shape[0] != nx or target.shape[1] != ny:
            raise ValueError("Plane buffer dimensions do not match the shell and grid")
        if nx == 0 or ny == 0:
            return
        cdef float min_x = distance_to_interval(X, axis_bounds(&xx[0], nx))
        cdef float min_y = distance_to_interval(Y, axis_bounds(&yy[0], ny))
        cdef float min_z = fabs(z_a - Z)

        cdef float min_rr = min_x*min_x + min_y*min_y + min_z*min_z
        if min_rr >= logthresh:
            return

        cdef float sqrt3 = sqrt(3.0)
        cdef float sqrt5 = sqrt(5.0)
        cdef float sqrt15 = sqrt(15.0)
        cdef float sqrt_35_5 = sqrt(35.0/5.0)
        cdef float sqrt_35_3 = sqrt(35.0/3.0)
        cdef float sqrt_35_1 = sqrt(35.0/1.0)

        cdef float dx, dy, dz, radial

        cdef const DTYPE_t[::1] denormed_coeffs = self.radial_view
        cdef const DTYPE_t[::1] exponents = self.exponent_view

        cdef ndarray[DTYPE_t, ndim=2] expyy = np.zeros([ny, nexp], dtype=DTYPE)
        cdef ndarray[DTYPE_t, ndim=1] expxx = np.zeros(nexp, dtype=DTYPE)

        for iy in range(ny):
            dy = yy[iy] - Y
            for iexp in range(nexp):
                expyy[iy, iexp] = exp(-exponents[iexp]*dy*dy)

        dz = z_a - Z
        for ix in range(nx):
            dx = xx[ix] - X
            for iexp in range(nexp):
                expxx[iexp] = exp(-exponents[iexp] * (dx*dx + dz*dz)) * denormed_coeffs[iexp]

            for iy in range(ny):
                dy = yy[iy] - Y

                if (dx*dx + dy*dy + dz*dz >= logthresh):
                    continue

                radial = 0.0
                for iexp in range(nexp):
                    radial += expxx[iexp] * expyy[iy, iexp]

                if l == 0:
                    target[ix,iy] += radial * coeff[0]
                elif l == 1:
                    target[ix,iy] += radial * (
                              coeff[0] * dx
                            + coeff[1] * dy
                            + coeff[2] * dz)
                elif l == 2:
                    target[ix,iy] += radial * (
                              coeff[0] * dx*dx
                            + coeff[1] * dy*dy
                            + coeff[2] * dz*dz
                            + coeff[3] * dx*dy * sqrt3
                            + coeff[4] * dx*dz * sqrt3
                            + coeff[5] * dy*dz * sqrt3
                            )
                elif l == 3:
                    target[ix,iy] += radial * (
                              coeff[0] * dx*dx*dx
                            + coeff[1] * dy*dy*dy
                            + coeff[2] * dz*dz*dz
                            + coeff[3] * dx*dy*dy * sqrt5
                            + coeff[4] * dx*dx*dy * sqrt5
                            + coeff[5] * dx*dx*dz * sqrt5
                            + coeff[6] * dx*dz*dz * sqrt5
                            + coeff[7] * dy*dz*dz * sqrt5
                            + coeff[8] * dy*dy*dz * sqrt5
                            + coeff[9] * dx*dy*dz * sqrt15
                            )
                elif l == 4:
                    target[ix,iy] += radial * (
                              coeff[0] * dx*dx*dx*dx
                            + coeff[1] * dy*dy*dy*dy
                            + coeff[2] * dz*dz*dz*dz
                            + coeff[3] * dx*dx*dx*dy * sqrt_35_5
                            + coeff[4] * dx*dx*dx*dz * sqrt_35_5
                            + coeff[5] * dy*dy*dy*dx * sqrt_35_5
                            + coeff[6] * dy*dy*dy*dz * sqrt_35_5
                            + coeff[7] * dz*dz*dz*dx * sqrt_35_5
                            + coeff[8] * dz*dz*dz*dy * sqrt_35_5
                            + coeff[9] * dx*dx*dy*dy * sqrt_35_3
                            + coeff[10] * dx*dx*dz*dz * sqrt_35_3
                            + coeff[11] * dy*dy*dz*dz * sqrt_35_3
                            + coeff[12] * dx*dx*dy*dz * sqrt_35_1
                            + coeff[13] * dy*dy*dx*dz * sqrt_35_1
                            + coeff[14] * dz*dz*dx*dy * sqrt_35_1
                            )

    @cython.boundscheck(CYDEBUG)
    @cython.wraparound(CYDEBUG)
    @cython.initializedcheck(CYDEBUG)
    @cython.nonecheck(CYDEBUG)
    cpdef bint compute_chi_plane(self, const DTYPE_t[::1] xx,
            const DTYPE_t[::1] yy, float z_a,
            DTYPE_t[:,:,::1] target, float logmxcoeff=0.0):

        cdef float x = self.X
        cdef float y = self.Y
        cdef float z = self.Z

        cdef float logthresh = self.logthr + self.screen_log_weight * logmxcoeff

        cdef Py_ssize_t nx = len(xx)
        cdef Py_ssize_t ny = len(yy)

        cdef Py_ssize_t nexp = self.exponent_view.shape[0]

        cdef int l = self.l
        cdef int shellsize = self.size

        if target.shape[0] < shellsize or target.shape[1] != nx or target.shape[2] != ny:
            raise ValueError("Plane buffer dimensions do not match the shell and grid")
        if nx == 0 or ny == 0:
            return False
        cdef float min_x = distance_to_interval(x, axis_bounds(&xx[0], nx))
        cdef float min_y = distance_to_interval(y, axis_bounds(&yy[0], ny))
        cdef float min_z = fabs(z_a - z)

        cdef float min_rr = min_x*min_x + min_y*min_y + min_z*min_z
        if min_rr >= logthresh:
            return False

        cdef float sqrt3 = sqrt(3.0)
        cdef float sqrt5 = sqrt(5.0)
        cdef float sqrt15 = sqrt(15.0)
        cdef float sqrt_35_5 = sqrt(35.0/5.0)
        cdef float sqrt_35_3 = sqrt(35.0/3.0)
        cdef float sqrt_35_1 = sqrt(35.0/1.0)

        cdef float dx, dy, dz, radial

        cdef const DTYPE_t[::1] denormed_coeffs = self.radial_view
        cdef const DTYPE_t[::1] exponents = self.exponent_view

        cdef ndarray[DTYPE_t, ndim=2] expyy = np.zeros([ny, nexp], dtype=DTYPE)
        cdef ndarray[DTYPE_t, ndim=1] expxx = np.zeros(nexp, dtype=DTYPE)

        for iy in range(ny):
            dy = yy[iy] - y
            for iexp in range(nexp):
                expyy[iy, iexp] = exp(-exponents[iexp]*dy*dy)

        dz = z_a - z

        cdef int icoeff = 0
        for ix in range(nx):
            dx = xx[ix] - x
            for iexp in range(nexp):
                expxx[iexp] = exp(-exponents[iexp] * (dx*dx + dz*dz)) * denormed_coeffs[iexp]

            for iy in range(ny):
                dy = yy[iy] - y

                if (dx*dx + dy*dy + dz*dz >= logthresh):
                    continue

                radial = 0.0
                for iexp in range(nexp):
                    radial += expxx[iexp] * expyy[iy, iexp]

                if l == 0:
                    target[icoeff,ix,iy] = radial
                elif l == 1:
                    target[icoeff,ix,iy] = radial * dx
                    target[icoeff+1,ix,iy] = radial * dy
                    target[icoeff+2,ix,iy] = radial * dz
                elif l == 2:
                    target[icoeff,ix,iy] = radial * dx*dx
                    target[icoeff+1,ix,iy] = radial * dy*dy
                    target[icoeff+2,ix,iy] = radial * dz*dz
                    target[icoeff+3,ix,iy] = radial * dx*dy * sqrt3
                    target[icoeff+4,ix,iy] = radial * dx*dz * sqrt3
                    target[icoeff+5,ix,iy] = radial * dy*dz * sqrt3
                elif l == 3:
                    target[icoeff,ix,iy] = radial * dx*dx*dx
                    target[icoeff+1,ix,iy] = radial * dy*dy*dy
                    target[icoeff+2,ix,iy] = radial * dz*dz*dz
                    target[icoeff+3,ix,iy] = radial * dx*dy*dy * sqrt5
                    target[icoeff+4,ix,iy] = radial * dx*dx*dy * sqrt5
                    target[icoeff+5,ix,iy] = radial * dx*dx*dz * sqrt5
                    target[icoeff+6,ix,iy] = radial * dx*dz*dz * sqrt5
                    target[icoeff+7,ix,iy] = radial * dy*dz*dz * sqrt5
                    target[icoeff+8,ix,iy] = radial * dy*dy*dz * sqrt5
                    target[icoeff+9,ix,iy] = radial * dx*dy*dz * sqrt15
                elif l == 4:
                    target[icoeff,ix,iy] = radial * dx*dx*dx*dx
                    target[icoeff+1,ix,iy] = radial * dy*dy*dy*dy
                    target[icoeff+2,ix,iy] = radial * dz*dz*dz*dz
                    target[icoeff+3,ix,iy] = radial * dx*dx*dx*dy * sqrt_35_5
                    target[icoeff+4,ix,iy] = radial * dx*dx*dx*dz * sqrt_35_5
                    target[icoeff+5,ix,iy] = radial * dy*dy*dy*dx * sqrt_35_5
                    target[icoeff+6,ix,iy] = radial * dy*dy*dy*dz * sqrt_35_5
                    target[icoeff+7,ix,iy] = radial * dz*dz*dz*dx * sqrt_35_5
                    target[icoeff+8,ix,iy] = radial * dz*dz*dz*dy * sqrt_35_5
                    target[icoeff+9,ix,iy] = radial * dx*dx*dy*dy * sqrt_35_3
                    target[icoeff+10,ix,iy] = radial * dx*dx*dz*dz * sqrt_35_3
                    target[icoeff+11,ix,iy] = radial * dy*dy*dz*dz * sqrt_35_3
                    target[icoeff+12,ix,iy] = radial * dx*dx*dy*dz * sqrt_35_1
                    target[icoeff+13,ix,iy] = radial * dy*dy*dx*dz * sqrt_35_1
                    target[icoeff+14,ix,iy] = radial * dz*dz*dx*dy * sqrt_35_1
        return True


def gaussian_overlap(axyz, aexp, apoly, bxyz, bexp, bpoly):
    """returns overlap integral for gaussians centered at xyz with exponents exp
    and angular momentum poly"""
    abpoly = [k + n for k, n in zip(apoly, bpoly)]
    if any([i % 2 == 1 for i in abpoly]):
        return 0.0
    abpoly = np.array(abpoly, dtype=DTYPE) + 1.
    fac = np.power(aexp + bexp, -0.5 * abpoly, dtype=DTYPE)
    gam = gamma2n(abpoly)
    return np.prod(fac * gam)


@cython.boundscheck(CYDEBUG)
@cython.wraparound(CYDEBUG)
@cython.initializedcheck(CYDEBUG)
@cython.nonecheck(CYDEBUG)
cdef void add_box_values(Shell shell,
        Py_ssize_t nx, const DTYPE_t* xx,
        Py_ssize_t ny, const DTYPE_t* yy,
        Py_ssize_t nz, const DTYPE_t* zz,
        DTYPE_t* coeff,
        DTYPE_t* target,
        float logmxcoeff, const Interval* bounds, DTYPE_t* scratch) noexcept nogil:
    cdef float X = shell.X
    cdef float Y = shell.Y
    cdef float Z = shell.Z

    cdef float logthresh = shell.logthr + shell.screen_log_weight * logmxcoeff

    cdef Py_ssize_t ix
    cdef Py_ssize_t iy
    cdef Py_ssize_t iz
    cdef Py_ssize_t iexp

    cdef Py_ssize_t nexp = shell.exponent_view.shape[0]

    cdef int l = shell.l

    cdef float min_x = distance_to_interval(X, bounds[0])
    cdef float min_y = distance_to_interval(Y, bounds[1])
    cdef float min_z = distance_to_interval(Z, bounds[2])
    cdef float min_rr = min_x*min_x + min_y*min_y + min_z*min_z
    if min_rr >= logthresh:
        return

    cdef float sqrt3 = sqrt(3.0)
    cdef float sqrt5 = sqrt(5.0)
    cdef float sqrt15 = sqrt(15.0)
    cdef float sqrt_35_5 = sqrt(35.0/5.0)
    cdef float sqrt_35_3 = sqrt(35.0/3.0)
    cdef float sqrt_35_1 = sqrt(35.0/1.0)

    cdef float dx, dy, dz, radial, radial_xy, angular_xy, angular_z
    cdef Py_ssize_t ixyz

    cdef const DTYPE_t* denormed_coeffs = &shell.radial_view[0]
    cdef const DTYPE_t* exponents = &shell.exponent_view[0]

    cdef DTYPE_t* expxx = scratch
    cdef DTYPE_t* expyy = expxx + nexp
    cdef DTYPE_t* expzz = expyy + ny*nexp

    cdef float xthresh, ythresh, zthresh

    for iy in range(ny):
        dy = yy[iy] - Y
        for iexp in range(nexp):
            expyy[iy*nexp + iexp] = exp(-exponents[iexp]*dy*dy)

    for iz in range(nz):
        dz = zz[iz] - Z
        for iexp in range(nexp):
            expzz[iz*nexp + iexp] = exp(-exponents[iexp]*dz*dz)

    for ix in range(nx):
        dx = xx[ix] - X

        xthresh = dx * dx + min_y*min_y + min_z*min_z
        if xthresh > logthresh:
            continue

        for iexp in range(nexp):
            expxx[iexp] = exp(-exponents[iexp] * dx*dx) * denormed_coeffs[iexp]

        for iy in range(ny):
            dy = yy[iy] - Y

            ythresh = dx * dx + dy * dy + min_z*min_z
            if ythresh > logthresh:
                continue

            # Most outer s and p shells have one primitive. Factor the x/y terms
            # once per row so the inner loop walks contiguous z samples.
            if nexp == 1 and l <= 1:
                radial_xy = expxx[0] * expyy[iy]
                angular_xy = coeff[0] if l == 0 else coeff[0] * dx + coeff[1] * dy
                angular_z = 0 if l == 0 else coeff[2]
                for iz in range(nz):
                    dz = zz[iz] - Z
                    if dx * dx + dy * dy + dz * dz < logthresh:
                        target[(ix*ny + iy)*nz + iz] += (radial_xy * expzz[iz]) * (angular_xy + angular_z * dz)
                continue

            for iz in range(nz):
                dz = zz[iz] - Z

                zthresh = dx * dx + dy * dy + dz*dz
                if zthresh >= logthresh:
                    continue

                ixyz = iz + iy*nz + ix*nz*ny

                radial = 0.0
                for iexp in range(nexp):
                    radial += expxx[iexp] * expyy[iy*nexp + iexp] * expzz[iz*nexp + iexp]

                if l == 0:
                    target[ixyz] += radial * coeff[0]
                elif l == 1:
                    target[ixyz] += radial * (
                              coeff[0] * dx
                            + coeff[1] * dy
                            + coeff[2] * dz)
                elif l == 2:
                    target[ixyz] += radial * (
                              coeff[0] * dx*dx
                            + coeff[1] * dy*dy
                            + coeff[2] * dz*dz
                            + coeff[3] * dx*dy * sqrt3
                            + coeff[4] * dx*dz * sqrt3
                            + coeff[5] * dy*dz * sqrt3
                            )
                elif l == 3:
                    target[ixyz] += radial * (
                              coeff[0] * dx*dx*dx
                            + coeff[1] * dy*dy*dy
                            + coeff[2] * dz*dz*dz
                            + coeff[3] * dx*dy*dy * sqrt5
                            + coeff[4] * dx*dx*dy * sqrt5
                            + coeff[5] * dx*dx*dz * sqrt5
                            + coeff[6] * dx*dz*dz * sqrt5
                            + coeff[7] * dy*dz*dz * sqrt5
                            + coeff[8] * dy*dy*dz * sqrt5
                            + coeff[9] * dx*dy*dz * sqrt15
                            )
                elif l == 4:
                    target[ixyz] += radial * (
                              coeff[0] * dx*dx*dx*dx
                            + coeff[1] * dy*dy*dy*dy
                            + coeff[2] * dz*dz*dz*dz
                            + coeff[3] * dx*dx*dx*dy * sqrt_35_5
                            + coeff[4] * dx*dx*dx*dz * sqrt_35_5
                            + coeff[5] * dy*dy*dy*dx * sqrt_35_5
                            + coeff[6] * dy*dy*dy*dz * sqrt_35_5
                            + coeff[7] * dz*dz*dz*dx * sqrt_35_5
                            + coeff[8] * dz*dz*dz*dy * sqrt_35_5
                            + coeff[9] * dx*dx*dy*dy * sqrt_35_3
                            + coeff[10] * dx*dx*dz*dz * sqrt_35_3
                            + coeff[11] * dy*dy*dz*dz * sqrt_35_3
                            + coeff[12] * dx*dx*dy*dz * sqrt_35_1
                            + coeff[13] * dy*dy*dx*dz * sqrt_35_1
                            + coeff[14] * dz*dz*dx*dy * sqrt_35_1
                            )

@cython.boundscheck(CYDEBUG)
@cython.wraparound(CYDEBUG)
@cython.initializedcheck(CYDEBUG)
@cython.nonecheck(CYDEBUG)
cdef void compute_chi_box_shell(Shell shell,
        Py_ssize_t nx, const DTYPE_t* xx,
        Py_ssize_t ny, const DTYPE_t* yy,
        Py_ssize_t nz, const DTYPE_t* zz,
        DTYPE_t* target,
        Py_ssize_t sx, Py_ssize_t sy, Py_ssize_t sz,
        float logmxcoeff, const Interval* bounds, DTYPE_t* scratch) noexcept nogil:
    """Computes the value of the shell on a grid

    target is a 3D array of shape (nx, ny, nz, shellsize)
    sx, sy, sz are strides, such that target[ix,iy,iz,ishell] is at
    target[ix*sx + iy*sy + iz*sz + ishell]
    """
    cdef float X = shell.X
    cdef float Y = shell.Y
    cdef float Z = shell.Z

    cdef float logthresh = shell.logthr + shell.screen_log_weight * logmxcoeff

    cdef Py_ssize_t ix
    cdef Py_ssize_t iy
    cdef Py_ssize_t iz
    cdef Py_ssize_t iexp
    cdef Py_ssize_t target_ixyz

    cdef Py_ssize_t nexp = shell.exponent_view.shape[0]

    cdef int l = shell.l

    cdef float min_x = distance_to_interval(X, bounds[0])
    cdef float min_y = distance_to_interval(Y, bounds[1])
    cdef float min_z = distance_to_interval(Z, bounds[2])
    cdef float min_rr = min_x*min_x + min_y*min_y + min_z*min_z
    if min_rr >= logthresh:
        return

    cdef float sqrt3 = sqrt(3.0)
    cdef float sqrt5 = sqrt(5.0)
    cdef float sqrt15 = sqrt(15.0)
    cdef float sqrt_35_5 = sqrt(35.0/5.0)
    cdef float sqrt_35_3 = sqrt(35.0/3.0)
    cdef float sqrt_35_1 = sqrt(35.0/1.0)

    cdef float dx, dy, dz, radial

    cdef const DTYPE_t* denormed_coeffs = &shell.radial_view[0]
    cdef const DTYPE_t* exponents = &shell.exponent_view[0]

    cdef DTYPE_t* expxx = scratch
    cdef DTYPE_t* expyy = expxx + nexp
    cdef DTYPE_t* expzz = expyy + ny*nexp

    cdef float xthresh, ythresh, zthresh


    for iy in range(ny):
        dy = yy[iy] - Y
        for iexp in range(nexp):
            expyy[iy*nexp + iexp] = exp(-exponents[iexp]*dy*dy)

    for iz in range(nz):
        dz = zz[iz] - Z
        for iexp in range(nexp):
            expzz[iz*nexp + iexp] = exp(-exponents[iexp]*dz*dz)

    for ix in range(nx):
        dx = xx[ix] - X

        xthresh = dx * dx + min_y*min_y + min_z*min_z
        if xthresh > logthresh:
            continue

        for iexp in range(nexp):
            expxx[iexp] = exp(-exponents[iexp] * dx*dx) * denormed_coeffs[iexp]

        for iy in range(ny):
            dy = yy[iy] - Y

            ythresh = dx * dx + dy * dy + min_z*min_z
            if ythresh > logthresh:
                continue

            for iz in range(nz):
                dz = zz[iz] - Z

                zthresh = dx * dx + dy * dy + dz*dz
                if zthresh >= logthresh:
                    continue

                radial = 0.0
                for iexp in range(nexp):
                    radial += expxx[iexp] * expyy[iy*nexp + iexp] * expzz[iz*nexp + iexp]

                target_ixyz = ix*sx + iy*sy + iz*sz

                if l == 0:
                    target[target_ixyz] = radial
                elif l == 1:
                    target[target_ixyz] = radial * dx
                    target[target_ixyz+1] = radial * dy
                    target[target_ixyz+2] = radial * dz
                elif l == 2:
                    target[target_ixyz] = radial * dx*dx
                    target[target_ixyz+1] = radial * dy*dy
                    target[target_ixyz+2] = radial * dz*dz
                    target[target_ixyz+3] = radial * dx*dy * sqrt3
                    target[target_ixyz+4] = radial * dx*dz * sqrt3
                    target[target_ixyz+5] = radial * dy*dz * sqrt3
                elif l == 3:
                    target[target_ixyz] = radial * dx*dx*dx
                    target[target_ixyz+1] = radial * dy*dy*dy
                    target[target_ixyz+2] = radial * dz*dz*dz
                    target[target_ixyz+3] = radial * dx*dy*dy * sqrt5
                    target[target_ixyz+4] = radial * dx*dx*dy * sqrt5
                    target[target_ixyz+5] = radial * dx*dx*dz * sqrt5
                    target[target_ixyz+6] = radial * dx*dz*dz * sqrt5
                    target[target_ixyz+7] = radial * dy*dz*dz * sqrt5
                    target[target_ixyz+8] = radial * dy*dy*dz * sqrt5
                    target[target_ixyz+9] = radial * dx*dy*dz * sqrt15
                elif l == 4:
                    target[target_ixyz] = radial * dx*dx*dx*dx
                    target[target_ixyz+1] = radial * dy*dy*dy*dy
                    target[target_ixyz+2] = radial * dz*dz*dz*dz
                    target[target_ixyz+3] = radial * dx*dx*dx*dy * sqrt_35_5
                    target[target_ixyz+4] = radial * dx*dx*dx*dz * sqrt_35_5
                    target[target_ixyz+5] = radial * dy*dy*dy*dx * sqrt_35_5
                    target[target_ixyz+6] = radial * dy*dy*dy*dz * sqrt_35_5
                    target[target_ixyz+7] = radial * dz*dz*dz*dx * sqrt_35_5
                    target[target_ixyz+8] = radial * dz*dz*dz*dy * sqrt_35_5
                    target[target_ixyz+9] = radial * dx*dx*dy*dy * sqrt_35_3
                    target[target_ixyz+10] = radial * dx*dx*dz*dz * sqrt_35_3
                    target[target_ixyz+11] = radial * dy*dy*dz*dz * sqrt_35_3
                    target[target_ixyz+12] = radial * dx*dx*dy*dz * sqrt_35_1
                    target[target_ixyz+13] = radial * dy*dy*dx*dz * sqrt_35_1
                    target[target_ixyz+14] = radial * dz*dz*dx*dy * sqrt_35_1


cdef tuple checked_shells(shells, Py_ssize_t nao):
    owners = tuple(shells)
    if not owners:
        raise ValueError("Expected at least one basis shell")
    for item in owners:
        if not isinstance(item, Shell):
            raise TypeError("Expected a native basis shell")
        if (<Shell>item).radial_view is None or (<Shell>item).start > nao - (<Shell>item).size:
            raise ValueError("Shell AO range exceeds the coefficient array")
    return owners

cdef class OrbitalCalculater:
    """Computes value of orbital at real space points"""
    cdef tuple shells
    cdef object coeff
    cdef DTYPE_t[:] logmxcoeff
    cdef Py_ssize_t nshells
    cdef Py_ssize_t max_primitives
    cdef DTYPE_t[::1] _coeff_c

    def __init__(self, shells, coeff):
        self.coeff = np.array(coeff, dtype=DTYPE, order="C", copy=True)
        if self.coeff.ndim != 1 or not np.isfinite(self.coeff).all():
            raise ValueError("Expected a finite vector of orbital coefficients")
        self._coeff_c = self.coeff
        self.shells = checked_shells(shells, len(self.coeff))
        self.nshells = len(self.shells)
        self.max_primitives = max((<Shell>sh).exponent_view.shape[0] for sh in self.shells)

        # use very small number to replace zeros
        self.logmxcoeff = np.log(np.array([
            max(np.max(abs(self.coeff[sh.start:sh.start + sh.size])), 1e-30)
            for sh in self.shells], dtype=DTYPE))

    def bounding_box(self, thr=1.0e-5):
        """returns lower and upper limits of box (in bohr) that should fully contain orbital"""
        # find lower bounds
        p0 = [np.min([sh.center[ixyz] - sh.bounding_box_size(thr, lmx)
                      for sh, lmx in zip(self.shells, self.logmxcoeff)]) for ixyz in range(3)]

        # find upper bounds
        p1 = [np.max([sh.center[ixyz] + sh.bounding_box_size(thr, lmx)
                      for sh, lmx in zip(self.shells, self.logmxcoeff)]) for ixyz in range(3)]

        return p0, p1

    cpdef value(self, x, y, z):
        """Returns value of the orbital at specified point"""
        out = 0.0
        for sh, lmx in zip(self.shells, self.logmxcoeff):
            val = sh.value(x, y, z, self.coeff[sh.start:sh.start+sh.size], logmxcoeff=lmx)
            out += val
        return out

    cpdef plane_values(self, xvals, yvals, z_a):
        out = np.zeros([len(xvals), len(yvals)], dtype=DTYPE)
        xvals, yvals = np.ascontiguousarray(xvals, dtype=DTYPE), np.ascontiguousarray(yvals, dtype=DTYPE)
        for sh, lmx in zip(self.shells, self.logmxcoeff):
            sh.add_plane_values(xvals, yvals, z_a, self.coeff[sh.start:sh.start+sh.size],
                    out, logmxcoeff=lmx)
        return out

    @cython.boundscheck(CYDEBUG)
    @cython.wraparound(CYDEBUG)
    @cython.initializedcheck(CYDEBUG)
    @cython.nonecheck(CYDEBUG)
    cpdef box_values(self, ndarray[DTYPE_t, ndim=1] xvals,
            ndarray[DTYPE_t, ndim=1] yvals,
            ndarray[DTYPE_t, ndim=1] zvals):
        if xvals is None or yvals is None or zvals is None:
            raise ValueError("Expected three grid axes")
        cdef Py_ssize_t nshells = self.nshells
        cdef Py_ssize_t ish, nx, ny, nz

        nx = len(xvals)
        ny = len(yvals)
        nz = len(zvals)

        cdef ndarray[DTYPE_t, ndim=3, mode="c"] out = np.zeros([nx, ny, nz], dtype=DTYPE)
        cdef DTYPE_t[::1] coeff_view = self._coeff_c

        if not out.size:
            return out
        xvals = np.ascontiguousarray(xvals)
        yvals = np.ascontiguousarray(yvals)
        zvals = np.ascontiguousarray(zvals)
        cdef Interval bounds[3]
        bounds[0] = axis_bounds(&xvals[0], nx)
        bounds[1] = axis_bounds(&yvals[0], ny)
        bounds[2] = axis_bounds(&zvals[0], nz)
        cdef ndarray[DTYPE_t, ndim=2] scratch = np.empty((ny + nz + 1, self.max_primitives), dtype=DTYPE)
        cdef Shell sh
        for ish in range(nshells):
            sh = self.shells[ish]
            add_box_values(sh, nx, &xvals[0], ny, &yvals[0], nz, &zvals[0],
                           &coeff_view[sh.start], &out[0,0,0], self.logmxcoeff[ish], bounds, &scratch[0,0])

        return out

    def isovalue_containing_proportion(self, values=[0.90], resolution=0.2*ang2bohr, box=None):
        if box is None:
            p0, p1 = self.bounding_box(1e-4)
        else:
            p0, p1 = box

        npoints = [ int(math.ceil((b - a)/resolution)) for a, b in zip(p0, p1) ]

        xvals, yvals, zvals = [ np.linspace(a, b, num=n, endpoint=True, dtype=DTYPE) for a, b, n in zip(p0, p1, npoints) ]
        dV = (xvals[1] - xvals[0]) * (yvals[1] - yvals[0]) * (zvals[1] - zvals[0])
        boxvalues = self.box_values(xvals, yvals, zvals).reshape(-1)

        return isovalue_containing_proportion(values, boxvalues, dV)

cdef class DensityCalculater:
    """Computes value of densities at real space points"""

    cdef:
        tuple shells
        np.ndarray density
        public double nelec
        np.ndarray logmxorb
        Py_ssize_t max_primitives
        Py_ssize_t nshells

    def __init__(self, list shells, np.ndarray[DTYPE_t, ndim=2] density, double nelec):
        if density is None or density.shape[0] != density.shape[1]:
            raise ValueError("Expected a square AO density matrix")
        if not np.isfinite(density).all():
            raise ValueError("Expected finite AO density coefficients")
        self.shells = checked_shells(shells, density.shape[0])
        self.max_primitives = max((<Shell>sh).exponent_view.shape[0] for sh in self.shells)
        self.density = density.copy(order="C")
        self.nelec = nelec

        self.nshells = len(self.shells)
        weights = density_shell_weights(self.shells, self.density)
        self.logmxorb = np.log(np.maximum(weights, 1.e-30)).astype(DTYPE)

    def bounding_box(self, DTYPE_t thr=1.0e-5):
        """returns lower and upper limits of box (in bohr) that should fully contain density"""
        cdef:
            list p0, p1
            int ixyz

        # find lower bounds
        p0 = [np.min([sh.center[ixyz] - sh.bounding_box_size(thr, lmx)
               for sh, lmx in zip(self.shells, self.logmxorb)]) for ixyz in range(3)]

        # find upper bounds
        p1 = [np.max([sh.center[ixyz] + sh.bounding_box_size(thr, lmx)
               for sh, lmx in zip(self.shells, self.logmxorb)]) for ixyz in range(3)]

        return p0, p1

    cpdef DTYPE_t value(self, DTYPE_t x, DTYPE_t y, DTYPE_t z):
        """Sample the same screened AO contraction as a one-pixel density plane."""
        return self.plane_values(np.asarray([x], dtype=DTYPE), np.asarray([y], dtype=DTYPE), z)[0, 0]

    @cython.boundscheck(CYDEBUG)
    @cython.wraparound(CYDEBUG)
    @cython.initializedcheck(CYDEBUG)
    @cython.nonecheck(CYDEBUG)
    def plane_values(self, np.ndarray[DTYPE_t, ndim=1] xvals,
                    np.ndarray[DTYPE_t, ndim=1] yvals,
                    DTYPE_t z_a):
        if xvals is None or yvals is None:
            raise ValueError("Expected two plane axes")
        if not len(xvals) or not len(yvals):
            return np.zeros((len(xvals), len(yvals)), dtype=DTYPE)
        xvals, yvals = np.ascontiguousarray(xvals), np.ascontiguousarray(yvals)
        cdef:
            np.ndarray[DTYPE_t, ndim=2] out = np.zeros([len(xvals), len(yvals)], dtype=DTYPE)
            DTYPE_t[::1] oview = out.reshape(len(xvals) * len(yvals))
            Py_ssize_t i, j, ish, jsh, ixy
            Py_ssize_t nao = self.density.shape[0]
            np.ndarray[DTYPE_t, ndim=3, mode="c"] chivals = np.zeros([nao, len(xvals), len(yvals)], dtype=DTYPE)
            DTYPE_t[:,::1] chivals_f = chivals.reshape(nao, out.size)
            float[:,:] rho = self.density
            Shell ishell, jshell
            int nshls = len(self.shells)
            Py_ssize_t nxy = len(xvals) * len(yvals)
            Py_ssize_t istart, iend, jstart, jend

            bint[:] contributes = np.zeros([nshls], dtype=np.int32)
            bint cont

        for ish in range(nshls):
            ishell = self.shells[ish]
            istart = ishell.start
            iend = istart + ishell.size
            cont = ishell.compute_chi_plane(xvals, yvals, z_a, chivals[istart:iend],
                                  logmxcoeff=self.logmxorb[ish])
            contributes[ish] = cont

        for ish in range(nshls):
            if not contributes[ish]:
                continue
            ishell = self.shells[ish]
            istart = ishell.start
            iend = istart + ishell.size
            for jsh in range(nshls):
                if not contributes[jsh]:
                    continue

                jshell = self.shells[jsh]
                jstart = jshell.start
                jend = jstart + jshell.size
                for i in range(istart, iend):
                    for j in range(jstart, jend):
                        for ixy in range(nxy):
                            oview[ixy] = oview[ixy] + chivals_f[i,ixy] * chivals_f[j,ixy] * rho[i,j]
        return out

    @cython.boundscheck(CYDEBUG)
    @cython.wraparound(CYDEBUG)
    @cython.initializedcheck(CYDEBUG)
    @cython.nonecheck(CYDEBUG)
    cpdef box_values(self, np.ndarray[DTYPE_t, ndim=1] xvals,
                  np.ndarray[DTYPE_t, ndim=1] yvals,
                  np.ndarray[DTYPE_t, ndim=1] zvals):
        if xvals is None or yvals is None or zvals is None:
            raise ValueError("Expected three grid axes")
        cdef:
            np.ndarray[DTYPE_t, ndim=3] out = np.zeros([len(xvals), len(yvals), len(zvals)], dtype=DTYPE)
            Py_ssize_t nshells = self.nshells
            Py_ssize_t nx = len(xvals)
            Py_ssize_t ny = len(yvals)
            Py_ssize_t nz = len(zvals)
            Py_ssize_t nxyz = out.size
            Py_ssize_t ixyz, ish, i, j, istart, isize, offset
            float otmp

            DTYPE_t[::1] logmxorb = self.logmxorb
            DTYPE_t[::1] oview = out.reshape(nxyz)

            float logthresh, min_x, min_y, min_z, min_rr
            Interval bounds[3]
            Shell sh
            np.ndarray[np.intp_t, ndim=1] active_shells_arr, compact_offs
            Py_ssize_t n_active = 0
            Py_ssize_t nactive_ao = 0

        if not nxyz:
            return out
        xvals = np.ascontiguousarray(xvals)
        yvals = np.ascontiguousarray(yvals)
        zvals = np.ascontiguousarray(zvals)
        bounds[0] = axis_bounds(&xvals[0], nx)
        bounds[1] = axis_bounds(&yvals[0], ny)
        bounds[2] = axis_bounds(&zvals[0], nz)
        active_shells_arr = np.empty(nshells, dtype=np.intp)
        compact_offs = np.empty(nshells, dtype=np.intp)
        for ish in range(nshells):
            sh = self.shells[ish]
            logthresh = sh.logthr + sh.screen_log_weight * logmxorb[ish]
            min_x = distance_to_interval(sh.X, bounds[0])
            min_y = distance_to_interval(sh.Y, bounds[1])
            min_z = distance_to_interval(sh.Z, bounds[2])
            min_rr = min_x*min_x + min_y*min_y + min_z*min_z
            if min_rr < logthresh:
                active_shells_arr[n_active] = ish
                compact_offs[n_active] = nactive_ao
                n_active += 1
                nactive_ao += sh.size

        if nactive_ao == 0:
            return out

        # Large response/density grids would otherwise materialize two full
        # (grid points × active AOs) matrices. Keep complete z rows when possible.
        cdef Py_ssize_t capacity = max(1, _DENSITY_WORKSPACE_BYTES // (2 * sizeof(DTYPE_t) * nactive_ao))
        cdef Py_ssize_t bx, by, bz, ix, iy, iz
        if nxyz > capacity:
            bz = min(nz, capacity)
            by = min(ny, capacity // bz)
            bx = min(nx, capacity // (by * bz))
            for ix in range(0, nx, bx):
                for iy in range(0, ny, by):
                    for iz in range(0, nz, bz):
                        out[ix:ix+bx, iy:iy+by, iz:iz+bz] = self.box_values(
                            xvals[ix:ix+bx], yvals[iy:iy+by], zvals[iz:iz+bz])
            return out

        cdef np.ndarray[DTYPE_t, ndim=4] chivals = np.zeros([nx, ny, nz, nactive_ao], dtype=DTYPE)

        cdef ndarray[DTYPE_t, ndim=2] scratch = np.empty((ny + nz + 1, self.max_primitives), dtype=DTYPE)
        for i in range(n_active):
            ish = active_shells_arr[i]
            offset = compact_offs[i]
            sh = self.shells[ish]
            compute_chi_box_shell(sh, nx, &xvals[0], ny, &yvals[0], nz, &zvals[0],
                    &chivals[0,0,0,offset], nactive_ao*ny*nz, nactive_ao*nz, nactive_ao, logmxorb[ish], bounds, &scratch[0,0])

        active_ao_list = []
        for i in range(n_active):
            ish = active_shells_arr[i]
            sh = self.shells[ish]
            istart = sh.start
            isize = sh.size
            for j in range(istart, istart + isize):
                active_ao_list.append(j)

        idx = np.array(active_ao_list, dtype=np.intp)
        cdef np.ndarray[DTYPE_t, ndim=2] rho_compact = np.ascontiguousarray(
            self.density[np.ix_(idx, idx)])

        cdef DTYPE_t[:,::1] chivals_f = chivals.reshape(nxyz, nactive_ao)
        phi = np.dot(chivals_f, rho_compact)
        cdef DTYPE_t[:,::1] phiview = phi
        for ixyz in range(nxyz):
            otmp = 0.0
            for i in range(nactive_ao):
                otmp += chivals_f[ixyz, i] * phiview[ixyz, i]
            oview[ixyz] = otmp

        return out

    def isovalue_containing_proportion(self, list values=[0.90], DTYPE_t resolution=0.2*ang2bohr, tuple box=None):
        cdef:
            list p0, p1, npoints
            np.ndarray[DTYPE_t, ndim=1] xvals, yvals, zvals
            DTYPE_t dV
            np.ndarray boxvalues

        if box is None:
            p0, p1 = self.bounding_box(1e-7)
        else:
            p0, p1 = box

        npoints = [int(ceil((b - a)/resolution)) for a, b in zip(p0, p1)]

        xvals, yvals, zvals = [np.linspace(a, b, num=n, endpoint=True, dtype=DTYPE)
                              for a, b, n in zip(p0, p1, npoints)]
        dV = (xvals[1] - xvals[0]) * (yvals[1] - yvals[0]) * (zvals[1] - zvals[0])
        boxvalues = self.box_values(xvals, yvals, zvals).reshape(-1)

        return isovalue_containing_proportion(values, boxvalues, dV, square=False)
