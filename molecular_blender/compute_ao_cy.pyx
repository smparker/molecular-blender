#distutils: language = c++
#cython: cdivision=True
#cython: profile=False
# -*- coding: utf-8 -*-
#
#  Molecular Blender
#  Filename: orbitals.py
#  Copyright (C) 2017 Shane Parker, Joshua Szekely
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
from libc.stdlib cimport malloc, free

DTYPE = np.float32
ctypedef np.float32_t DTYPE_t

DEF CYDEBUG = False

from .constants import ang2bohr
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

cdef class Shell:
    """ Fast Cython container for shell data"""
    cdef public float X, Y, Z
    cdef public int l, start, size
    cdef public DTYPE_t[::1] center
    cdef public object exponents
    cdef public object prim_coeffs
    cdef public object norms
    cdef public object denormed_coeffs
    cdef public float mxx
    cdef public float most_diffuse
    cdef public float logthr
    cdef public float mxnorm

    def __init__(self, center, int l,
            exponents, coeff,
            int start, float thr=1e-6, float mxx=1.0):
        self.center = np.array(center, dtype=DTYPE)

        self.X = center[0]
        self.Y = center[1]
        self.Z = center[2]

        self.l = int(l)

        self.start = start
        self.size = (self.l + 1) * (self.l + 2) // 2

        self.exponents = np.array(exponents, dtype=DTYPE)
        self.prim_coeffs = np.array(coeff, dtype=DTYPE)

        self.norms = 1.0 / np.sqrt(self.__norms(), dtype=DTYPE)
        self.denormed_coeffs = np.array(coeff, dtype=DTYPE) / np.sqrt(self.__shellnorms(), dtype=DTYPE)

        # most diffuse exponent in shell
        cdef int imin = np.argmin(self.exponents)
        self.mxx = mxx  # save max x just for use later
        self.most_diffuse = np.min(self.exponents)
        self.logthr = np.log(thr * mxx**self.l, dtype=DTYPE)
        self.mxnorm = np.max(self.norms[:, imin]) * abs(self.prim_coeffs[imin])

        self.logthr -= np.log(self.mxnorm)
        self.logthr /= -self.most_diffuse

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
        xx = self.l * (np.log(self.mxx) - np.log(thr / self.mxnorm) + logmxcoeff) / self.most_diffuse
        if xx < 0:
            return 0.0
        return np.sqrt(xx, dtype=DTYPE)

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

        cdef float logthresh = self.logthr

        cdef int nexp = len(self.exponents)
        cdef int l = self.l
        cdef int shellsize = self.size

        cdef DTYPE_t[:] exponents = self.exponents
        cdef DTYPE_t[:] denormed_coeffs = self.denormed_coeffs

        cdef float rr = dx*dx + dy*dy + dz*dz

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
    cpdef void add_plane_values(self, DTYPE_t[::1] xx,
            DTYPE_t[::1] yy, float z_a,
            DTYPE_t[::1] coeff,
            DTYPE_t[:,::1] target, float logmxcoeff=0.0):

        cdef float X = self.X
        cdef float Y = self.Y
        cdef float Z = self.Z

        cdef float logthresh = self.logthr - logmxcoeff

        cdef int nx = len(xx)
        cdef int ny = len(yy)

        cdef int nexp = len(self.exponents)

        cdef int l = self.l
        cdef int shellsize = self.size

        # check the four corners to find a max value
        cdef float x0 = xx[0] - X
        cdef float x1 = xx[nx-1] - X
        cdef float y0 = yy[0] - Y
        cdef float y1 = yy[ny-1] - Y

        cdef float min_x = min(abs(x0),abs(x1))
        if (x0*x1 < 0):
            min_x = 0.0
        cdef float min_y = min(abs(y0),abs(y1))
        if (y0*y1 < 0):
            min_y = 0.0
        cdef float min_z = z_a

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

        cdef DTYPE_t[::1] denormed_coeffs = self.denormed_coeffs
        cdef DTYPE_t[::1] exponents = self.exponents

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
    cpdef bint compute_chi_plane(self, DTYPE_t[::1] xx,
            DTYPE_t[::1] yy, float z_a,
            DTYPE_t[:,:,::1] target, float logmxcoeff=0.0):

        cdef float x = self.X
        cdef float y = self.Y
        cdef float z = self.Z

        cdef float logthresh = self.logthr - logmxcoeff

        cdef int nx = len(xx)
        cdef int ny = len(yy)

        cdef int nexp = len(self.exponents)

        cdef int l = self.l
        cdef int shellsize = self.size

        # check the four corners to find a max value
        cdef float x0 = xx[0] - x
        cdef float x1 = xx[nx-1] - x
        cdef float y0 = yy[0] - y
        cdef float y1 = yy[ny-1] - y

        cdef float min_x = min(abs(x0),abs(x1))
        if (x0*x1 < 0):
            min_x = 0.0
        cdef float min_y = min(abs(y0),abs(y1))
        if (y0*y1 < 0):
            min_y = 0.0
        cdef float min_z = z_a

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

        cdef DTYPE_t[::1] denormed_coeffs = self.denormed_coeffs
        cdef DTYPE_t[::1] exponents = self.exponents

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


ctypedef struct CShell:
    float X, Y, Z
    int l
    int start, size
    int nexp
    DTYPE_t* exponents
    DTYPE_t* prim_coeffs
    DTYPE_t* norms
    DTYPE_t* denormed_coeffs
    float mxx
    float most_diffuse
    float logthr
    float mxnorm

@cython.boundscheck(CYDEBUG)
@cython.wraparound(CYDEBUG)
@cython.initializedcheck(CYDEBUG)
@cython.nonecheck(CYDEBUG)
cdef void add_box_values(CShell* shell,
        int nx, DTYPE_t* xx,
        int ny, DTYPE_t* yy,
        int nz, DTYPE_t* zz,
        DTYPE_t* coeff,
        DTYPE_t* target,
        float logmxcoeff) noexcept nogil:
    cdef float X = shell.X
    cdef float Y = shell.Y
    cdef float Z = shell.Z

    cdef float logthresh = shell.logthr - logmxcoeff

    cdef int ix
    cdef int iy
    cdef int iz
    cdef int iexp

    cdef int nexp = shell.nexp

    cdef int l = shell.l

    # check the four corners to find a max value
    cdef float x0 = xx[0] - X
    cdef float x1 = xx[nx-1] - X
    cdef float y0 = yy[0] - Y
    cdef float y1 = yy[ny-1] - Y
    cdef float z0 = zz[0] - Z
    cdef float z1 = zz[nz-1] - Z

    cdef float min_x = min(fabs(x0),fabs(x1))
    if (x0*x1 < 0):
        min_x = 0.0
    cdef float min_y = min(fabs(y0),fabs(y1))
    if (y0*y1 < 0):
        min_y = 0.0
    cdef float min_z = min(fabs(z0),fabs(z1))
    if (z0*z1 < 0):
        min_z = 0.0

    cdef float min_rr = min_x*min_x + min_y*min_y + min_z*min_z
    if min_rr >= logthresh:
        return

    cdef float sqrt3 = sqrt(3.0)
    cdef float sqrt5 = sqrt(5.0)
    cdef float sqrt15 = sqrt(15.0)
    cdef float sqrt_35_5 = sqrt(35.0/5.0)
    cdef float sqrt_35_3 = sqrt(35.0/3.0)
    cdef float sqrt_35_1 = sqrt(35.0/1.0)

    cdef float dx, dy, dz, rr, radial

    cdef DTYPE_t* denormed_coeffs = shell.denormed_coeffs
    cdef DTYPE_t* exponents = shell.exponents

    cdef DTYPE_t* expzz
    cdef DTYPE_t* expyy
    cdef DTYPE_t* expxx

    expzz = <DTYPE_t*>malloc(nz*nexp*sizeof(DTYPE_t))
    expyy = <DTYPE_t*>malloc(ny*nexp*sizeof(DTYPE_t))
    expxx = <DTYPE_t*>malloc(nexp*sizeof(DTYPE_t))

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
                if zthresh > logthresh:
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
    free(expxx)
    free(expyy)
    free(expzz)

@cython.boundscheck(CYDEBUG)
@cython.wraparound(CYDEBUG)
@cython.initializedcheck(CYDEBUG)
@cython.nonecheck(CYDEBUG)
cdef bint compute_chi_box_shell(CShell* shell,
        int nx, DTYPE_t* xx,
        int ny, DTYPE_t* yy,
        int nz, DTYPE_t* zz,
        DTYPE_t* target,
        int sx, int sy, int sz,
        float logmxcoeff) noexcept nogil:
    """Computes the value of the shell on a grid

    target is a 3D array of shape (nx, ny, nz, shellsize)
    sx, sy, sz are strides, such that target[ix,iy,iz,ishell] is at
    target[ix*sx + iy*sy + iz*sz + ishell]
    """
    cdef bint out = False

    cdef float X = shell.X
    cdef float Y = shell.Y
    cdef float Z = shell.Z

    cdef float logthresh = shell.logthr - logmxcoeff

    cdef int ix
    cdef int iy
    cdef int iz
    cdef int iexp
    cdef int target_ixyz

    cdef int nexp = shell.nexp

    cdef int l = shell.l

    # check the four corners to find a max value
    cdef float x0 = xx[0] - X
    cdef float x1 = xx[nx-1] - X
    cdef float y0 = yy[0] - Y
    cdef float y1 = yy[ny-1] - Y
    cdef float z0 = zz[0] - Z
    cdef float z1 = zz[nz-1] - Z

    cdef float min_x = min(fabs(x0),fabs(x1))
    if (x0*x1 < 0):
        min_x = 0.0
    cdef float min_y = min(fabs(y0),fabs(y1))
    if (y0*y1 < 0):
        min_y = 0.0
    cdef float min_z = min(fabs(z0),fabs(z1))
    if (z0*z1 < 0):
        min_z = 0.0

    cdef float min_rr = min_x*min_x + min_y*min_y + min_z*min_z
    if min_rr >= logthresh:
        return out

    cdef float sqrt3 = sqrt(3.0)
    cdef float sqrt5 = sqrt(5.0)
    cdef float sqrt15 = sqrt(15.0)
    cdef float sqrt_35_5 = sqrt(35.0/5.0)
    cdef float sqrt_35_3 = sqrt(35.0/3.0)
    cdef float sqrt_35_1 = sqrt(35.0/1.0)

    cdef float dx, dy, dz, rr, radial

    cdef DTYPE_t* denormed_coeffs = shell.denormed_coeffs
    cdef DTYPE_t* exponents = shell.exponents

    cdef DTYPE_t* expzz
    cdef DTYPE_t* expyy
    cdef DTYPE_t* expxx

    expzz = <DTYPE_t*>malloc(nz*nexp*sizeof(DTYPE_t))
    expyy = <DTYPE_t*>malloc(ny*nexp*sizeof(DTYPE_t))
    expxx = <DTYPE_t*>malloc(nexp*sizeof(DTYPE_t))

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
                if zthresh > logthresh:
                    continue

                out = True

                ixyz = iz + iy*nz + ix*nz*ny

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

    free(expxx)
    free(expyy)
    free(expzz)

    return out

cdef class OrbitalCalculater:
    """Computes value of orbital at real space points"""
    cdef object shells
    cdef object coeff
    cdef DTYPE_t[:] logmxcoeff
    cdef int nshells
    cdef CShell* cshells
    cdef Shell sh

    def __init__(self, shells, coeff):
        self.shells = shells
        self.coeff = np.array(coeff, dtype=DTYPE)

        cdef ndarray[DTYPE_t, ndim=1] sh_exponents
        cdef ndarray[DTYPE_t, ndim=1] sh_prim_coeffs
        cdef ndarray[DTYPE_t, ndim=2] sh_norms
        cdef ndarray[DTYPE_t, ndim=1] sh_denormed_coeffs

        self.nshells = len(shells)
        self.cshells = <CShell*>malloc(self.nshells * sizeof(CShell))
        for ish in range(self.nshells):
            sh = shells[ish]
            sh_exponents = sh.exponents
            sh_prim_coeffs = sh.prim_coeffs
            sh_norms = sh.norms
            sh_denormed_coeffs = sh.denormed_coeffs

            self.cshells[ish].X = sh.X
            self.cshells[ish].Y = sh.Y
            self.cshells[ish].Z = sh.Z
            self.cshells[ish].l = sh.l
            self.cshells[ish].start = sh.start
            self.cshells[ish].size = sh.size
            self.cshells[ish].nexp = len(sh.exponents)
            self.cshells[ish].exponents = &sh_exponents[0]
            self.cshells[ish].prim_coeffs = &sh_prim_coeffs[0]
            self.cshells[ish].norms = &sh_norms[0,0]
            self.cshells[ish].denormed_coeffs = &sh_denormed_coeffs[0]
            self.cshells[ish].mxx = sh.mxx
            self.cshells[ish].most_diffuse = sh.most_diffuse
            self.cshells[ish].logthr = sh.logthr
            self.cshells[ish].mxnorm = sh.mxnorm

        # use very small number to replace zeros
        self.logmxcoeff = np.log(np.array([
            max(np.max(abs(coeff[sh.start:sh.start + sh.size])), 1e-30)
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
        for sh, lmx in zip(self.shells, self.logmxcoeff):
            sh.add_plane_values(xvals, yvals, z_a, self.coeff[sh.start:sh.start+sh.size],
                    out, logmxcoeff=lmx)
        return out

    @cython.boundscheck(CYDEBUG)
    @cython.wraparound(CYDEBUG)
    @cython.initializedcheck(CYDEBUG)
    @cython.nonecheck(CYDEBUG)
    cpdef box_values(self, ndarray[DTYPE_t, ndim=1, mode="c"] xvals,
            ndarray[DTYPE_t, ndim=1, mode="c"] yvals,
            ndarray[DTYPE_t, ndim=1, mode="c"] zvals):
        cdef int nshells = self.nshells
        cdef int ish, nx, ny, nz

        nx = len(xvals)
        ny = len(yvals)
        nz = len(zvals)

        cdef ndarray[DTYPE_t, ndim=3, mode="c"] out = np.zeros([nx, ny, nz], dtype=DTYPE)
        cdef ndarray[DTYPE_t, ndim=1, mode="c"] coeff = np.array(self.coeff, dtype=DTYPE)

        for ish in range(nshells):
            add_box_values(&self.cshells[ish], nx, &xvals[0], ny, &yvals[0], nz, &zvals[0],
                    &coeff[self.cshells[ish].start], &out[0,0,0], self.logmxcoeff[ish])

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
        list shells
        np.ndarray density
        int nelec
        np.ndarray shellmxdens
        np.ndarray logmxdens
        np.ndarray shellmxorb
        np.ndarray logmxorb
        CShell* cshells
        int nshells

    def __init__(self, list shells, np.ndarray[DTYPE_t, ndim=2] density, int nelec):
        self.shells = shells
        self.density = density
        self.nelec = nelec

        # use very small number to replace zeros
        self.nshells = len(self.shells)
        # shellmxdens[i, j] is the maximum value in density of shell i, j
        self.shellmxdens = np.zeros([self.nshells, self.nshells], dtype=DTYPE)

        cdef:
            int i, j
            int istart, iend, jstart, jend
            DTYPE_t shelldens

        for i in range(self.nshells):
            istart = shells[i].start
            iend = istart + shells[i].size
            for j in range(self.nshells):
                jstart = shells[j].start
                jend = jstart + shells[j].size
                shelldens = np.max(np.abs(self.density[istart:iend, jstart:jend]))
                self.shellmxdens[i,j] = shelldens

        self.logmxdens = np.log(self.shellmxdens + 1e-30)

        # shellmxorb[i] of shell i in entire density
        self.shellmxorb = np.array([np.max(self.shellmxdens[i,:]) for i in range(self.nshells)], dtype=DTYPE)
        # log(shellmxorb)
        self.logmxorb = np.log(self.shellmxorb + 1e-30)

        cdef ndarray[DTYPE_t, ndim=1] sh_exponents
        cdef ndarray[DTYPE_t, ndim=1] sh_prim_coeffs
        cdef ndarray[DTYPE_t, ndim=2] sh_norms
        cdef ndarray[DTYPE_t, ndim=1] sh_denormed_coeffs

        self.nshells = len(shells)
        self.cshells = <CShell*>malloc(self.nshells * sizeof(CShell))
        for ish in range(self.nshells):
            sh = shells[ish]
            sh_exponents = sh.exponents
            sh_prim_coeffs = sh.prim_coeffs
            sh_norms = sh.norms
            sh_denormed_coeffs = sh.denormed_coeffs

            self.cshells[ish].X = sh.X
            self.cshells[ish].Y = sh.Y
            self.cshells[ish].Z = sh.Z
            self.cshells[ish].l = sh.l
            self.cshells[ish].start = sh.start
            self.cshells[ish].size = sh.size
            self.cshells[ish].nexp = len(sh.exponents)
            self.cshells[ish].exponents = &sh_exponents[0]
            self.cshells[ish].prim_coeffs = &sh_prim_coeffs[0]
            self.cshells[ish].norms = &sh_norms[0,0]
            self.cshells[ish].denormed_coeffs = &sh_denormed_coeffs[0]
            self.cshells[ish].mxx = sh.mxx
            self.cshells[ish].most_diffuse = sh.most_diffuse
            self.cshells[ish].logthr = sh.logthr
            self.cshells[ish].mxnorm = sh.mxnorm

    def bounding_box(self, DTYPE_t thr=1.0e-5):
        """returns lower and upper limits of box (in bohr) that should fully contain orbital"""
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
        """Returns value of the density at specified point"""
        cdef:
            DTYPE_t out = 0.0
            list nonzero_shells = []
            int nao = sum([sh.size for sh in self.shells])
            np.ndarray[DTYPE_t, ndim=1] phivals = np.zeros(nao, dtype=DTYPE)
            int i, j
            np.ndarray[DTYPE_t, ndim=1] ivals

        for i, ish in enumerate(self.shells):
            ivals = ish.value(x, y, z, logmxcoeff=self.logmxorb[i])
            if ivals is not None:
                phivals[ish.start:ish.start + ish.size] = ivals
                nonzero_shells.append([i, ish])

        for i, ish in nonzero_shells:
            for j, jsh in nonzero_shells:
                out += np.einsum("p,q,pq->",
                               phivals[ish.start:ish.start + ish.size],
                               phivals[jsh.start:jsh.start + jsh.size],
                               self.density[ish.start:ish.start + ish.size,
                                          jsh.start:jsh.start + jsh.size],
                                 optimize=True)
        return out

    @cython.boundscheck(CYDEBUG)
    @cython.wraparound(CYDEBUG)
    @cython.initializedcheck(CYDEBUG)
    @cython.nonecheck(CYDEBUG)
    def plane_values(self, np.ndarray[DTYPE_t, ndim=1] xvals,
                    np.ndarray[DTYPE_t, ndim=1] yvals,
                    DTYPE_t z_a):
        cdef:
            np.ndarray[DTYPE_t, ndim=2] out = np.zeros([len(xvals), len(yvals)], dtype=DTYPE)
            DTYPE_t[::1] oview = out.reshape(len(xvals) * len(yvals))
            int i, j, ish, jsh, ixy
            int nao = sum([sh.size for sh in self.shells])
            float[:,:,:] chivals = np.zeros([nao, len(xvals), len(yvals)], dtype=DTYPE)
            float[:,:] chivals_f = (<float[:nao,:len(xvals)*len(yvals)]>(&chivals[0,0,0]))
            float[:,:] rho = self.density
            list nonzero_shells = []
            np.ndarray ivals
            int nshls = len(self.shells)
            int nxy = len(xvals) * len(yvals)
            int istart, iend, jstart, jend

            float[:] xxvals = xvals
            float[:] yyvals = yvals
            bint[:] contributes = np.zeros([nshls], dtype=np.int32)
            bint cont

        for ish in range(nshls):
            ishell = self.shells[ish]
            istart = self.cshells[ish].start
            iend = istart + self.cshells[ish].size
            cont = ishell.compute_chi_plane(xxvals, yyvals, z_a, chivals[istart:iend],
                                  logmxcoeff=self.logmxorb[ish])
            contributes[ish] = cont

        for ish in range(nshls):
            if not contributes[ish]:
                continue
            istart = self.cshells[ish].start
            iend = istart + self.cshells[ish].size
            for jsh in range(nshls):
                if not contributes[jsh]:
                    continue

                jstart = self.cshells[jsh].start
                jend = jstart + self.cshells[jsh].size
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
        cdef:
            np.ndarray[DTYPE_t, ndim=3] out = np.zeros([len(xvals), len(yvals), len(zvals)], dtype=DTYPE)
            int nshells = self.nshells
            DTYPE_t z_a
            int nx = len(xvals)
            int ny = len(yvals)
            int nz = len(zvals)
            int nxyz = nx * ny * nz
            int ixyz, ish, jsh, i, j, istart, iend, jstart, jend, isize, jsize
            int nao = sum([sh.size for sh in self.shells])

            DTYPE_t[:,::1] rho = self.density
            # notice that the shape of chivals is different from the one in plane_values
            DTYPE_t[:,:,:,::1] chivals = np.zeros([nx, ny, nz, nao], dtype=DTYPE)
            DTYPE_t[:,::1] chivals_f = (<float[:nx*ny*nz,:nao]>(&chivals[0,0,0,0]))

            DTYPE_t[::1] logmxorb = self.logmxorb
            #bint[:] contributes = np.zeros([nshells], dtype=np.int32)
            DTYPE_t[::1] oview = out.reshape(nx * ny * nz)
            #float rhoij

            float otmp;

        # TODO consider batching into smaller cubes
        for ish in range(nshells):
            istart = self.cshells[ish].start
            compute_chi_box_shell(&self.cshells[ish], nx, &xvals[0], ny, &yvals[0], nz, &zvals[0],
                    &chivals[0,0,0,istart], nao*ny*nz, nao*nz, nao, logmxorb[ish])

        phi = np.dot(chivals_f, rho)
        cdef DTYPE_t[:,:] phiview = phi
        for ixyz in range(nxyz):
            otmp = 0.0
            for i in range(nao):
                otmp += phiview[ixyz, i] * chivals_f[ixyz, i]
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
