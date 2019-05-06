# cython: profile=True
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

import numpy as np
cimport numpy as np

from numpy cimport ndarray
from libc.math cimport sqrt, exp

DTYPE = np.float32
ctypedef np.float32_t DTYPE_t

import math

from .constants import ang2bohr

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
    cdef float X, Y, Z
    cdef public int l, start, size
    cdef public DTYPE_t[:] center
    cdef DTYPE_t[:] exponents
    cdef DTYPE_t[:] prim_coeffs
    cdef DTYPE_t[:,:] norms
    cdef DTYPE_t[:] denormed_coeffs
    cdef float mxx
    cdef float most_diffuse
    cdef float logthr
    cdef float mxnorm

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

    def value(self, DTYPE_t x, DTYPE_t y, DTYPE_t z, ndarray[DTYPE_t, ndim=1] coeff,
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

    @cython.boundscheck(False)
    @cython.wraparound(False)
    cpdef add_plane_values(self, ndarray[DTYPE_t, ndim=1] xx,
            ndarray[DTYPE_t, ndim=1] yy, float z_a,
            ndarray[DTYPE_t, ndim=1] coeff,
            ndarray[DTYPE_t, ndim=2] target, float logmxcoeff=0.0):

        cdef float X = self.X
        cdef float Y = self.Y
        cdef float Z = self.Z

        cdef float sqrt3 = sqrt(3.0)
        cdef float sqrt5 = sqrt(5.0)
        cdef float sqrt15 = sqrt(15.0)
        cdef float sqrt_35_5 = sqrt(35.0/5.0)
        cdef float sqrt_35_3 = sqrt(35.0/3.0)
        cdef float sqrt_35_1 = sqrt(35.0/1.0)

        cdef float logthresh = self.logthr

        cdef int nx = len(xx)
        cdef int ny = len(yy)

        cdef int nexp = len(self.exponents)

        cdef int l = self.l
        cdef int shellsize = self.size

        # check the four corners to find a max value
        cdef float x0 = xx[0]
        cdef float x1 = xx[nx-1]
        cdef float y0 = yy[0]
        cdef float y1 = yy[ny-1]

        cdef float min_x = min(abs(x0),abs(x1))
        if (x0*x1 < 0):
            min_x = 0.0
        cdef float min_y = min(abs(y0),abs(y1))
        if (y0*y1 < 0):
            min_y = 0.0
        cdef float min_z = z_a

        cdef float min_rr = min_x*min_x + min_y*min_y + min_z*min_z
        if min_rr >= logthresh - logmxcoeff:
            return

        cdef float dx, dy, dz, rr, radial

        cdef DTYPE_t[:] denormed_coeffs = self.denormed_coeffs
        cdef DTYPE_t[:] exponents = self.exponents

        cdef ndarray[DTYPE_t, ndim=2] expyy = np.zeros([ny, nexp], dtype=DTYPE)
        cdef ndarray[DTYPE_t, ndim=1] expxx = np.zeros(nexp, dtype=DTYPE)

        for iy in range(ny):
            dy = yy[iy] - Y
            for iexp in range(nexp):
                expyy[iy,iexp] = exp(-exponents[iexp]*dy*dy)

        dz = z_a - Z
        for ix in range(nx):
            dx = xx[ix] - X
            for iexp in range(nexp):
                expxx[iexp] = exp(-exponents[iexp] * (dx*dx + dz*dz)) * denormed_coeffs[iexp]

            for iy in range(ny):
                dy = yy[iy] - Y

                radial = 0.0
                for iexp in range(nexp):
                    radial += expxx[iexp] * expyy[iy,iexp]

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


class MOData(object):
    """Organizes data for computing orbitals in real space"""

    def __init__(self, shells, coeff, nocc):
        assert sum([sh.size for sh in shells]) == coeff.shape[1]

        self.shells = shells
        self.coeff = coeff
        self.nocc = nocc

    def homo(self):
        return self.nocc

    def get_orbital(self, iorb):
        """Returns OrbitalCalculater for orbital"""
        i = (iorb if iorb > 0 else self.nocc + iorb) - 1
        return OrbitalCalculater(self.shells, self.coeff[i, :])

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

        nocc = sum([float(x["occup"]) > 0.0 for x in geo["mo"]])
        nmo = len(geo["mo"])
        nao = sum([sh.size for sh in shells])

        coeff = np.zeros([nmo, nao], dtype=DTYPE)
        for i in range(nmo):
            coeff[i, :] = geo["mo"][i]["coeff"]

        return cls(shells, coeff, nocc)


class OrbitalCalculater(object):
    """Computes value of orbital at real space points"""

    def __init__(self, shells, coeff):
        self.shells = shells
        self.coeff = np.array(coeff, dtype=DTYPE)

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

    def value(self, x, y, z):
        """Returns value of the orbital at specified point"""
        out = 0.0
        for sh, lmx in zip(self.shells, self.logmxcoeff):
            val = sh.value(x, y, z, self.coeff[sh.start:sh.start+sh.size], logmxcoeff=lmx)
            out += val
        return out

    def plane_values(self, xvals, yvals, z_a):
        out = np.zeros([len(xvals), len(yvals)], dtype=DTYPE)
        for sh, lmx in zip(self.shells, self.logmxcoeff):
            sh.add_plane_values(xvals, yvals, z_a, self.coeff[sh.start:sh.start+sh.size],
                    out, logmxcoeff=lmx)
        return out

    def box_values(self, xvals, yvals, zvals):
        out = np.zeros([len(xvals), len(yvals), len(zvals)], dtype=DTYPE)

        for k, z_a in enumerate(zvals):
            out[:,:,k] = self.plane_values(xvals, yvals, z_a)

        return out

    def isovalue_containing_proportion(self, values=[0.90], resolution=0.2*ang2bohr, box=None):
        if box is None:
            p0, p1 = self.bounding_box(1e-4)

        npoints = [ int(math.ceil((b - a)/resolution)) for a, b in zip(p0, p1) ]

        xvals, yvals, zvals = [ np.linspace(a, b, num=n, endpoint=True, dtype=DTYPE) for a, b, n in zip(p0, p1, npoints) ]
        dV = (xvals[1] - xvals[0]) * (yvals[1] - yvals[0]) * (zvals[1] - zvals[0])
        boxvalues = self.box_values(xvals, yvals, zvals).reshape(-1)

        return [0.1]
        return isovalue_containing_proportion(values, boxvalues, dV)
