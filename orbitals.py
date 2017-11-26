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

from .constants import ang2bohr, bohr2ang
from .util import stopwatch
from .periodictable import symbols

import numpy as np

def gamma2n(n2):
    """Returns Gamma(n/2), only for half-integers"""
    return np.array([ gamma2n_impl(int(round(x))) for x in n2 ])

def gamma2n_impl(n2):
    """Returns gamma2n_impl(n) where n is input as 2*n=n2"""
    if n2 <= 0:
        raise Exception("gamma2n_impl function cannot accept negative numbers")
    elif n2 == 2:
        return 1
    elif n2 == 1:
        return np.sqrt(np.pi)
    else:
        return (0.5*n2-1.) * gamma2n_impl(n2-2)

class Shell(object):
    """Container for shell data"""
    def __init__(self, center, l, exponents, coeff, start = 0, thr = 1e-6, mxx = 1.0):
        self.center = np.array(center)
        self.X, self.Y, self.Z = center[0], center[1], center[2]
        self.l = int(l)
        self.exponents = np.array(exponents)
        self.coeff = np.array(coeff)

        self.start = start
        self.size = (self.l+1)*(self.l+2)//2

        self.norms = 1.0/np.sqrt(self.__norms())
        self.__poly = [self.__poly_s, self.__poly_p, self.__poly_d, self.__poly_f, self.__poly_g][l]

        # most diffuse exponent in shell
        imin = np.argmin(self.exponents)
        self.mxx = mxx # save max x just for use later
        self.diffuse = np.min(self.exponents)
        self.logthr = np.log(thr * mxx**self.l)
        self.mxnorm = np.max(self.norms[:,imin]) * abs(self.coeff[imin])

        self.logthr -= np.log(self.mxnorm)
        self.logthr /= -self.diffuse

    shorder = [ [""],
            ["x", "y", "z"],
            ["xx", "yy", "zz", "xy", "xz", "yz"],
            ["xxx", "yyy", "zzz", "xyy", "xxy", "xxz", "xzz", "yzz", "yyz", "xyz"],
            ["xxxx", "yyyy", "zzzz", "xxxy", "xxxz", "yyyx", "yyyz", "zzzx", "zzzy",
                "xxyy", "xxzz", "yyzz", "xxyz", "yyxz", "zzxy"]
            ]

    def __norms(self):
        """Returns [nL, nexp] array of norms"""
        out = np.zeros([self.size, len(self.exponents)])
        for ixyz, xyz in enumerate(self.shorder[self.l]):
            cxyz = [xyz.count("x"), xyz.count("y"), xyz.count("z")]
            for ig, ex in enumerate(self.exponents):
                out[ixyz,ig] = gaussian_overlap(self.center, ex, cxyz, self.center, ex, cxyz)

        return out

    def __poly_s(self, X, Y, Z):
        """n/a"""
        return 1.0

    def __poly_p(self, X, Y, Z):
        """x, y, z"""
        return np.array([X, Y, Z])

    def __poly_d(self, X, Y, Z):
        """xx, yy, zz, xy, xz, yz"""
        return np.array([X*X, Y*Y, Z*Z, X*Y, X*Z, Y*Z])

    def __poly_f(self, X, Y, Z):
        """xxx, yyy, zzz, xyy, xxy, xxz, xzz, yzz, yyz, xyz"""
        return np.array([X*X*X, Y*Y*Y, Z*Z*Z, X*Y*Y, X*X*Y, X*Z*Z, Y*Z*Z, Y*Y*Z, X*Y*Z])

    def __poly_g(self, X, Y, Z):
        """xxxx yyyy zzzz xxxy xxxz yyyx yyyz zzzx zzzy xxyy xxzz yyzz xxyz yyxz zzxy"""
        return np.array([ X*X*X*X, Y*Y*Y*Y, Z*Z*Z*Z, X*X*X*Y, X*X*X*Z, Y*Y*Y*X, Y*Y*Y*Z,
            Z*Z*Z*X, Z*Z*Z*Y, X*X*Y*Y, X*X*Z*Z, Y*Y*Z*Z, X*X*Y*Z, Y*Y*X*Z, Z*Z*X*Y])

    def bounding_box_size(self, thr, logmxcoeff = 0.0):
        """returns half the edgelength of a box outside of which the shell is guaranteed to be below the given threshold"""
        xx = self.l * np.log(self.mxx) / self.diffuse - np.log(thr/self.mxnorm) / self.diffuse + logmxcoeff / self.diffuse
        if xx < 0:
            return 0.0
        else:
            return np.sqrt(xx)

    def value(self, x, y, z, logmxcoeff = 0.0):
        """compute values of shell at point, screening out contributions below shell threshold"""
        px, py, pz = x - self.X, y - self.Y, z - self.Z
        rr = px*px + py*py + pz*pz

        if rr < self.logthr - logmxcoeff:
            expz = np.exp(-self.exponents * rr) * self.coeff
            base = np.dot(self.norms, expz)
            return base * self.__poly(px, py, pz)
        else: # screen out tight functions
            return None

def gaussian_overlap(ax, aexp, ak, bx, bexp, bk):
    """returns overlap integral for gaussians centered at x with exponents exp
    and angular momentum k"""
    assert((ax == bx).all()) # too lazy to do two-center overlap for now
    kn = [ k + n for k, n in zip(ak,bk) ]
    if any([i%2 == 1 for i in kn]):
        return 0.0
    else:
        kn = np.array(kn) + 1.
        fac = np.power(aexp+bexp, -0.5*kn)
        gam = gamma2n(kn)
        return np.prod(fac*gam)

class MOData(object):
    """Organizes data for computing orbitals in real space"""
    def __init__(self, shells, coeff, nocc):
        assert(sum([sh.size for sh in shells]) == coeff.shape[1])

        self.shells = shells
        self.coeff = coeff
        self.nocc = nocc

    def get_orbital(self, iorb):
        """Returns OrbitalCalculater for orbital"""
        i = (iorb if iorb > 0 else self.nocc+iorb) - 1
        return OrbitalCalculater(self.shells, self.coeff[i,:])

    @classmethod
    def from_dict(cls, geo):
        """Builds MOData out of the dict object available in MolecularBlender"""
        shells = []
        ico = 0
        for a, ash in zip(geo["atoms"], geo["basis"]):
            xyz = [ x * ang2bohr for x in a["position"] ] # make sure centers are given in bohr
            for sh in ash:
                newshell = Shell(xyz, "spdfg".index(sh["shell"]), sh["exponents"], sh["contractions"], start=ico)
                ico += newshell.size
                shells.append(newshell)

        nocc = sum([ float(x["occup"]) > 0.0 for x in geo["mo"] ])
        nmo = len(geo["mo"])
        nao = sum([ sh.size for sh in shells ])

        coeff = np.zeros([nmo, nao])
        for i in range(nmo):
            coeff[i,:] = geo["mo"][i]["coeff"]

        return cls(shells, coeff, nocc)

class OrbitalCalculater(object):
    """Computes value of orbital at real space points"""

    def __init__(self, shells, coeff):
        self.shells = shells
        self.coeff = coeff

        # use very small number to replace zeros
        self.logmxcoeff = np.log(np.array([
            max(np.max(abs(self.coeff[sh.start:sh.start+sh.size])), 1e-30) for sh in self.shells ]))

    def bounding_box(self, thr = 1.0e-5):
        """returns lower and upper limits of box that should fully contain orbital"""
        # find lower bounds
        p0 = [ np.min([ sh.center[ixyz] - sh.bounding_box_size(thr, lmx) for sh, lmx in zip(self.shells, self.logmxcoeff) ]) for ixyz in range(3) ]

        # find upper bounds
        p1 = [ np.max([ sh.center[ixyz] + sh.bounding_box_size(thr, lmx) for sh, lmx in zip(self.shells, self.logmxcoeff) ]) for ixyz in range(3) ]

        return p0, p1

    def value(self, x, y, z):
        """Returns value of the orbital at specified point"""
        out = 0.0
        for sh, lmx in zip(self.shells, self.logmxcoeff):
            val = sh.value(x, y, z, logmxcoeff = lmx)
            if val is not None:
                out += np.dot(val, self.coeff[sh.start:sh.start+sh.size])
        return out
