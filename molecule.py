#
#  Molecular Blender
#  Filename: molecule.py
#  Copyright (C) 2014 Shane Parker, Joshua Szekely
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

from .periodictable import generate_table
from .orbitals import MOData

from .util import stopwatch

import numpy as np
import mathutils

elements = {}  # container for elements


class Snapshot(object):
    def __init__(self, position, charge, gradient):
        self.position = mathutils.Vector(position)
        self.charge = float(charge)
        self.gradient = mathutils.Vector(position)

    @classmethod
    def from_dict(cls, inp):
        return cls(inp["position"], inp["charge"], inp["gradient"])


class Atom(object):
    """Collect information on single atom"""

    def __init__(self, symbol, position, index, name="", charge=0.0, gradient=[0.0, 0.0, 0.0],
                 trajectory=[], hidden=False):
        self.el = elements[symbol]
        self.position = mathutils.Vector(position)
        self.index = index
        self.name = name
        self.charge = charge
        self.gradient = mathutils.Vector(gradient)
        self.trajectory = trajectory
        self.hidden = hidden

    @classmethod
    def from_dict(cls, inp):
        return cls(inp["symbol"], inp["position"], inp["index"], charge=inp["charge"],
                   gradient=inp["gradient"], trajectory=[
                       Snapshot.from_dict(x) for x in inp["trajectory"]],
                   hidden=inp["hidden"])


class Bond(object):
    """Join two atoms in a bond"""

    def __init__(self, iatom, jatom, style="", name=""):
        self.iatom = iatom
        self.jatom = jatom
        self.name = name
        self.bevelname = ""
        if style == "":
            self.style = (jatom if iatom.el.vdw >
                          jatom.el.vdw else iatom).el.symbol
        else:
            self.style = style
        self.threshold = 1.2 * (iatom.el.covalent + jatom.el.covalent)

    def is_bonded(self, distance=None):
        dist = distance if (distance is not None) else (
            self.iatom.position - self.jatom.position).length
        return dist <= self.threshold

    def make_name(self, basename, split=False):
        """builds name for a molecule from basename and connected atoms"""
        iname = self.iatom.name.split('_')[-1]
        jname = self.jatom.name.split('_')[-1]
        if split:
            return ("%s_%s-%s_a" % (basename, iname, jname),
                    "%s_%s-%s_b" % (basename, iname, jname))
        else:
            return basename + "_" + iname + "-" + jname


class VolumeData(object):
    """Stores volumetric information from, for example, a cube file"""

    def __init__(self, origin, axes, nres, data):
        self.origin = origin
        self.axes = axes
        self.nres = nres
        self.data = data

        self.axis_norms = np.array(
            [np.linalg.norm(self.axes[:, i]) for i in range(3)])

    @classmethod
    def from_dict(cls, inp):
        return cls(inp["origin"], inp["axes"], inp["nres"], inp["data"])


class Molecule(object):
    """Atoms and bonds form a molecule"""

    def __init__(self, name, atoms, volume=None, orbitals=None):
        self.name = name
        self.atoms = atoms
        self.bonds = []
        self.materials = {}
        self.bond_materials = {}
        self.chgoff = 1.0
        self.chgfac = 1.0
        self.volume = volume
        self.orbitals = orbitals

    @classmethod
    def from_dict(cls, name, inp):
        vol = VolumeData.from_dict(inp["volume"]) if "volume" in inp else None
        mo = MOData.from_dict(inp) if "basis" in inp else None
        return cls(name, [Atom.from_dict(x) for x in inp["atoms"]], volume=vol, orbitals=mo)

    def COM(self):
        """Computes center of mass from atoms list"""
        total_mass = float(0.0)
        out = mathutils.Vector((0.0, 0.0, 0.0))
        for i in self.atoms:
            out += i.position * i.el.mass
            total_mass += i.el.mass
        out /= total_mass
        return out

    @stopwatch("determine_bonding")
    def determine_bonding(self, options):
        """builds list of bonded pairs"""
        atoms = self.atoms
        natoms = len(atoms)
        bondstyle = "universal" if options["universal_bonds"] else ""
        if (options["bonds"]):
            for i in range(natoms):
                for j in range(i):
                    # make a Bond object no matter what, then use it to check whether to keep
                    bond = Bond(atoms[i], atoms[j], bondstyle)
                    if (options["plot_type"] == "frame" or options["animate_bonds"] == "staticfirst"):
                        # only use the first frame to look for bonded pairs
                        vec = atoms[i].position - atoms[j].position
                        if bond.is_bonded(vec.length):
                            self.bonds.append(bond)
                    elif (options["animate_bonds"] in ["staticall", "dynamic"]):
                        # search through entire trajectory for a bond
                        for va, vb in zip(atoms[i].trajectory, atoms[j].trajectory):
                            if (bond.is_bonded((va.position - vb.position).length)):
                                self.bonds.append(bond)
                                break

    def bond_mask(self, options):
        """Construct mask determining whether a bond should be drawn at each frame"""
        outmask = {}
        natoms = len(self.atoms)
        for bond in self.bonds:
            iatom, jatom = bond.iatom, bond.jatom

            pairmask = []
            for va, vb in zip(iatom.trajectory, jatom.trajectory):
                vec = va.position - vb.position
                pairmask.append(bond.is_bonded(vec.length))

            outmask[(iatom.index, jatom.index)] = pairmask
        return outmask

    def scale(self, chg):
        thr = 1.0 - self.chgoff
        out = self.chgfac * (chg - thr) + 1.0 if chg > thr else 0.0
        return mathutils.Vector((out, out, out))
