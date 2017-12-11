# -*- coding: utf-8 -*-
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

"""Classes to encapsulate molecule data"""

import numpy as np
import mathutils
from .orbitals import MOData
from .util import stopwatch
from .aromatics import find_planar_cycles
from .periodictable import elements


class Snapshot(object):
    """Single snapshot of atomic position, charge, gradient during trajectory"""

    def __init__(self, position, charge, gradient):
        """Create Snapshot from position, charge, gradient"""
        self.position = mathutils.Vector(position)
        self.charge = float(charge)
        self.gradient = mathutils.Vector(gradient)

    @classmethod
    def from_dict(cls, inp):
        """Build Snapshot from dictionary"""
        return cls(inp["position"], inp["charge"], inp["gradient"])


class Atom(object):
    """Collect information on single atom"""

    def __init__(self, symbol, position, index, name="", charge=0.0, gradient=None,
                 trajectory=None, hidden=False):
        """Create Atom"""
        self.element = elements[symbol]
        self.position = mathutils.Vector(position)
        self.index = index
        self.name = name
        self.charge = charge
        self.gradient = mathutils.Vector(gradient if gradient is not None else [0.0, 0.0, 0.0])
        self.trajectory = trajectory if trajectory is not None else []
        self.hidden = hidden

    def __repr__(self):
        """String representation"""
        return "%s (%d keyframes)" % (self.element.symbol, len(self.trajectory))

    @classmethod
    def from_dict(cls, inp):
        """Build Atom from dictionary"""
        return cls(inp["symbol"], inp["position"], inp["index"], charge=inp["charge"],
                   gradient=inp["gradient"], trajectory=[
                       Snapshot.from_dict(x) for x in inp["trajectory"]],
                   hidden=inp["hidden"])


class Bond(object):
    """Join two atoms in a bond"""

    def __init__(self, iatom, jatom, style="", name=""):
        """Create Bond"""
        self.iatom = iatom
        self.jatom = jatom
        self.name = name
        self.bevelname = ""
        if style == "":
            self.style = (jatom if iatom.element.vdw >
                          jatom.element.vdw else iatom).element.symbol
        else:
            self.style = style
        self.threshold = 1.2 * (iatom.element.covalent + jatom.element.covalent)

    def is_bonded(self, distance=None):
        """Returns whether the pair of atoms are bonded at given distance"""
        dist = distance if (distance is not None) else (
            self.iatom.position - self.jatom.position).length
        return dist <= self.threshold

    def make_names(self, basename, split=False):
        """builds name for a molecule from basename and connected atoms"""
        iname = self.iatom.name.split('_')[-1]
        jname = self.jatom.name.split('_')[-1]
        if split:
            return ["%s_%s-%s_a" % (basename, iname, jname),
                    "%s_%s-%s_b" % (basename, iname, jname)]
        return [ basename + "_" + iname + "-" + jname ]


class VolumeData(object):
    """Stores volumetric information from, for example, a cube file"""

    def __init__(self, origin, axes, nres, data):
        """Create Volumetric data storage"""
        self.origin = origin
        self.axes = axes
        self.nres = nres
        self.data = data

        self.axis_norms = np.array(
            [np.linalg.norm(self.axes[:, i]) for i in range(3)])

    @classmethod
    def from_dict(cls, inp):
        """Build VolumeData from dictionary"""
        return cls(inp["origin"], inp["axes"], inp["nres"], inp["data"])


class Molecule(object):
    """Atoms and bonds form a molecule"""

    def __init__(self, name, atoms, volume=None, orbitals=None):
        """Create Molecule"""
        self.name = name
        self.atoms = atoms
        self.bonds = []
        self.rings = []
        self.materials = {}
        self.bond_materials = {}
        self.chgoff = 1.0
        self.chgfac = 1.0
        self.volume = volume
        self.orbitals = orbitals

    @classmethod
    def from_dict(cls, name, inp):
        """Build Molecule from dictionary"""
        vol = VolumeData.from_dict(inp["volume"]) if "volume" in inp else None
        modata = MOData.from_dict(inp) if "basis" in inp else None
        return cls(name, [Atom.from_dict(x) for x in inp["atoms"]], volume=vol, orbitals=modata)

    def center_of_mass(self):
        """Computes center of mass from atoms list"""
        total_mass = 0.0
        out = mathutils.Vector((0.0, 0.0, 0.0))
        for i in self.atoms:
            out += i.position * i.element.mass
            total_mass += i.element.mass
        out /= total_mass
        return out

    @stopwatch("determine_bonding")
    def determine_bonding(self, options):
        """builds list of bonded pairs"""
        if not options["bonds"]:
            return

        atoms = self.atoms
        natoms = len(atoms)
        bondstyle = "universal" if options["universal_bonds"] else ""
        static = options["plot_type"] == "frame" or options["animate_bonds"] == "staticfirst"

        for i in range(natoms):
            for j in range(i):
                # make a Bond object no matter what, then use it to check whether to keep
                bond = Bond(atoms[i], atoms[j], bondstyle)
                distance = (atoms[i].position - atoms[j].position).length if static \
                    else min([(iatom.position - jatom.position).length for \
                                    iatom, jatom in zip(atoms[i].trajectory, atoms[j].trajectory)])
                if bond.is_bonded(distance):
                    self.bonds.append(bond)

        if options["find_aromatic"]:
            self.rings = find_planar_cycles(self)

    def bond_mask(self, options):
        """Construct mask determining whether a bond should be drawn at each frame"""
        outmask = {}
        for bond in self.bonds:
            iatom, jatom = bond.iatom, bond.jatom

            pairmask = []
            for isnap, jsnap in zip(iatom.trajectory, jatom.trajectory):
                vec = isnap.position - jsnap.position
                pairmask.append(bond.is_bonded(vec.length))

            outmask[(iatom.index, jatom.index)] = pairmask
        return outmask

    def scale(self, chg):
        """Returns sphere scale as a function of charge"""
        thr = 1.0 - self.chgoff
        out = self.chgfac * (chg - thr) + 1.0 if chg > thr else 0.0
        return mathutils.Vector((out, out, out))
