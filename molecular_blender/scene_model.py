# Copyright (C) 2014-2025 Shane Parker, Joshua Szekely
# SPDX-License-Identifier: GPL-3.0-or-later

"""Blender-facing molecule with presentation names, materials and native vectors."""
import math
from mathutils import Vector

from . import molecular_data as data


class Snapshot(data.Snapshot):
    vector = Vector


class Atom(data.Atom):
    vector = Vector
    snapshot_type = Snapshot

    def __init__(self, symbol, position, index, name="", charge=0.0, gradient=None,
                 trajectory=None, hidden=False, formal_charge=0):
        super().__init__(symbol, position, index, charge, gradient, trajectory, hidden, formal_charge)
        self.name = name


class Bond(data.Bond):
    def __init__(self, iatom, jatom, style="", name="", mapping_function="erf", order=1, direction=None):
        super().__init__(iatom, jatom, order, direction)
        self.name = name
        self.bevelname = ""
        self.style = style or min((iatom, jatom), key=lambda atom: atom.element.vdw).element.symbol
        if mapping_function != "erf":
            raise ValueError("Expected the erf bond mapping function")

    def mapping_function(self, distance):
        return math.erfc((distance - 1) * 15 * self.threshold) / 2

    def make_names(self, basename, split=False):
        """builds name for a molecule from basename and connected atoms"""
        iname = self.iatom.name.split('_')[-1]
        jname = self.jatom.name.split('_')[-1]
        if split:
            return ["%s_%s-%s_a" % (basename, iname, jname),
                    "%s_%s-%s_b" % (basename, iname, jname)]
        return [ basename + "_" + iname + "-" + jname ]


class Molecule(data.Molecule):
    atom_type = Atom
    bond_type = Bond

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.materials = {}
        self.bond_materials = {}
        self.chgoff = self.chgfac = 1.0

    def determine_bonding(self, options):
        super().determine_bonding(options)
        if options["universal_bonds"]:
            for bond in self.bonds:
                bond.style = "universal"

    def center_of_mass(self):
        return Vector(super().center_of_mass())

    def scale(self, chg):
        """Returns sphere scale as a function of charge"""
        thr = 1.0 - self.chgoff
        out = self.chgfac * (chg - thr) + 1.0 if chg > thr else 0.0
        return Vector((out, out, out))
