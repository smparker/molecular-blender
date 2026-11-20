# -*- coding: utf-8 -*-
#
#  Molecular Blender
#  Filename: molecule.py
#  Copyright (C) 2014-2025 Shane Parker, Joshua Szekely
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

from itertools import product
import math
import numpy as np
import mathutils
from .orbitals import MOData
from .options import DEFAULTS
from .bond_perception import perceive_bond_orders, BondPerceptionError
from .util import stopwatch
from .aromatics import find_planar_cycles
from .periodictable import elements
from .containing_isovalues import isovalue_containing_proportion


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
                 trajectory=None, hidden=False, formal_charge=0):
        """Create Atom"""
        self.element = elements[symbol]
        self.position = mathutils.Vector(position)
        self.index = index
        self.name = name
        self.charge = charge
        self.gradient = mathutils.Vector(gradient if gradient is not None else [0.0, 0.0, 0.0])
        self.trajectory = trajectory if trajectory is not None else []
        self.hidden = hidden
        self.formal_charge = formal_charge

    def __repr__(self):
        """String representation"""
        return "%s (%d keyframes)" % (self.element.symbol, len(self.trajectory))

    @classmethod
    def from_dict(cls, inp):
        """Build Atom from dictionary"""
        return cls(inp["symbol"], inp["position"], inp["index"], charge=inp["charge"],
                   gradient=inp["gradient"], trajectory=[
                       Snapshot.from_dict(x) for x in inp["trajectory"]],
                   hidden=inp["hidden"], formal_charge=inp.get("formal_charge", 0))


class Bond(object):
    """Join two atoms in a bond"""

    cutoff_scale = 1.2

    def __init__(self, iatom, jatom, style="", name="", mapping_function="erf",
                 order=1, direction=None):
        """Create Bond"""
        self.iatom = iatom
        self.jatom = jatom
        self.name = name
        self.bevelname = ""
        # Explicit drawing metadata; bond-order perception is intentionally separate.
        self.order = order
        self.direction = direction
        if style == "":
            self.style = (jatom if iatom.element.vdw >
                          jatom.element.vdw else iatom).element.symbol
        else:
            self.style = style
        self.threshold = self.cutoff_scale * (iatom.element.covalent + jatom.element.covalent)

        if mapping_function == "erf":
            def f(x):
                return math.erfc((x-1)*15*self.threshold)/2
        self.mapping_function = f

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

    def __init__(self, origin, axes, nres, data, datasets=None, labels=None, dataset_index=0, comments=None, field_type="unknown", units=None, linear_combination=None, resample=None):
        """Create Volumetric data storage"""
        self.origin = origin
        self.axes = axes
        self.nres = nres
        self.comments = tuple(comments or ())
        self.field_type = field_type
        self.units = units
        self.linear_combination = linear_combination
        self.resample = resample
        self.datasets = np.asarray(data)[..., None] if datasets is None else np.asarray(datasets)
        self.labels = list(range(self.datasets.shape[-1])) if labels is None else list(labels)
        if self.datasets.shape != tuple(nres) + (len(self.labels),):
            raise ValueError('Cube dataset shape and labels must match the grid')
        self.select_dataset(dataset_index)

        self.axis_norms = np.linalg.norm(self.axes, axis=1).astype(np.float32)
        self.dvol = abs(float(np.linalg.det(np.asarray(self.axes, dtype=float))))

    def select_dataset(self, index):
        """Select a zero-based dataset while preserving all other fields."""
        if isinstance(index, bool) or not isinstance(index, (int, np.integer)) or not 0 <= index < len(self.labels):
            raise ValueError('Cube dataset index is out of range')
        self.dataset_index = int(index)
        self.data = np.ascontiguousarray(self.datasets[..., index])

    def isovalue_containing_proportion(self, values=[0.90], field_type='auto'):
        """Untyped scalar fields integrate magnitude; orbital amplitude integrates its square."""
        kind = self.field_type if field_type == 'auto' else field_type
        if kind not in ('orbital', 'density', 'signed', 'unknown'):
            raise ValueError('Cube field type must be auto, orbital, density, or signed')
        return isovalue_containing_proportion(values, self.data, self.dvol, square=(kind == 'orbital'))

    @classmethod
    def from_dict(cls, inp):
        """Build VolumeData from dictionary"""
        return cls(inp["origin"], inp["axes"], inp["nres"], inp["data"],
                   inp.get('datasets'), inp.get('labels'), inp.get('dataset_index', 0), inp.get('comments'), inp.get('field_type', 'unknown'),
                   inp.get('units'), inp.get('linear_combination'), inp.get('resample'))


_ADJACENT_CELLS = tuple(product((-1, 0, 1), repeat=3))


def bonded_pairs(atoms, static):
    """Return positional pairs in the original i-major order, over any requested frame."""
    radii = [a.element.covalent for a in atoms]
    # Leave a little slack for mathutils' float32 distance at a cell boundary.
    width = max(1., 2 * Bond.cutoff_scale * max(radii, default=0.) * 1.000001)
    count = 0 if static else max((len(a.trajectory) for a in atoms), default=0)
    pairs = set()
    for frame in range(max(count, 1)):
        cells = {}
        for i, atom in enumerate(atoms):
            if count and frame >= len(atom.trajectory):
                continue
            pos = atom.trajectory[frame].position if count else atom.position
            if not all(math.isfinite(v) for v in pos):
                continue
            x, y, z = (math.floor(v / width) for v in pos)
            for dx, dy, dz in _ADJACENT_CELLS:
                for j, previous in cells.get((x+dx, y+dy, z+dz), ()):
                    if (pos - previous).length <= Bond.cutoff_scale * (radii[i] + radii[j]):
                        pairs.add((i, j))
            cells.setdefault((x, y, z), []).append((i, pos))
    return sorted(pairs)


class Molecule(object):
    """Atoms and bonds form a molecule"""

    def __init__(self, name, atoms, volume=None, orbitals=None,
            volume_trajectory=None, orbitals_trajectory=None, bonds=None, frame_numbers=None):
        """Create Molecule"""
        self.name = name
        self.atoms = atoms
        self.bonds = []
        self.explicit_bonds = bonds
        self.rings = []
        self.materials = {}
        self.bond_materials = {}
        self.chgoff = 1.0
        self.chgfac = 1.0
        self.volume = volume
        self.volume_trajectory = volume_trajectory
        self.orbitals = orbitals
        self.orbitals_trajectory = orbitals_trajectory
        self.frame_numbers = frame_numbers

    @classmethod
    def from_dict(cls, name, inp):
        """Build Molecule from dictionary"""
        vol = VolumeData.from_dict(inp["volume"]) if "volume" in inp else None
        vol_traj = [ VolumeData.from_dict(vol) for vol
                in inp["volume_trajectory"] ] if "volume_trajectory" in inp else None
        modata = MOData.from_dict(inp) if "basis" in inp else None
        mo_traj = [ MOData.from_dict(mo) for mo
                in inp["orbital_trajectory" ] ] if "orbital_trajectory" in inp else None
        return cls(name, [Atom.from_dict(x) for x in inp["atoms"]], volume=vol, orbitals=modata,
            volume_trajectory=vol_traj, orbitals_trajectory=mo_traj,
            bonds=inp.get("bonds"), frame_numbers=inp.get("frame_numbers"))

    def center_of_mass(self):
        """Computes center of mass from atoms list"""
        total_mass = 0.0
        out = mathutils.Vector((0.0, 0.0, 0.0))
        for i in self.atoms:
            out += i.position * i.element.mass
            total_mass += i.element.mass
        out /= max(total_mass, 1.0e-10)
        return out

    @stopwatch("determine_bonding")
    def determine_bonding(self, options):
        """builds list of bonded pairs"""
        self.bonds = []
        if not options["bonds"]:
            return

        atoms = self.atoms
        bondstyle = "universal" if options["universal_bonds"] else ""
        if self.explicit_bonds is not None:
            for entry in self.explicit_bonds:
                i, j = entry["atoms"]
                self.bonds.append(Bond(atoms[i], atoms[j], bondstyle, order=entry["order"]))
            self._orient_multiple_bonds()
            if options["find_aromatic"]:
                self.rings = find_planar_cycles(self)
            return
        perceive = options.get("perceive_bond_orders", False)
        if (perceive and options["plot_type"] == "animate"
                and options["animate_bonds"] != "staticfirst"):
            raise BondPerceptionError("Perception requires a static graph; use staticfirst for animations")
        static = options["plot_type"] == "frame" or options["animate_bonds"] == "staticfirst"

        self.bonds = [Bond(atoms[i], atoms[j], bondstyle) for i, j in bonded_pairs(atoms, static)]

        if perceive:
            indices = {id(atom): i for i, atom in enumerate(atoms)}
            edges = [(indices[id(b.iatom)], indices[id(b.jatom)]) for b in self.bonds]
            lengths = [(b.iatom.position - b.jatom.position).length /
                       (b.iatom.element.covalent + b.jatom.element.covalent)
                       for b in self.bonds]
            orders, charges = perceive_bond_orders(
                [a.element.symbol for a in atoms], edges,
                charge=options.get("molecular_charge", 0), bond_lengths=lengths,
                atom_ids=[a.index for a in atoms],
                max_iterations=options.get("bond_order_max_iterations", 100000))
            for bond, order in zip(self.bonds, orders):
                bond.order = order
            for atom, charge in zip(atoms, charges):
                atom.formal_charge = charge
            self._orient_multiple_bonds()

        if options["find_aromatic"]:
            self.rings = find_planar_cycles(self)

    def _orient_multiple_bonds(self):
        """Orient strands using neighbors from the supplied or perceived graph."""
        atoms = self.atoms
        # Use the strongest off-axis neighbor without changing connectivity.
        neighbors = {atom.index: [] for atom in atoms}
        for bond in self.bonds:
            neighbors[bond.iatom.index].append(bond.jatom)
            neighbors[bond.jatom.index].append(bond.iatom)
        for bond in self.bonds:
            if bond.order == 1:
                continue
            axis = bond.jatom.position - bond.iatom.position
            if axis.length < 1e-8:
                raise ValueError("Multiple bond has coincident endpoints")
            axis.normalize()
            candidates = []
            for atom, other in ((bond.iatom, bond.jatom), (bond.jatom, bond.iatom)):
                for neighbor in neighbors[atom.index]:
                    if neighbor is other:
                        continue
                    delta = neighbor.position - atom.position
                    candidates.append(delta - axis * delta.dot(axis))
            if candidates:
                direction = max(candidates, key=lambda v: v.length_squared)
                if direction.length > 1e-8:
                    bond.direction = direction

    def bond_mask(self, options):
        """Construct mask determining whether a bond should be drawn at each frame"""
        outmask = {}
        break_scale = 1 + options.get("bond_break_margin", DEFAULTS["bond_break_margin"]) / 100
        for bond in self.bonds:
            iatom, jatom = bond.iatom, bond.jatom

            pairmask = []
            bonded = False
            break_distance = bond.threshold * break_scale
            for isnap, jsnap in zip(iatom.trajectory, jatom.trajectory):
                vec = isnap.position - jsnap.position
                bonded = vec.length <= (break_distance if bonded else bond.threshold)
                pairmask.append(bonded)

            outmask[(iatom.index, jatom.index)] = pairmask
        return outmask

    def scale(self, chg):
        """Returns sphere scale as a function of charge"""
        thr = 1.0 - self.chgoff
        out = self.chgfac * (chg - thr) + 1.0 if chg > thr else 0.0
        return mathutils.Vector((out, out, out))
