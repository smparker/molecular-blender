# -*- coding: utf-8 -*-
#
#  Molecular Blender
#  Filename: molecular_data.py
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

"""Scientific model. Positions are in angstroms. Atom.index is a source ID.

Explicit edges and ring cycles address zero-based positions in Molecule.atoms.
"""

from __future__ import annotations

from itertools import product
import math
import numpy as np
from typing import Sequence, TYPE_CHECKING

from .bond_perception import AtomOffset, AtomId

if TYPE_CHECKING:
    from numpy.typing import ArrayLike, NDArray
from .orbitals import MOData
from .options import DEFAULTS
from .bond_perception import perceive_bond_orders, BondPerceptionError
from .util import stopwatch
from .aromatics import find_planar_cycles
from .periodictable import elements
from .containing_isovalues import isovalue_containing_proportion


def coordinate(values: ArrayLike) -> NDArray[np.float32]:
    """Source-space xyz triple in angstroms, matching stored Blender precision."""
    xyz = np.asarray(values, dtype=np.float32)
    if xyz.shape != (3,):
        raise ValueError("Expected three Cartesian components")
    return xyz.copy()


class Snapshot(object):
    """Single snapshot of atomic position, charge, gradient during trajectory"""

    vector = staticmethod(coordinate)

    def __init__(self, position, charge, gradient):
        """Create Snapshot from position, charge, gradient"""
        self.position = self.vector(position)
        self.charge = float(charge)
        self.gradient = self.vector(gradient)

    @classmethod
    def from_dict(cls, inp):
        """Build Snapshot from dictionary"""
        return cls(inp["position"], inp["charge"], inp["gradient"])


class Atom(object):
    """Collect information on single atom"""

    snapshot_type = Snapshot
    vector = staticmethod(coordinate)

    def __init__(self, symbol: str, position: ArrayLike, index: AtomId, charge=0.0, gradient=None,
                 trajectory=None, hidden=False, formal_charge=0):
        """Create Atom"""
        self.element = elements[symbol]
        self.position = self.vector(position)
        self.index = AtomId(index)
        self.charge = charge
        self.gradient = self.vector(gradient if gradient is not None else [0.0, 0.0, 0.0])
        self.trajectory = trajectory if trajectory is not None else []
        self.hidden = hidden
        self.formal_charge = formal_charge
        self.frame_charges = ()

    def __repr__(self):
        """String representation"""
        return "%s (%d keyframes)" % (self.element.symbol, len(self.trajectory))

    @classmethod
    def from_dict(cls, inp):
        """Build Atom from dictionary"""
        return cls(inp["symbol"], inp["position"], inp["index"], charge=inp["charge"],
                   gradient=inp["gradient"], trajectory=[
                       cls.snapshot_type.from_dict(x) for x in inp["trajectory"]],
                   hidden=inp["hidden"], formal_charge=inp.get("formal_charge", 0))


class Bond(object):
    """Join two atoms in a bond"""

    cutoff_scale = 1.2

    def __init__(self, iatom: Atom, jatom: Atom, order: int = 1, direction=None):
        self.iatom = iatom
        self.jatom = jatom
        self.order = order
        self.frame_orders = ()
        self.direction = direction
        self.threshold = self.cutoff_scale * (iatom.element.covalent + jatom.element.covalent)

    def is_bonded(self, distance=None):
        """Returns whether the pair of atoms are bonded at given distance"""
        dist = distance if (distance is not None) else math.dist(
            self.iatom.position, self.jatom.position)
        return dist <= self.threshold

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

    def scalar_field(self, field_type="auto"):
        """Snapshot the selected dataset and its provenance without copying the grid."""
        from .scalar_fields import CubeField
        grid = dict(origin=self.origin, axes=self.axes, data=self.data, nres=self.nres,
                    labels=self.labels, dataset_index=self.dataset_index, units=self.units,
                    field_type=self.field_type, source=getattr(self, 'source', ''))
        return CubeField(grid, field_type=field_type)

    def select_dataset(self, index: int) -> None:
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
        volume = cls(inp["origin"], inp["axes"], inp["nres"], inp["data"],
                   inp.get('datasets'), inp.get('labels'), inp.get('dataset_index', 0), inp.get('comments'), inp.get('field_type', 'unknown'),
                   inp.get('units'), inp.get('linear_combination'), inp.get('resample'))
        volume.source = inp.get('source', '')
        return volume


_ADJACENT_CELLS = tuple(product((-1, 0, 1), repeat=3))


def bonded_pairs(atoms: Sequence[Atom], static: bool) -> list[tuple[AtomOffset, AtomOffset]]:
    """Return positional pairs in the original i-major order, over any requested frame."""
    radii = [a.element.covalent for a in atoms]
    # Leave slack for float32 coordinates at cell boundaries.
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
                    if math.dist(pos, previous) <= Bond.cutoff_scale * (radii[i] + radii[j]):
                        pairs.add((AtomOffset(i), AtomOffset(j)))
            cells.setdefault((x, y, z), []).append((i, pos))
    return sorted(pairs)


class Molecule(object):
    """Atoms and bonds form a molecule"""

    atom_type = Atom
    bond_type = Bond

    def __init__(self, name, atoms, volume=None, orbitals=None,
            volume_trajectory=None, orbitals_trajectory=None, bonds=None, frame_numbers=None):
        """Create Molecule"""
        self.name = name
        self.atoms = atoms
        if len({atom.index for atom in atoms}) != len(atoms):
            raise ValueError("Atom source IDs must be unique within a molecule")
        self.bonds = []
        self.explicit_bonds = bonds
        self.rings = []
        self.volume = volume
        self.volume_trajectory = volume_trajectory
        self.orbitals = orbitals
        self.orbitals_trajectory = orbitals_trajectory
        self.frame_numbers = frame_numbers
        self.perception = []

    @classmethod
    def from_dict(cls, name, inp):
        """Build Molecule from dictionary"""
        vol = VolumeData.from_dict(inp["volume"]) if "volume" in inp else None
        vol_traj = [ VolumeData.from_dict(vol) for vol
                in inp["volume_trajectory"] ] if "volume_trajectory" in inp else None
        modata = MOData.from_dict(inp) if inp.get("mo") else None
        mo_traj = [ MOData.from_dict(mo) for mo
                in inp["orbital_trajectory" ] ] if "orbital_trajectory" in inp else None
        return cls(name, [cls.atom_type.from_dict(x) for x in inp["atoms"]], volume=vol, orbitals=modata,
            volume_trajectory=vol_traj, orbitals_trajectory=mo_traj,
            bonds=inp.get("bonds"), frame_numbers=inp.get("frame_numbers"))

    def center_of_mass(self) -> NDArray:
        """Computes center of mass from atoms list"""
        total_mass = 0.0
        out = np.zeros(3, dtype=np.float64)
        for i in self.atoms:
            out += np.asarray(i.position) * i.element.mass
            total_mass += i.element.mass
        out /= max(total_mass, 1.0e-10)
        return out

    @stopwatch("determine_bonding")
    def determine_bonding(self, options: dict) -> None:
        """builds list of bonded pairs"""
        atoms = self.atoms
        if not options["bonds"]:
            self.bonds, self.rings, self.perception = [], [], []
            for atom in atoms:
                atom.frame_charges = ()
            return

        if self.explicit_bonds is not None:
            explicit = []
            for entry in self.explicit_bonds:
                i, j = entry["atoms"]
                if not (0 <= i < len(atoms) and 0 <= j < len(atoms)) or i == j:
                    raise ValueError(f"Invalid atom offsets {i}, {j} in explicit bond")
                explicit.append(self.bond_type(atoms[i], atoms[j], order=entry["order"]))
            self.bonds, self.rings, self.perception = explicit, [], []
            for atom in atoms:
                atom.frame_charges = ()
            self._orient_multiple_bonds()
            if options["find_aromatic"]:
                self.rings = find_planar_cycles(self)
            return
        perceive = options.get("perceive_bond_orders", False)
        animated = options["plot_type"] == "animate"
        mode = options["animate_bonds"]
        if perceive and animated and mode == "staticall":
            raise BondPerceptionError("Choose dynamic bonds for per-frame perception or staticfirst for one fixed graph")
        edges = bonded_pairs(atoms, not animated or mode == "staticfirst")
        bonds = [self.bond_type(atoms[i], atoms[j]) for i, j in edges]

        reports = []
        if perceive:
            dynamic = animated and mode == "dynamic"
            if dynamic and (not atoms or not atoms[0].trajectory or
                            any(len(a.trajectory) != len(atoms[0].trajectory) for a in atoms)):
                raise BondPerceptionError('Dynamic perception needs aligned trajectory snapshots for every atom')
            masks = self.bond_mask(options, bonds) if dynamic else {}
            symbols = [a.element.symbol for a in atoms]
            identifiers = [a.index for a in atoms]
            radii = [bond.iatom.element.covalent + bond.jatom.element.covalent for bond in bonds]
            samples, charges_by_frame = [], []
            for frame in range(len(atoms[0].trajectory) if dynamic else 1):
                active = [k for k, bond in enumerate(bonds) if not dynamic or
                          masks[(bond.iatom.index, bond.jatom.index)][frame]]
                positions = [a.trajectory[frame].position if dynamic else a.position for a in atoms]
                lengths = [math.dist(positions[edges[k][0]], positions[edges[k][1]]) / radii[k] for k in active]
                report = {}
                try:
                    orders, charges = perceive_bond_orders(symbols, [edges[k] for k in active],
                        charge=options.get("molecular_charge", 0), bond_lengths=lengths, atom_ids=identifiers,
                        max_iterations=options.get("bond_order_max_iterations", 100000), diagnostics=report)
                except BondPerceptionError as exc:
                    if not dynamic:
                        raise
                    source_frame = self.frame_numbers[frame] if self.frame_numbers else frame
                    raise BondPerceptionError(f'Input frame {source_frame}: {exc}', reason=exc.reason,
                                              diagnostics=dict(exc.diagnostics, frame=source_frame)) from exc
                frame_orders = [0] * len(bonds)
                for k, order in zip(active, orders):
                    frame_orders[k] = order
                samples.append(frame_orders)
                charges_by_frame.append(charges)
                reports.append(report)
            for k, bond in enumerate(bonds):
                if dynamic:
                    bond.frame_orders = tuple(sample[k] for sample in samples)
                bond.order = next(sample[k] for sample in samples if sample[k])
            for i, atom in enumerate(atoms):
                atom.formal_charge = charges_by_frame[0][i]
                atom.frame_charges = tuple(sample[i] for sample in charges_by_frame) if dynamic else ()
        else:
            for atom in atoms:
                atom.frame_charges = ()

        self.bonds, self.rings, self.perception = bonds, [], reports
        if perceive:
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
            if max(bond.frame_orders, default=bond.order) <= 1:
                continue
            axis = np.asarray(bond.jatom.position) - bond.iatom.position
            length = np.linalg.norm(axis)
            if length < 1e-8:
                if bond.frame_orders:
                    continue
                raise ValueError(f"Multiple bond has coincident endpoints at atoms {bond.iatom.index}, {bond.jatom.index}")
            axis /= length
            candidates = []
            for atom, other in ((bond.iatom, bond.jatom), (bond.jatom, bond.iatom)):
                for neighbor in neighbors[atom.index]:
                    if neighbor is other:
                        continue
                    delta = np.asarray(neighbor.position) - atom.position
                    candidates.append(delta - axis * delta.dot(axis))
            if candidates:
                direction = max(candidates, key=lambda v: v.dot(v))
                if np.linalg.norm(direction) > 1e-8:
                    bond.direction = self.atom_type.vector(direction)

    def bond_mask(self, options: dict, bonds=None) -> dict[tuple[AtomId, AtomId], list[bool]]:
        """Construct mask determining whether a bond should be drawn at each frame"""
        outmask = {}
        break_scale = 1 + options.get("bond_break_margin", DEFAULTS["bond_break_margin"]) / 100
        for bond in self.bonds if bonds is None else bonds:
            iatom, jatom = bond.iatom, bond.jatom

            pairmask = []
            bonded = False
            break_distance = bond.threshold * break_scale
            for isnap, jsnap in zip(iatom.trajectory, jatom.trajectory):
                bonded = math.dist(isnap.position, jsnap.position) <= (break_distance if bonded else bond.threshold)
                pairmask.append(bonded)

            outmask[(iatom.index, jatom.index)] = pairmask
        return outmask
