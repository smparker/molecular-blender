# -*- coding: utf-8 -*-
#
#  Molecular Blender
#  Filename: orbitals.py
#  Copyright (C) 2017-2025 Shane Parker, Joshua Szekely
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
from __future__ import annotations

import re
from typing import Sequence
import numpy as np

from .compute_ao import *

"""Classes and functions to compute orbital values in real space."""


def _spin_weights(occupations, spins):
    labels = np.array(spins if spins is not None else ['Alpha'] * len(occupations))
    if np.any(occupations > 1):
        majority = np.minimum(occupations, 1.)
        minority = np.maximum(occupations - 1., 0.)
        return np.where(labels == 'Beta', minority, majority), np.where(labels == 'Beta', majority, minority)
    return np.where(labels == 'Alpha', occupations, 0.), np.where(labels == 'Beta', occupations, 0.)


def sampled_orbital_overlap(previous, current):
    """Estimate normalized signed spatial overlap on a shared midpoint grid in Bohr."""
    lo_a, hi_a = previous.bounding_box()
    lo_b, hi_b = current.bounding_box()
    if np.any(np.minimum(hi_a, hi_b) <= np.maximum(lo_a, lo_b)):
        return 0.0
    lower = np.minimum(lo_a, lo_b)
    span = np.maximum(hi_a, hi_b) - lower
    # Bound the cost independently of the render resolution. Each grid cell
    # has equal volume, which cancels in the normalized inner product.
    counts = np.clip(np.ceil(span / 0.3), 8, 80).astype(int)
    axes = [(origin + (np.arange(n) + 0.5) * width / n).astype(np.float32)
            for origin, width, n in zip(lower, span, counts)]
    left = previous.box_values(*axes).astype(np.float64).ravel()
    right = current.box_values(*axes).astype(np.float64).ravel()
    norm = np.linalg.norm(left) * np.linalg.norm(right)
    return float(np.clip(np.dot(left, right) / norm, -1, 1)) if norm > 0 else 0.0


def frontier_index(occupations: Sequence[float], spins: Sequence[str] | None,
                   energies: Sequence[float | None], kind: str, spin: str | None = None,
                   offset: int = 0) -> int:
    """Resolve a one-based frontier index from orbital metadata alone."""
    if kind not in ('homo', 'lumo') or isinstance(offset, bool) or not isinstance(offset, (int, np.integer)) or offset < 0:
        raise ValueError('Invalid frontier kind or offset')
    if spin not in (None, 'alpha', 'beta'):
        raise ValueError('Orbital spin must be alpha, beta, or omitted')
    weights = np.asarray(occupations)
    capacity = 2. if np.any(weights > 1) else 1.
    indices = list(range(len(weights)))
    if spin is not None:
        if capacity == 1:
            labels = spins or ('Alpha',) * len(weights)
            indices = [i for i in indices if labels[i].lower() == spin]
        weights = _spin_weights(weights, spins)[0 if spin == 'alpha' else 1]
        capacity = 1.
    # Use input order for the whole channel when an energy is missing.
    if all(energies[i] is not None and np.isfinite(energies[i]) for i in indices):
        indices.sort(key=lambda i: (energies[i], i))
    if kind == 'homo':
        candidates = [i for i in reversed(indices) if weights[i] > capacity / 2]
    else:
        candidates = [i for i in indices if weights[i] <= capacity / 2]
    if offset >= len(candidates):
        raise ValueError('Orbital expression {}{} is out of range: no matching {} orbital'.format(
            (spin + ':') if spin else '', kind, 'occupied' if kind == 'homo' else 'empty'))
    return candidates[offset] + 1


class MOData:
    """Organizes data for computing orbitals in real space"""

    def __init__(self, shells, coeff, nocc, occupations, spins=None, densities=None,
                 energies=None, symmetries=None):
        assert sum([sh.size for sh in shells]) == coeff.shape[1]

        self.shells = shells
        self.coeff = coeff
        self.nocc = nocc
        self.nvir = coeff.shape[0] - nocc
        self.occupations, self.spins = self._validate_populations(occupations, spins, coeff.shape[0])
        if isinstance(nocc, (bool, np.bool_)) or not isinstance(nocc, (int, np.integer)) or nocc != np.count_nonzero(self.occupations > 0):
            raise ValueError('Occupied orbital count must match the positive occupations')
        # Metadata follows input orbital order; absent values remain explicit.
        nmo = coeff.shape[0]
        self.energies = tuple([None] * nmo if energies is None else energies)
        self.symmetries = tuple([None] * nmo if symmetries is None else symmetries)
        if len(self.energies) != nmo or len(self.symmetries) != nmo:
            raise ValueError('Provide one energy and symmetry entry per orbital')
        # map of density names to AO densities for things like polarizabilities
        self.densities = {} if densities is None else densities

    @staticmethod
    def _validate_populations(occupations, spins, nmo):
        try:
            occupations = np.asarray(occupations, dtype=np.float64)
        except (TypeError, ValueError) as exc:
            raise ValueError('Orbital occupations must be numeric') from exc
        if occupations.shape != (nmo,):
            raise ValueError('Provide one occupation per orbital')
        if not np.all(np.isfinite(occupations)) or np.any((occupations < 0) | (occupations > 2)):
            raise ValueError('Orbital occupations must be finite and between 0 and 2')
        if spins is None:
            return occupations, None
        try:
            correct_length = len(spins) == nmo
        except TypeError:
            correct_length = False
        if not correct_length:
            raise ValueError('Provide one spin label per orbital, or omit all spin labels')
        labels = []
        for index, label in enumerate(spins):
            if label is None or (isinstance(label, str) and not label.strip()):
                labels.append(None)
            elif isinstance(label, str) and label.strip().lower() in ('alpha', 'beta'):
                labels.append(label.strip().capitalize())
            else:
                raise ValueError('Invalid spin label at orbital {}: {!r}'.format(index+1, label))
        if all(label is None for label in labels):
            return occupations, None
        if any(label is None for label in labels):
            raise ValueError('Spin labels are incomplete; provide Alpha/Beta for every orbital or omit all labels')
        return occupations, tuple(labels)

    def frontier(self, kind='homo', spin=None, offset=0):
        """Return a one-based input index using the half-occupation heuristic."""
        occupations, spins = self.occupations, self.spins
        if spin is not None:
            occupations, spins = self._validate_populations(occupations, spins, self.coeff.shape[0])
        return frontier_index(occupations, spins, self.energies, kind, spin, offset)

    def homo(self, spin=None):
        return self.frontier('homo', spin)

    def lumo(self, spin=None):
        return self.frontier('lumo', spin)

    def resolve_orbital(self, expression: str | int, spin: str | None = None) -> int:
        """Resolve e.g. beta:lumo+1; explicit positive numbers remain input indices."""
        text = str(expression).strip().lower()
        match = re.fullmatch(r'(?:(alpha|beta):)?(homo|lumo)(?:\s*([-+])\s*(\d+))?', text)
        if match:
            channel, kind, sign, offset = match.groups()
            if sign is not None and sign != ('-' if kind == 'homo' else '+'):
                raise ValueError('Use homo-N or lumo+N')
            return self.frontier(kind, channel or spin, int(offset or 0))
        try:
            number = int(text)
        except ValueError as exc:
            raise ValueError('Could not understand orbital {!r}'.format(expression)) from exc
        self._orbital_index(number, spin)
        return number if number > 0 else self.frontier('homo', spin, -number)

    def _orbital_index(self, iorb, spin=None):
        if isinstance(iorb, (bool, np.bool_)) or not isinstance(iorb, (int, np.integer)):
            raise ValueError('Orbital index must be an integer')
        if iorb <= 0:
            if self.nocc == 0:
                raise ValueError('HOMO-relative orbital selection requires an occupied orbital')
            try:
                iorb = self.frontier('homo', spin, -int(iorb))
            except ValueError as exc:
                raise ValueError('Orbital index {} is out of range for occupied orbitals'.format(iorb)) from exc
        if not 1 <= iorb <= self.coeff.shape[0]:
            raise ValueError('Orbital index {} is out of range; available absolute indices are 1..{}'.format(
                iorb, self.coeff.shape[0]))
        return int(iorb)-1

    def get_orbital(self, iorb, spin=None, *, flip_phase=False):
        """Return a Bohr-coordinate evaluator for a one-based orbital or HOMO-relative index."""
        coeff = self.coeff[self._orbital_index(iorb, spin), :]
        return OrbitalCalculator(self.shells, -coeff if flip_phase else coeff)

    @property
    def occupation_scale(self):
        """Visualization heuristic: infer spin (1) or spatial (2) occupations."""
        return 2.0 if np.any(self.occupations > 1) else 1.0

    def spin_populations(self):
        """Conserve actual populations, filling the majority spin first for spatial MOs.

        Spatial orbitals default to alpha majority; a Beta label reverses the
        majority for that orbital. Missing spin labels default to alpha. These
        are visualization defaults, not inferred wavefunction provenance.
        """
        occupations, spins = self._validate_populations(
            self.occupations, self.spins, self.coeff.shape[0])
        return _spin_weights(occupations, spins)

    def density_matrix(self, densityname: str = 'density') -> tuple[np.ndarray, float]:
        """Return the float32 (nao, nao) matrix and its electron population."""
        if not isinstance(densityname, str) or not self.is_density(densityname):
            raise ValueError('Unknown density {!r}; choose from {}'.format(
                densityname, ', '.join(['density', 'alpha', 'beta', 'spin'] + sorted(self.densities))))
        occupations, _ = self._validate_populations(self.occupations, self.spins, self.coeff.shape[0])
        if densityname in self.densities and densityname not in ('density', 'alpha', 'beta', 'spin'):
            return self.densities[densityname], float(np.sum(occupations))
        weights = occupations
        if densityname in ('alpha', 'beta', 'spin'):
            alpha, beta = self.spin_populations()
            weights = {'alpha': alpha, 'beta': beta, 'spin': alpha-beta}[densityname]
        matrix = np.einsum('pu,p,pv->uv', self.coeff, weights, self.coeff).astype(np.float32)
        population = np.sum(occupations) if densityname == 'spin' else np.sum(weights)
        return matrix, float(population)

    def get_density(self, densityname='density'):
        """Return a Bohr-coordinate evaluator for the requested AO density matrix."""
        matrix, population = self.density_matrix(densityname)
        return DensityCalculator(self.shells, matrix, population)

    def scalar_field(self, selection, spin=None, *, flip_phase=False):
        """Resolve a named density or MO to an angstrom-coordinate ScalarField."""
        from .scalar_fields import WavefunctionField
        source = getattr(self, 'source', '')
        if self.is_density(selection):
            return WavefunctionField(self.get_density(selection), source=source, selection=selection)
        index = self.resolve_orbital(selection, spin)
        return WavefunctionField(self.get_orbital(index, flip_phase=flip_phase), source=source,
                                 selection=str(index), orbital=True)

    def is_density(self, densityname):
        return densityname in ["density", "spin", "alpha", "beta"] or densityname in self.densities.keys()

    def response_field(self, eigenpair):
        fields = [name for name in self.densities if name.startswith('response ')]
        if not fields:
            raise ValueError('Response density requires a Molden manifest entry with a DIPL file')
        name = f'response {eigenpair}' if eigenpair else fields[0]
        if name not in fields:
            raise ValueError(f'DIPL eigenpair {eigenpair} is unavailable')
        return name

    @classmethod
    def from_dict(cls, geo):
        """Builds MOData out of the dict object available in MolecularBlender"""
        denorm = geo.get("denorm", True)
        shells = []
        ico = 0
        for atom, ash in zip(geo["atoms"], geo["basis"]):
            # make sure centers are given in bohr
            xyz = [x * ang2bohr for x in atom["position"]]
            for sh in ash:
                newshell = Shell(xyz, "spdfg".index(
                    sh["shell"]), sh["exponents"], sh["contractions"], start=ico,
                                 denorm=denorm)
                ico += newshell.size
                shells.append(newshell)

        occupations, spins = cls._validate_populations(
            [x.get('occup') for x in geo['mo']],
            [x.get('spin') for x in geo['mo']], len(geo['mo']))
        nocc = np.sum(occupations > 0.0)
        nmo = len(occupations)
        nao = sum([sh.size for sh in shells])

        coeff = np.zeros([nmo, nao], dtype=np.float32)
        for i in range(nmo):
            coeff[i, :] = geo["mo"][i]["coeff"]

        densities = {}
        if "polar" in geo:
            # turbomole polarizabilities in occ-virt MO are stored in geo
            nvir = nmo - nocc
            for i, pair in enumerate(geo["polar"]):
                vec = pair["vector"].reshape(nocc, nvir)
                dens = np.einsum("iu,ia,av->uv", coeff[:nocc, :], vec, coeff[nocc:, :])
                matrix = 2.0 * (dens + dens.T)
                densities[f"polar{i:2d}"] = matrix
                densities[f"response {pair.get('index', i + 1)}"] = matrix

        model = cls(shells, coeff, nocc, occupations, spins, densities,
                   energies=[float(orb['ene']) if 'ene' in orb else None for orb in geo['mo']],
                   symmetries=[orb.get('sym') for orb in geo['mo']])
        model.source = geo.get('source', '')
        return model


def orbital_metadata(model):
    if model is None:
        return []
    return [dict(index=i+1, energy=float(e) if e is not None and np.isfinite(e) else None,
                 symmetry=sym, occupation=float(occ), spin=model.spins[i] if model.spins is not None else None)
            for i, (e, sym, occ) in enumerate(zip(model.energies, model.symmetries, model.occupations))]


def response_metadata(data):
    """Keep the DIPL IDs and values for a single geometry, without the response vectors."""
    if data.get('orbital_trajectory'):
        return []
    return [[p['index'], p['eigenvalue']] for p in data.get('polar', ())]
