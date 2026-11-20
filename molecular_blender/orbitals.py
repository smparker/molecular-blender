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
import re
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


def frontier_index(occupations, spins, energies, kind, spin=None, offset=0):
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

    def resolve_orbital(self, expression, spin=None):
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

    def get_orbital(self, iorb, spin=None):
        return OrbitalCalculater(self.shells, self.coeff[self._orbital_index(iorb, spin), :])

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

    def density_matrix(self, densityname='density'):
        """Return the AO matrix and population for a known field."""
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
        matrix, population = self.density_matrix(densityname)
        return DensityCalculater(self.shells, matrix, population)

    def is_density(self, densityname):
        return densityname in ["density", "spin", "alpha", "beta"] or densityname in self.densities.keys()

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
                eigenvalue = pair["eigenvalue"]
                if False:
                    vec = pair["vector"].reshape(nvir, nocc)
                    dens = np.einsum("au,ai,iv->uv", coeff[nocc:, :], vec, coeff[:nocc, :])
                    densities[f"polar{i:2d}"] = 2.0 * (dens + dens.T)
                else:
                    vec = pair["vector"].reshape(nocc, nvir)
                    dens = np.einsum("iu,ia,av->uv", coeff[:nocc, :], vec, coeff[nocc:, :])
                    densities[f"polar{i:2d}"] = 2.0*(dens + dens.T)

        return cls(shells, coeff, nocc, occupations, spins, densities,
                   energies=[float(orb['ene']) if 'ene' in orb else None for orb in geo['mo']],
                   symmetries=[orb.get('sym') for orb in geo['mo']])
