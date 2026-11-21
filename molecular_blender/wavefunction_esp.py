"""All-electron molecular electrostatic potential from Cartesian Gaussian products.

McMurchie–Davidson Hermite coefficients contract the total AO density once.
Each point then needs only Boys functions and Coulomb auxiliary integrals.
Internal coordinates are Bohr and potentials are Hartree/e.
"""
from functools import lru_cache
import math

import numpy as np

from .compute_ao_py import Shell
from .constants import ang2bohr
from .orbitals import MOData
from .periodictable import symbols

try:
    from .esp_kernel import coulomb_sum
except ModuleNotFoundError as exc:
    if exc.name != __package__ + '.esp_kernel':
        raise
    coulomb_sum = None

HARTREE_IN_VOLTS = 27.211386245981
POWERS = [[tuple(label.count(axis) for axis in 'xyz') for label in order] for order in Shell.shorder]


def _odd_factorial(l):
    return math.prod(range(1, 2*l, 2))


def nuclear_charges(data, population, charge):
    if isinstance(charge, bool) or not isinstance(charge, (int, np.integer)):
        raise ValueError('Enter the integer net molecular charge')
    if data.get('has_ecp'):
        raise ValueError('Wavefunction ESP requires an all-electron basis. The file declares core or ECP data.')
    charges = []
    for atom in data['atoms']:
        z = atom.get('nuclear_charge')
        if z is None or z <= 0 or z != symbols.index(atom['symbol']):
            raise ValueError('Wavefunction ESP requires full nuclear charges for every atom. Check for ghost centers or ECPs.')
        charges.append(z)
    expected = sum(charges) - charge
    if not charges or population <= 0 or abs(population - expected) > .0001:
        raise ValueError(f'Occupations sum to {population:g} electrons. Net charge {charge:+d} requires '
                         f'{expected:g}. Check the charge, occupied orbitals and all-electron basis.')
    return np.array(charges)


@lru_cache(maxsize=9)
def _indices(order):
    powers = [(x, y, level-x-y) for level in range(order+1) for x in range(level+1)
              for y in range(level-x+1)]
    lookup = {p: i for i, p in enumerate(powers)}
    # For each derivative, store its axis, first and second predecessor, and multiplicity.
    rules = [(0, 0, 0, 0)]
    for powers3 in powers[1:]:
        p = list(powers3)
        axis = next(i for i in range(3) if p[i])
        p[axis] -= 1
        parent, factor = lookup[tuple(p)], p[axis]
        p[axis] -= 1
        rules.append((axis, parent, lookup[tuple(p)] if factor else 0, factor))
    counts = [math.comb(n+3, 3) for n in range(order+1)]
    return lookup, np.asarray(rules, dtype=np.int32), np.asarray(counts, dtype=np.int32)


def _hermite_axis(a, b, shift_a, shift_b, half_inverse):
    # E(i+1,j,t) = E(i,j,t-1)/(2p) + (P-A) E(i,j,t) + (t+1) E(i,j,t+1).
    values = {(0, 0): [np.ones_like(half_inverse)]}
    for i in range(a+1):
        for j in range(b+1):
            if i+j == 0:
                continue
            previous = values[i, j-1] if j else values[i-1, j]
            shift = shift_b if j else shift_a
            row = [v * shift for v in previous] + [np.zeros_like(shift)]
            for t, v in enumerate(previous):
                row[t+1] += v * half_inverse
                if t:
                    row[t-1] += t*v
            values[i, j] = row
    return values


def _primitive_shells(data):
    centers, exps, coeffs, starts, orders = [], [], [], [], []
    start = 0
    for atom, basis in zip(data['atoms'], data['basis']):
        for shell in basis:
            l = 'spdfg'.index(shell['shell'])
            a = np.asarray(shell['exponents'])
            c = np.asarray(shell['contractions'])
            if data['denorm']:
                c = c * (2*a/np.pi)**.75 * np.sqrt((4*a)**l / _odd_factorial(l))
            centers.extend([atom['position']] * len(a))
            exps.extend(a)
            coeffs.extend(c)
            starts.extend([start] * len(a))
            orders.extend([l] * len(a))
            start += len(POWERS[l])
    return np.asarray(centers)*ang2bohr, np.array(exps), np.array(coeffs), np.array(starts), np.array(orders)


def _pair_blocks(data, density, overlap):
    centers, exps, radial, starts, orders = _primitive_shells(data)
    angular = [[math.sqrt(_odd_factorial(l) / math.prod(_odd_factorial(k) for k in p))
                for p in row] for l, row in enumerate(POWERS)]
    for la in sorted(set(orders)):
        shell_a = np.flatnonzero(orders == la)
        for lb in sorted(set(orders[orders <= la])):
            shell_b = np.flatnonzero(orders == lb)
            if la == lb:
                rows, cols = np.tril_indices(len(shell_a))
                first, second = shell_a[rows], shell_b[cols]
            else:
                first, second = np.repeat(shell_a, len(shell_b)), np.tile(shell_b, len(shell_a))
            degree = int(la+lb)
            indices, _, _ = _indices(degree)
            for k in range(0, len(first), 512):
                i, j = first[k:k+512], second[k:k+512]
                a, b = exps[i], exps[j]
                p = a+b
                displacement = centers[j] - centers[i]
                pa = displacement * (b/p)[:, None]
                pb = -displacement * (a/p)[:, None]
                center = centers[i] + pa
                factor = radial[i]*radial[j]*np.where(i == j, 1., 2.)
                factor *= np.exp(-(a*b/p) * np.einsum('ij,ij->i', displacement, displacement))
                axes = [_hermite_axis(la, lb, ax, bx, .5/p) for ax, bx in zip(pa.T, pb.T)]
                volume = (np.pi/p)**1.5
                e = np.zeros((len(indices), len(i)))
                for u, left in enumerate(POWERS[la]):
                    for v, right in enumerate(POWERS[lb]):
                        radial_weight = factor * angular[la][u] * angular[lb][v]
                        row, col = starts[i]+u, starts[j]+v
                        weights = radial_weight * density[row, col]
                        xrow, yrow, zrow = [axis[m, n] for axis, m, n in zip(axes, left, right)]
                        np.add.at(overlap, (row, col), volume * radial_weight * xrow[0]*yrow[0]*zrow[0])
                        for x, ex in enumerate(xrow):
                            wx = weights*ex
                            for y, ey in enumerate(yrow):
                                wxy = wx*ey
                                for z, ez in enumerate(zrow):
                                    e[indices[x, y, z]] += wxy*ez
                # Discard only exact zero products, including floating-point underflow.
                keep = e.any(axis=0)
                if np.any(keep):
                    yield degree, p[keep], center[keep], np.ascontiguousarray((e[:, keep] * (2*np.pi/p[keep])).T)


@lru_cache(maxsize=1)
def _boys_table():
    # F_n(T) = integral_0^1 t^(2n) exp(-T t^2) dt. Six Taylor corrections
    # around the nearest 1/16 grid point give double precision at small T.
    x, w = np.polynomial.legendre.leggauss(48)
    t2 = ((x+1)/2)**2
    return np.ascontiguousarray(np.exp(-np.arange(641)[:, None]/16*t2) @
                                ((w/2)[:, None] * t2[:, None]**np.arange(15)))


def _numpy_sum(points, degree, p, centers, coeff, rules, counts, table):
    # Portable fallback. Vectorize across primitive pairs to bound temporary memory.
    result = np.empty(len(points))
    for k, point in enumerate(points):
        delta = centers - point
        t = p * np.einsum('ij,ij->i', delta, delta)
        small = t < 40.
        grid = np.rint(t[small]*16).astype(int)
        dx = grid/16 - t[small]
        boys = np.empty((degree+1, len(p)))
        for n in range(degree+1):
            v = table[grid, n+6]
            for h in range(6, 0, -1):
                v = table[grid, n+h-1] + dx*v/h
            boys[n, small] = v
            if n:
                boys[n, ~small] = boys[n-1, ~small] * ((2*n-1)/2) / t[~small]
            else:
                boys[n, ~small] = .5*np.sqrt(np.pi/t[~small])
        aux = np.empty((len(rules), len(p)))
        for n in range(degree, -1, -1):
            for idx in range(counts[degree-n]-1, 0, -1):
                axis, prev, prev2, factor = rules[idx]
                aux[idx] = delta[:, axis]*aux[prev]
                if factor:
                    aux[idx] += factor*aux[prev2]
            aux[0] = (-2*p)**n * boys[n]
        result[k] = np.einsum('ij,ji', coeff, aux)
    return result


class WavefunctionESP:
    def __init__(self, data, charge=0, model=None):
        if not data.get('basis') or not data.get('mo'):
            raise ValueError('Wavefunction ESP needs a Molden basis and occupied orbitals')
        model = model if model is not None else MOData.from_dict(data)
        density, population = model.density_matrix('density')
        self.charges = nuclear_charges(data, population, charge)
        self.nuclei = np.array([atom['position'] for atom in data['atoms']])*ang2bohr
        overlap = np.zeros_like(density, dtype=float)
        self.blocks = list(_pair_blocks(data, density, overlap))
        # Only the zeroth Hermite term has a nonzero integral over all space.
        integrated = sum(np.sum(c[:, 0] * .5 * np.sqrt(np.pi/p)) for _, p, _, c in self.blocks)
        if not math.isclose(integrated, population, rel_tol=2e-5, abs_tol=2e-5):
            raise ValueError(f'Integrated density gives {integrated:.6f} electrons, expected {population:g}. '
                             'Check the Molden normalization convention and occupied orbitals.')

        overlap = (overlap + overlap.T)*.5
        spins = np.asarray(model.spins or ['Alpha']*len(model.occupations))
        for spin in np.unique(spins):
            occupied = model.coeff[(model.occupations > 0) & (spins == spin)]
            gram = occupied @ overlap @ occupied.T
            if not np.allclose(gram, np.eye(len(occupied)), atol=2e-4, rtol=0):
                raise ValueError(f'Occupied {spin} orbitals must be normalized and mutually orthogonal. '
                                 'Check the wavefunction and its Molden convention.')

    def potential_bohr(self, points):
        points = np.ascontiguousarray(points, dtype=float)
        if points.ndim != 2 or points.shape[1] != 3 or not np.isfinite(points).all():
            raise ValueError('ESP needs finite (N, 3) point coordinates')
        r = np.linalg.norm(points[:, None, :] - self.nuclei, axis=2)
        if np.any(r < 1e-8):
            raise ValueError('ESP is singular at a nucleus (within 1e-8 bohr)')
        values = (self.charges / r).sum(axis=1)
        kernel = coulomb_sum or _numpy_sum
        for degree, p, center, coeff in self.blocks:
            _, rules, counts = _indices(degree)
            values -= kernel(points, degree, p, center, coeff, rules, counts, _boys_table())
        return values

    def batch(self, points):
        """Potential in Hartree/e on source-coordinate vertices in angstroms."""
        return self.potential_bohr(np.asarray(points)*ang2bohr)
