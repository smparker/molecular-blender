"""Valence-constrained bond orders for an existing distance-based graph.

Inspired by RDKit DetermineBonds / xyz2mol's separation of connectivity, allowed
valences, unsaturation and charge validation. This is an independent constrained
search, not a port of their maximum-matching implementation:
https://github.com/rdkit/rdkit/blob/master/Code/GraphMol/DetermineBonds/DetermineBonds.cpp
https://github.com/jensengroup/xyz2mol

All atoms (including hydrogens) must be supplied. No radical, metal coordination,
stereochemistry or aromaticity model is inferred. Rings get integer Kekule orders.
"""
from __future__ import annotations

import math
from time import perf_counter
from typing import Hashable, Iterable, Mapping, NewType, Sequence

AtomId = NewType("AtomId", int)  # A source-file label, potentially sparse or one-based.
AtomOffset = NewType("AtomOffset", int)  # Position in the input atom sequence.


# (bond-order sum, formal charge). Main-group closed-shell states, including
# common electron-deficient B/C cations and expanded-valence P/S exceptions.
VALENCE_STATES = {
    'h': ((1, 0),),
    'b': ((3, 0), (4, -1)),
    'c': ((4, 0), (3, -1), (3, 1)),
    'n': ((3, 0), (4, 1), (2, -1)),
    'o': ((2, 0), (1, -1), (3, 1)),
    'f': ((1, 0), (0, -1)),
    'si': ((4, 0),),
    'p': ((3, 0), (5, 0), (4, 1)),
    's': ((2, 0), (4, 0), (6, 0), (1, -1), (3, 1)),
    'cl': ((1, 0), (0, -1)),
    'ge': ((4, 0),),
    'br': ((1, 0), (0, -1)),
    'i': ((1, 0), (0, -1)),
}


class BondPerceptionError(ValueError):
    """No supported, fully validated assignment could be obtained."""

    def __init__(self, message, *, reason=None, diagnostics=None, source=None):
        super().__init__(f"{source}: Bond perception: {message}" if source else message)
        self.path = source
        self.reason = reason
        self.diagnostics = {} if diagnostics is None else diagnostics


def perceive_bond_orders(symbols: Sequence[str], edges: Iterable[tuple[AtomOffset, AtomOffset]],
                         charge: int = 0, bond_lengths: Sequence[float] | None = None,
                         max_iterations: int = 100000, atom_ids: Sequence[Hashable] | None = None,
                         fixed_orders: Mapping[int, int] | None = None,
                         diagnostics: dict | None = None) -> tuple[list[int], list[int]]:
    """Return (orders in input edge order, formal charges in atom order).

    Every edge starts with order one and can rise to at most three. Propagation
    removes orders incompatible with the endpoint valences. Backtracking tries
    the most constrained edges first, using shorter bonds to break ties. Search
    uses increasing budgets for total absolute formal charge, so a neutral
    assignment wins over unnecessary charge separation. Equal-cost resonance
    forms are heuristic choices, not unique chemical determinations.

    Neither the inputs nor a molecule are mutated, including on failure. The
    node budget is shared across fragments and charge budgets. Disconnected
    graphs are solved separately and their charges combined to meet the total.
    Supply a diagnostics dict for elapsed seconds, nodes, fragment count,
    minimum charge cost and uniqueness. If the remaining nodes cannot settle
    uniqueness, the completed assignment is returned with ambiguity undetermined.
    bond_lengths may be normalized by covalent radii by the caller.
    atom_ids optionally supplies stable caller IDs for diagnostics; edges still
    use positional indices. fixed_orders maps input edge indices to locked
    integer orders (1..3); incompatible locks fail without a partial assignment.
    Errors distinguish supported-chemistry limitations,
    impossible connectivity, incompatible charge, and exhausted search.
    """
    started = perf_counter()
    if not isinstance(charge, int) or isinstance(charge, bool):
        raise ValueError('Molecular charge must be an integer')
    if (not isinstance(max_iterations, int) or isinstance(max_iterations, bool)
            or max_iterations < 1):
        raise ValueError('max_iterations must be a positive integer')
    symbols = [s.lower() for s in symbols]
    n = len(symbols)
    identifiers = list(range(n)) if atom_ids is None else list(atom_ids)
    if len(identifiers) != n or len(set(identifiers)) != n:
        raise ValueError('Provide one unique atom ID per symbol')
    iterations = 0

    def failure(reason, message, indices=None):
        indices = list(range(n)) if indices is None else list(indices)
        atoms = [dict(index=i, id=identifiers[i], element=symbols[i],
                      degree=len(incident[i]),
                      candidate_states=VALENCE_STATES.get(symbols[i], ())) for i in indices]
        details = '; '.join('atom {id} ({element}), degree {degree}, valence/charge states {candidate_states}'.format(**a)
                            for a in atoms[:8])
        if len(atoms) > 8:
            details += '; and {} more atoms'.format(len(atoms) - 8)
        diagnostics = dict(atoms=atoms, requested_charge=charge,
                           iterations=iterations, max_iterations=max_iterations)
        return BondPerceptionError(
            '{} Requested charge {}; search nodes {}/{}. {}'.format(
                message, charge, iterations, max_iterations, details),
            reason=reason, diagnostics=diagnostics)

    incident = [[] for _ in symbols]
    edges = list(edges)
    seen = set()
    for k, (i, j) in enumerate(edges):
        if (not isinstance(i, int) or not isinstance(j, int)
                or not 0 <= i < n or not 0 <= j < n or i == j):
            raise ValueError('Invalid bond endpoints')
        pair = tuple(sorted((i, j)))
        if pair in seen:
            raise ValueError('Duplicate bond')
        seen.add(pair)
        incident[i].append(k)
        incident[j].append(k)
    lengths = list(bond_lengths) if bond_lengths is not None else [1.0] * len(edges)
    if len(lengths) != len(edges) or any(not math.isfinite(x) or x <= 0 for x in lengths):
        raise ValueError('Provide one finite positive length per bond')

    states = []
    for i, symbol in enumerate(symbols):
        if symbol not in VALENCE_STATES:
            raise failure('unsupported_element', 'Unsupported element for bond-order perception.', [i])
        states.append(VALENCE_STATES[symbol])

    for i, bonds in enumerate(incident):
        degree = len(bonds)
        if degree > max(v for v, _ in states[i]):
            raise failure('bad_connectivity', 'No valid assignment: degree exceeds every supported valence; check connectivity.', [i])
        if 3 * degree < min(v for v, _ in states[i]):
            raise failure('insufficient_connectivity',
                          'No valid assignment: too few neighbors for a supported valence; check missing hydrogens and connectivity.', [i])

    # Each fragment owns disjoint edge variables. Only its total charge
    # couples it to the rest of the molecule.
    fragments = []
    visited = set()
    for atom in range(n):
        if atom in visited:
            continue
        stack = [atom]
        visited.add(atom)
        members, indices = [], set()
        while stack:
            i = stack.pop()
            members.append(i)
            for k in incident[i]:
                indices.add(k)
                a, b = edges[k]
                other = b if a == i else a
                if other not in visited:
                    stack.append(other)
                    visited.add(other)
        fragments.append((sorted(members), sorted(indices)))

    def propagate(domains, members, budget):
        """Prune by endpoint valence capacity until no further order changes."""
        while True:
            changed = False
            possible = []
            for i in members:
                bonds = incident[i]
                low = sum(domains[k][0] for k in bonds)
                high = sum(domains[k][-1] for k in bonds)
                choices = [(v, q) for v, q in states[i] if low <= v <= high and abs(q) <= budget]
                if not choices:
                    return None
                possible.append(choices)
                for k in bonds:
                    others_low = low - domains[k][0]
                    others_high = high - domains[k][-1]
                    keep = tuple(order for order in domains[k] if any(
                        others_low + order <= v <= others_high + order for v, _ in choices))
                    if not keep:
                        return None
                    if keep != domains[k]:
                        domains[k] = keep
                        changed = True
            if not changed:
                return possible

    def spend_node():
        nonlocal iterations
        if iterations >= max_iterations:
            raise failure('search_limit', 'Bond-order search limit reached.')
        iterations += 1

    # Retain a representative minimum-cost path and whether a second path
    # ties it. At atom level a path assigns charges, at fragment level it
    # selects a solution for each fragment.
    def extend(table, total, cost, path, count):
        previous = table.get(total)
        if previous is None or cost < previous[0]:
            table[total] = (cost, path, count)
        elif cost == previous[0]:
            table[total] = (cost, previous[1], min(2, count + previous[2]))

    def charge_paths(possible, budget, keep_path):
        totals = {0: (0, (), 1)}
        for choices in possible:
            following = {}
            for q in sorted({q for _, q in choices}, key=lambda q: (abs(q), q)):
                for total, (cost, path, count) in totals.items():
                    if cost + abs(q) <= budget:
                        extend(following, total + q, cost + abs(q),
                               path + (q,) if keep_path else (), count)
            totals = following
        return totals

    initial = {k: (1, 2, 3) for k in range(len(edges))}
    for index, order in (fixed_orders or {}).items():
        if (isinstance(index, bool) or not isinstance(index, int) or not 0 <= index < len(edges)
                or isinstance(order, bool) or not isinstance(order, int) or order not in (1, 2, 3)):
            raise ValueError('Fixed orders must map valid edge indices to orders 1, 2, or 3')
        initial[index] = (order,)

    local_domains, bounds = [], []
    for members, indices in fragments:
        domains = {k: initial[k] for k in indices}
        possible = propagate(domains, members, n)
        if possible is None:
            raise failure('unsatisfied_valence', 'No valid bond-order assignment; check connectivity and fixed orders.', members)
        bounds.append((sum(min(q for _, q in choices) for choices in possible),
                       sum(max(q for _, q in choices) for choices in possible)))
        local_domains.append(domains)
    minimum, maximum = (sum(v for v, _ in bounds), sum(v for _, v in bounds))
    if not minimum <= charge <= maximum:
        raise failure('incompatible_charge', 'No valid bond-order assignment: requested charge is incompatible with the candidate valence states. Check charge and include all hydrogens.')

    cache = {}
    def fragment_solution(f, target, budget, limit):
        key = (f, target, budget, limit)
        if key in cache:
            return cache[key]
        members, indices = fragments[f]
        def search(domains):
            spend_node()
            possible = propagate(domains, members, budget)
            if possible is None:
                return None
            free = [k for k, orders in domains.items() if len(orders) > 1]
            assignment = charge_paths(possible, budget, not free).get(target)
            if assignment is None:
                return None
            if not free:
                cost, charges, count = assignment
                return cost, ([domains[k][0] for k in indices], charges), count
            k = min(free, key=lambda k: (len(domains[k]), lengths[k], k))
            best = None
            for order in reversed(domains[k]):
                candidate = domains.copy()
                candidate[k] = (order,)
                result = search(candidate)
                if result is None:
                    continue
                if best is None or result[0] < best[0]:
                    best = result
                elif result[0] == best[0]:
                    best = (best[0], best[1], min(2, best[2] + result[2]))
                if best[0] == assignment[0] and best[2] >= limit:
                    break
            return best
        cache[key] = search(local_domains[f].copy())
        return cache[key]

    def combine(budget, limit):
        totals = {0: (0, (), 1)}
        left_min, left_max = minimum, maximum
        for f, (low, high) in enumerate(bounds):
            left_min -= low
            left_max -= high
            following = {}
            for total, (cost, path, count) in totals.items():
                available = budget - cost
                first = max(low, -available, charge - total - left_max)
                last = min(high, available, charge - total - left_min)
                for target in range(first, last + 1):
                    rest = charge - total - target
                    local_budget = min(len(fragments[f][0]), available - abs(rest))
                    if local_budget < abs(target):
                        continue
                    result = fragment_solution(f, target, local_budget, limit)
                    if result is not None:
                        expense, assignment, ways = result
                        extend(following, total + target, cost + expense,
                               path + (assignment,), min(2, count * ways))
            totals = following
            if not totals:
                break
        return totals.get(charge)

    try:
        for budget in range(abs(charge), n + 1, 2):
            result = combine(budget, 1)
            if result is None:
                continue
            cost, path, count = result
            if diagnostics is not None:
                ambiguity = 'multiple' if count > 1 else 'unique'
                if count == 1:
                    try:
                        if combine(cost, 2)[2] > 1:
                            ambiguity = 'multiple'
                    except (BondPerceptionError, RecursionError):
                        # A complete optimal assignment already exists. Never
                        # claim uniqueness when the remaining budget runs out.
                        ambiguity = 'undetermined'
                diagnostics.update(fragments=len(fragments), iterations=iterations,
                                   ambiguity=ambiguity, charge_cost=cost,
                                   seconds=perf_counter() - started)
            orders, charges = [0] * len(edges), [0] * n
            for (members, indices), (bond_orders, atom_charges) in zip(fragments, path):
                for k, order in zip(indices, bond_orders):
                    orders[k] = order
                for i, q in zip(members, atom_charges):
                    charges[i] = q
            return orders, charges
    except RecursionError as exc:
        raise failure('search_depth', 'Bond-order search is too deep; no assignment applied.') from exc
    raise failure('unsatisfied_valence', 'No valid bond-order assignment. Check connectivity, hydrogens and charge; radicals and coordination bonds are not supported.')
