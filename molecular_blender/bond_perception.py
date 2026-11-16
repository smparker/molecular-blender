"""Valence-constrained bond orders for an existing distance-based graph.

Inspired by RDKit DetermineBonds / xyz2mol's separation of connectivity, allowed
valences, unsaturation and charge validation. This is an independent constrained
search, not a port of their maximum-matching implementation:
https://github.com/rdkit/rdkit/blob/master/Code/GraphMol/DetermineBonds/DetermineBonds.cpp
https://github.com/jensengroup/xyz2mol

All atoms (including hydrogens) must be supplied. No radical, metal coordination,
stereochemistry or aromaticity model is inferred. Rings get integer Kekule orders.
"""
import math


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

    def __init__(self, message, *, reason=None, diagnostics=None):
        super().__init__(message)
        self.reason = reason
        self.diagnostics = {} if diagnostics is None else diagnostics


def perceive_bond_orders(symbols, edges, charge=0, bond_lengths=None,
                         max_iterations=100000, atom_ids=None, fixed_orders=None):
    """Return (orders in input edge order, formal charges in atom order).

    Every edge starts with order one and can rise to at most three. Propagation
    removes orders incompatible with the endpoint valences. Backtracking tries
    the most constrained edges first, using shorter bonds to break ties. Search
    uses increasing budgets for total absolute formal charge, so a neutral
    assignment wins over unnecessary charge separation. Equal-cost resonance
    forms are heuristic choices, not unique chemical determinations.

    Neither the inputs nor a molecule are mutated, including on failure. The
    node budget covers all charge budgets; exhaustion never returns a partial
    assignment. bond_lengths may be normalized by covalent radii by the caller.
    atom_ids optionally supplies stable caller IDs for diagnostics; edges still
    use positional indices. fixed_orders maps input edge indices to locked
    integer orders (1..3); incompatible locks fail without a partial assignment.
    Errors distinguish supported-chemistry limitations,
    impossible connectivity, incompatible charge, and exhausted search.
    """
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

    def propagate(domains):
        """Prune by endpoint valence capacity until no further order changes."""
        while True:
            changed = False
            possible = []
            for i, bonds in enumerate(incident):
                low = sum(domains[k][0] for k in bonds)
                high = sum(domains[k][-1] for k in bonds)
                choices = [(v, q) for v, q in states[i] if low <= v <= high]
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

    def charges_for(possible):
        # Dynamic programming gives both a lower bound on charge separation and
        # a concrete charge assignment when the bond domains are singletons.
        totals = {0: (0, ())}
        for choices in possible:
            next_totals = {}
            for total, (cost, charges) in totals.items():
                for q in sorted({q for _, q in choices}, key=lambda q: (abs(q), q)):
                    target = total + q
                    candidate = (cost + abs(q), charges + (q,))
                    if target not in next_totals or candidate[0] < next_totals[target][0]:
                        next_totals[target] = candidate
            totals = next_totals
        return totals.get(charge)

    def search(domains, budget, failed):
        nonlocal iterations
        if iterations >= max_iterations:
            raise failure('search_limit', 'Bond-order search limit reached; no assignment applied.')
        iterations += 1
        possible = propagate(domains)
        if possible is None:
            return None
        key = tuple(domains)
        if key in failed:
            return None
        charges = charges_for(possible)
        if charges is None or charges[0] > budget:
            failed.add(key)
            return None
        free = [k for k, orders in enumerate(domains) if len(orders) > 1]
        if not free:
            return [orders[0] for orders in domains], list(charges[1])
        # Raising orders on short unsaturated edges is the educated first guess;
        # every alternative is still tried if that guess blocks a valid solution.
        k = min(free, key=lambda k: (len(domains[k]), lengths[k], k))
        for order in reversed(domains[k]):
            candidate = domains.copy()
            candidate[k] = (order,)
            result = search(candidate, budget, failed)
            if result is not None:
                return result
        failed.add(key)
        return None

    initial = [(1, 2, 3)] * len(edges)
    for index, order in (fixed_orders or {}).items():
        if (isinstance(index, bool) or not isinstance(index, int) or not 0 <= index < len(edges)
                or isinstance(order, bool) or not isinstance(order, int) or order not in (1, 2, 3)):
            raise ValueError('Fixed orders must map valid edge indices to orders 1, 2, or 3')
        initial[index] = (order,)
    possible = propagate(initial.copy())
    if possible is not None and charges_for(possible) is None:
        raise failure('incompatible_charge',
                      'No valid bond-order assignment: requested charge is incompatible with the candidate valence states. Check charge and include all hydrogens.')
    try:
        # All supported formal charges are -1, 0 or +1; cost and net charge have
        # the same parity. No bigger budget than n can yield a new assignment.
        for budget in range(abs(charge), n + 1, 2):
            result = search(initial.copy(), budget, set())
            if result is not None:
                return result
    except RecursionError as exc:
        raise failure('search_depth', 'Bond-order search is too deep; no assignment applied.') from exc
    raise failure('unsatisfied_valence',
        'No valid bond-order assignment for molecular charge {}. Check connectivity, '
        'include all hydrogens, and check charge; radicals and coordination bonds '
        'are not supported.'.format(charge))
