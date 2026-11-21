"""Time charge-conserving assignments for fused rings and disconnected fragments.

Run with python3 benchmarks/bond_perception.py. Pass an older bond_perception.py
as an optional argument to time its public API on the same graphs.
"""
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path
from statistics import median
import sys
from time import perf_counter

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from molecular_blender.bond_perception import perceive_bond_orders


def acene(rings):
    # Two rails with rungs at alternating positions form linearly fused hexagons.
    width = 2*rings + 1
    size = 2*width
    edges = [(side+i, side+i+1) for side in (0, width) for i in range(width-1)]
    edges.extend((i, i+width) for i in range(0, width, 2))
    symbols = ['c']*size
    for i in range(size):
        if sum(i in edge for edge in edges) == 2:
            edges.append((i, len(symbols)))
            symbols.append('h')
    return symbols, edges


def copies(graph, count):
    atoms, bonds = graph
    return (atoms*count, [(i+k*len(atoms), j+k*len(atoms)) for k in range(count) for i, j in bonds])


if __name__ == '__main__':
    before = len(sys.argv) > 1
    if before:
        spec = spec_from_file_location('baseline_perception', sys.argv[1])
        baseline = module_from_spec(spec)
        spec.loader.exec_module(baseline)
        perceive_bond_orders = baseline.perceive_bond_orders
    for label, (atoms, edges) in [('benzene', acene(1)), ('six fused rings', acene(6)),
                                   ('eight separate benzenes', copies(acene(1), 8))]:
        times = []
        stats = {}
        for trial in range(4):
            now = perf_counter()
            orders, charges = perceive_bond_orders(atoms, edges, **({} if before else {'diagnostics':stats}))
            times.append(perf_counter()-now)
            assert sum(charges) == 0 and orders.count(2) == atoms.count('c')//2
        print(label, 'ms', round(1000*median(times[1:]), 3), stats, flush=True)
