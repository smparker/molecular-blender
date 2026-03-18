"""Tests for molecular_blender.aromatics"""

import unittest
import numpy as np

from molecular_blender.aromatics import (
    create_graph, depth_first, find_specific_cycles,
    find_unique_cycles, ordered_cycle, is_planar,
)


class _Vec:
    """Minimal vector class mimicking mathutils.Vector for testing is_planar."""

    def __init__(self, xyz):
        self._v = np.array(xyz, dtype=float)

    def __sub__(self, other):
        return _Vec(self._v - other._v)

    def normalize(self):
        n = np.linalg.norm(self._v)
        if n > 0:
            self._v /= n

    def dot(self, other):
        return float(np.dot(self._v, other._v))

    def cross(self, other):
        return _Vec(np.cross(self._v, other._v))


class _MockAtom:
    def __init__(self, x, y, z):
        self.position = _Vec([x, y, z])


class TestCreateGraph(unittest.TestCase):
    def test_simple_chain(self):
        graph = create_graph([0, 1, 2], [(0, 1), (1, 2)])
        assert 1 in graph[0]
        assert 0 in graph[1]
        assert 2 in graph[1]
        assert 1 in graph[2]
        assert 2 not in graph[0]

    def test_triangle(self):
        graph = create_graph([0, 1, 2], [(0, 1), (1, 2), (2, 0)])
        for i in range(3):
            assert len(graph[i]) == 2

    def test_empty(self):
        graph = create_graph([0, 1], [])
        assert len(graph[0]) == 0
        assert len(graph[1]) == 0


class TestDepthFirst(unittest.TestCase):
    def test_triangle_start_to_neighbor(self):
        """depth_first finds paths from start to goal (different vertices)."""
        graph = create_graph([0, 1, 2], [(0, 1), (1, 2), (2, 0)])
        paths = list(depth_first(graph, 0, 1))
        assert len(paths) > 0
        for path in paths:
            assert path[0] == 0
            assert path[-1] == 1

    def test_no_path(self):
        """No path between disconnected vertices."""
        graph = create_graph([0, 1, 2], [(0, 1)])
        paths = list(depth_first(graph, 0, 2))
        assert len(paths) == 0

    def test_hexagon_paths(self):
        edges = [(i, (i + 1) % 6) for i in range(6)]
        graph = create_graph(list(range(6)), edges)
        paths = list(depth_first(graph, 0, 3))
        assert len(paths) > 0


class TestFindCycles(unittest.TestCase):
    def test_triangle(self):
        graph = create_graph([0, 1, 2], [(0, 1), (1, 2), (2, 0)])
        cycles = find_unique_cycles(graph)
        assert len(cycles) >= 1
        assert any(len(c) == 3 for c in cycles)

    def test_hexagon(self):
        edges = [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 0)]
        graph = create_graph(list(range(6)), edges)
        cycles = find_unique_cycles(graph)
        assert len(cycles) >= 1
        assert any(len(c) == 6 for c in cycles)

    def test_no_cycles(self):
        graph = create_graph([0, 1, 2, 3], [(0, 1), (1, 2), (2, 3)])
        cycles = find_unique_cycles(graph)
        assert len(cycles) == 0

    def test_two_fused_rings(self):
        # naphthalene-like: two fused 6-rings sharing edge 2-3
        edges = [
            (0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 0),
            (2, 6), (6, 7), (7, 8), (8, 3),
        ]
        graph = create_graph(list(range(9)), edges)
        cycles = find_unique_cycles(graph)
        # should find at least the two 6-membered and 5-membered rings
        assert len(cycles) >= 2


class TestOrderedCycle(unittest.TestCase):
    def test_triangle(self):
        edges = [(0, 1), (1, 2), (2, 0)]
        graph = create_graph([0, 1, 2], edges)
        cycle = {0, 1, 2}
        result = ordered_cycle(cycle, graph)
        assert set(result) == cycle
        assert len(result) == 3

    def test_hexagon(self):
        edges = [(i, (i + 1) % 6) for i in range(6)]
        graph = create_graph(list(range(6)), edges)
        cycle = set(range(6))
        result = ordered_cycle(cycle, graph)
        assert set(result) == cycle
        assert len(result) == 6


class TestFindSpecificCycles(unittest.TestCase):
    def test_triangle_from_vertex(self):
        edges = [(0, 1), (1, 2), (2, 0)]
        graph = create_graph([0, 1, 2], edges)
        cycles, fixed = find_specific_cycles(graph, 0)
        assert len(cycles) >= 1
        assert all(0 in c for c in cycles)


class TestIsPlanar(unittest.TestCase):
    def test_triangle_always_planar(self):
        atoms = [_MockAtom(0, 0, 0), _MockAtom(1, 0, 0), _MockAtom(0, 1, 0)]
        assert is_planar(atoms, {0, 1, 2}) is True

    def test_square_planar(self):
        atoms = [
            _MockAtom(0, 0, 0), _MockAtom(1, 0, 0),
            _MockAtom(1, 1, 0), _MockAtom(0, 1, 0),
        ]
        assert is_planar(atoms, [0, 1, 2, 3]) is True

    def test_hexagon_planar(self):
        import math
        atoms = [_MockAtom(math.cos(i * math.pi / 3), math.sin(i * math.pi / 3), 0)
                 for i in range(6)]
        assert is_planar(atoms, list(range(6))) is True

    def test_nonplanar(self):
        """Tetrahedron vertices are not coplanar."""
        atoms = [
            _MockAtom(1, 1, 1), _MockAtom(1, -1, -1),
            _MockAtom(-1, 1, -1), _MockAtom(-1, -1, 1),
        ]
        assert is_planar(atoms, [0, 1, 2, 3]) is False
