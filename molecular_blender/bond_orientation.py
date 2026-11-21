"""Orient an extruded bond profile in Blender's Z_UP curve frame."""
from collections import deque
import math

from mathutils import Vector


def profile_axes(tangent):
    xaxis = tangent.cross(Vector((0, 0, 1)))
    if xaxis.length < 1e-8:
        xaxis = Vector((1 if tangent.z > 0 else -1, 0, 0))
    xaxis.normalize()
    return xaxis, xaxis.cross(tangent)


def strand_tilt(tangent, direction):
    if direction is None:
        return 0.0
    separation = Vector(direction)
    separation -= tangent * separation.dot(tangent)
    if separation.length < 1e-8:
        raise ValueError("Bond direction must have a component perpendicular to the bond")
    xaxis, yaxis = profile_axes(tangent)
    return math.atan2(-separation.dot(yaxis), separation.dot(xaxis))


def _bond_ring(first, second, neighbors):
    """Find the smallest ring containing this edge, independently of aromatic styling."""
    parents = {first: None}
    queue = deque([first])
    while queue:
        atom = queue.popleft()
        for neighbor in neighbors[atom]:
            if neighbor in parents or (atom is first and neighbor is second):
                continue
            parents[neighbor] = atom
            if neighbor is second:
                ring = [second]
                while ring[-1] is not first:
                    ring.append(parents[ring[-1]])
                return ring
            queue.append(neighbor)
    return []


def trajectory_tilts(bond, neighbors):
    """Yield sampled profile angles. A pair of strands is unchanged by a half turn."""
    first, second = bond.iatom, bond.jatom
    ring = _bond_ring(first, second, neighbors)
    candidates = [(end, atom) for end, other in ((first, second), (second, first))
                  for atom in neighbors[end] if atom is not other]
    reference = None
    last_axis = last_direction = None
    angle = 0.0
    for frame, (a, b) in enumerate(zip(first.trajectory, second.trajectory)):
        axis = b.position - a.position
        if axis.length < 1e-8:
            yield angle
            continue
        axis.normalize()
        direction = None
        if ring:
            # Sum triangle area vectors relative to one corner for translation stability.
            offsets = [atom.trajectory[frame].position - a.position for atom in ring[:-1]]
            normal = sum((p.cross(q) for p, q in zip(offsets, offsets[1:])), Vector())
            scale = max(p.length_squared for p in offsets)
            in_plane = normal.cross(axis)
            if in_plane.length > 1e-4 * scale:
                direction = in_plane.normalized()
        elif candidates:
            # Keep the same reference atom when substituents have similar angles.
            best = 0.0
            for end, atom in ([reference] if reference else candidates):
                delta = atom.trajectory[frame].position - end.trajectory[frame].position
                projected = delta - axis * delta.dot(axis)
                if projected.length > max(best, 1e-4 * delta.length, 1e-8):
                    reference = (end, atom)
                    direction = projected.normalized()
                    best = projected.length
        if direction is None:
            # Transport the last orientation through linear or isolated geometries.
            direction = (last_axis.rotation_difference(axis) @ last_direction
                         if last_axis is not None else profile_axes(axis)[0])
        new_angle = strand_tilt(axis, direction)
        angle = new_angle if last_axis is None else angle + math.remainder(new_angle - angle, math.pi)
        last_direction, last_axis = direction, axis
        yield angle
