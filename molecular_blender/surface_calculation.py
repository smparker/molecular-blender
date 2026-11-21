# Copyright (C) 2014-2025 Shane Parker, Joshua Szekely
# SPDX-License-Identifier: GPL-3.0-or-later

"""Generate numerical isosurfaces and trajectory visibility without creating scene objects."""
from itertools import chain
import numpy as np

from .isosurfaces import field_isosurface
from .options import get_frame_mapper
from .orbitals import sampled_orbital_overlap
from .surface_geometry import index_field_triangles


def iter_surface_frames(molecule, options, progress_state=None, *, indexed=False):
    """Yield one snapshot's fields and visibility map, including empty trajectory frames."""
    if molecule.volume is None and molecule.orbitals is None:
        return
    fmap = get_frame_mapper(options, getattr(molecule, "frame_numbers", None))

    # Count total tasks for progress tracking
    if progress_state is not None:
        if molecule.volume is not None:
            total = 1 + (len(molecule.volume_trajectory) if molecule.volume_trajectory else 0)
        elif molecule.orbitals is not None:
            total = (1 if options["volume"] != 'orbital' else
                     len([x.strip() for x in options["orbital"].split(',')]))
            total *= 1 + len(molecule.orbitals_trajectory or [])
        else:
            total = 1
        progress_state["total"] = max(total, 1)
        progress_state["completed"] = 0

    def tick():
        if progress_state is not None:
            progress_state["completed"] = progress_state.get("completed", 0) + 1

    if molecule.volume is not None:  # volumetric data was read
        frames = [molecule.volume, *(molecule.volume_trajectory or ())]
        for vol in frames:
            vol.select_dataset(options.get("cube_dataset", 0))
        levels = options["isovalues"]
        if options["cumulative"]:
            reference = frames[0]
            if reference.linear_combination:
                reference = next((vol for vol in frames if np.any(vol.data)), None)
                if reference is None:
                    for _ in frames:
                        tick()
                    return
            levels = reference.isovalue_containing_proportion(levels, options.get("cube_field_type", "auto"))
        for i, vol in enumerate(frames):
            name = "cube-{:5d}".format(i) if i else "cube"
            surfaces = field_isosurface(vol.scalar_field(), levels, options["resolution"], name=name)
            if vol.linear_combination:
                surfaces = [s for s in surfaces if len(s['triangles'])]
            for surface in surfaces:
                surface['series'] = 'cube'
                surface['grid_spacing'] = float(min(vol.axis_norms))
                if indexed:
                    surface['vertices'], surface['faces'] = index_field_triangles(
                        surface.pop('triangles'), surface['grid_spacing'])
            mapping = {fmap(i): [s['name'] for s in surfaces]} if len(frames) > 1 else {}
            tick()
            yield surfaces, mapping

    elif molecule.orbitals is not None:  # orbital data was read
        orbitals = molecule.orbitals
        channel = options.get('orbital_spin', 'all')
        channel = None if channel == 'all' else channel
        if options['volume'] == 'response':
            orbnames = [orbitals.response_field(options.get('response_eigenpair', 0))]
        elif options['volume'] == 'orbital':
            orbnames = [x.strip().lower() for x in options['orbital'].split(',')]
        else:
            orbnames = [options['volume']]
        isovals = options["isovalues"]
        resolution = options["resolution"]

        align_phase = options.get('align_orbital_phase') and molecule.orbitals_trajectory
        previous = {}

        def field_surface(model, name, label, levels, frame):
            ambiguous = False
            if options['volume'] == 'response' and not model.is_density(name):
                raise ValueError(f'Blender frame {frame} is missing DIPL {name}')
            if model.is_density(name):
                if not np.any(model.density_matrix(name)[0]):
                    return [], levels
                scalar = model.scalar_field(name)
            else:
                index = model.resolve_orbital(name, channel)
                scalar = model.scalar_field(index)
                if align_phase:
                    spin = model.spins[index - 1] if model.spins else None
                    if name in previous:
                        last_field, last_spin = previous[name]
                        overlap = sampled_orbital_overlap(last_field, scalar.evaluator) if spin == last_spin else 0.0
                        ambiguous = abs(overlap) < 0.5
                        if overlap <= -0.5:
                            scalar = model.scalar_field(index, flip_phase=True)
                    current[name] = scalar.evaluator, spin
            if levels is None:
                levels = scalar.evaluator.isovalue_containing_proportion(isovals) if options['cumulative'] else isovals
            method = 'adaptive' if model.is_density(name) else 'tiled'
            surfaces = field_isosurface(scalar, levels, resolution, label, method=method)
            if ambiguous:
                for surface in surfaces:
                    surface['ambiguous_phase_frame'] = frame
            if indexed:
                for surface in surfaces:
                    surface['vertices'], surface['faces'] = index_field_triangles(surface.pop('triangles'), resolution)
            return surfaces, levels

        initial_levels = {}
        for i, snapshot in enumerate(chain((orbitals,), molecule.orbitals_trajectory or ())):
            surfaces = []
            current = {}
            for name in orbnames:
                label = f"{name}-{i:4d}" if i else name
                fields, levels = field_surface(snapshot, name, label, initial_levels.get(name), fmap(i))
                if i == 0:
                    initial_levels[name] = levels
                if i:
                    for field in fields:
                        field['series'] = name
                surfaces.extend(fields)
                tick()
            previous = current
            mapping = {fmap(i): [s['name'] for s in surfaces]} if molecule.orbitals_trajectory is not None else {}
            yield surfaces, mapping


def compute_surface_vertex_sets(molecule, options, progress_state=None, *, indexed=False):
    """Collect numerical results for callers that request the full calculation."""
    definitions, mapping = [], {}
    for surfaces, frames in iter_surface_frames(molecule, options, progress_state, indexed=indexed):
        definitions.extend(surfaces)
        mapping.update(frames)
    if molecule.orbitals is not None and molecule.orbitals_trajectory and options['volume'] == 'orbital':
        # The collecting API historically returns all snapshots of one orbital together.
        order = {name.strip().lower(): i for i, name in enumerate(options['orbital'].split(','))}
        definitions.sort(key=lambda field: order[field.get('series', field['name'])])
    return definitions, mapping
