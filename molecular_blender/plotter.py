# Copyright (C) 2014-2025 Shane Parker, Joshua Szekely
# SPDX-License-Identifier: GPL-3.0-or-later

"""Blender import driver and historical plotting entry points."""
from pathlib import Path

from .bond_perception import BondPerceptionError
from .diagnostics import InputError
from .geometry_import import geometry_import
from .importers import molecule_from_file
from .molecule import Molecule
from .options import default_options, get_frame_mapper, process_options, validated_options
from .scene_resources import mark_import_resources
from .scene_animation import AnimateMolecule
from .scene_geometry import (PlotMolecule, PlotWireFrame, bind_bond_endpoints, plot_gradients,
                             plot_rings, set_bond_order, _multiple_bond_profile)
from .scene_materials import (make_atom_materials, make_bond_materials, make_iso_materials,
                              make_materials, make_ring_materials)
from .scene_meshes import (SurfaceMeshStage, create_geometry, create_mesh, create_surface_mesh_data,
                           create_surface_meshes, draw_surfaces)
from .stylers import get_styler
from .surface_calculation import compute_surface_vertex_sets, iter_surface_frames
from .surface_geometry import field_geometry, field_weld_distance, index_field_triangles


def plot_prep(molecule, options):
    """Preprocess options"""
    process_options(options, cube=molecule.volume is not None)
    if molecule.volume is None and molecule.orbitals is None:
        options["isosurfaces"] = False
    if options["charges"] != "none":
        molecule.chgfac = options["charge_factor"]
        molecule.chgoff = options["charge_offset"]
    return options


def BlendMoleculeSetup(context, filename, **options):
    """Load molecule and create atoms/bonds. Returns (molecule, styler, options).

    Does not draw isosurfaces - call draw_surfaces separately if needed.
    """
    options = validated_options(options)
    name = Path(filename).stem
    input_data = molecule_from_file(filename, options)
    try:
        molecule = Molecule.from_dict(name, input_data)
        molecule.determine_bonding(options)
    except BondPerceptionError as exc:
        raise BondPerceptionError(str(exc), reason=exc.reason, diagnostics=exc.diagnostics, source=filename) from exc
    except ValueError as exc:
        raise InputError(filename, "Molecular model", str(exc)) from exc
    molecule.surface_input = input_data
    styler = get_styler(options)
    molecule.settings = dict(options)

    with geometry_import(context) as created:
        make_materials(molecule, styler, options)
        modes = input_data.get('vibrations')
        if options['object_type'] == 'wireframe' or input_data.get('frame_numbers'):
            modes = None
        if modes and options['hook_atoms'] == 'auto':
            options['hook_atoms'] = 'on'
        plot_prep(molecule, options)

        if options["object_type"] == "wireframe":
            PlotWireFrame(context, molecule, options)
            if options["plot_type"] == "animate":
                print("animation of wireframes not currently supported!")
        else:
            PlotMolecule(context, molecule, options)
            if options["plot_type"] == "animate":
                AnimateMolecule(context, molecule, options)

        if molecule.rings:
            plot_rings(context, molecule, options)

        from .scene_state import persist_molecule
        persist_molecule(molecule, filename, molecule.settings)
        if modes:
            from .vibrations import store_modes
            store_modes(molecule.scene_owner, modes)
        if molecule.settings.get('atom_scale', 1.0) != 1.0:
            from .scene_edit import restyle_molecule
            restyle_molecule(molecule.scene_owner, atom_scale=molecule.settings['atom_scale'],
                             bond_thickness=molecule.settings['bond_thickness'],
                             colors=molecule.settings['colors'],
                             surface_opacity=molecule.settings.get('surface_opacity', 1.0))
        mark_import_resources(molecule.scene_owner, created)
    return molecule, styler, options


def BlendMolecule(context, filename, **options):
    """basic driver that calls the appropriate plot functions"""
    molecule, styler, options = BlendMoleculeSetup(context, filename, **options)
    if options["isosurfaces"]:
        draw_surfaces(molecule, styler, context, options)
