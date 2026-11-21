# Copyright (C) 2014-2025 Shane Parker, Joshua Szekely
# SPDX-License-Identifier: GPL-3.0-or-later

"""Keyframe imported geometry and multiple-bond strand orientation."""
import bpy

from .bond_orientation import trajectory_tilts
from .options import get_frame_mapper


def AnimateMolecule(context, molecule, options):
    """Keyframes position of atoms and opacity of bonds"""
    frame_map = get_frame_mapper(options, molecule.frame_numbers)
    interpolation = {'hold': 'CONSTANT', 'linear': 'LINEAR'}.get(options.get('trajectory_interpolation'))

    # check to make sure trajectory information is stored for at least the first atom
    if not molecule.atoms[0].trajectory:
        raise Exception(
            "Trajectory must be read before calling AnimateMolecule")
    if all([atom.name == '' for atom in molecule.atoms]):
        raise Exception(
            "No names found for any atom in Molecule. Was PlotMolecule or PlotAtoms called?")
    for item in context.selected_objects:
        item.select_set(False)
    for atom in filter(lambda x: x.name != '', molecule.atoms):
        atom_obj = bpy.data.objects.get(atom.name)
        if atom_obj is None:
            raise Exception(
                "Atom " + atom.name + " not found. Was PlotMolecule or PlotAtoms called properly?")

        if atom.hidden:
            continue

        pluscharge = bpy.data.objects[atom.plus_charge] if (
            options["charges"] != "none") else ""
        negcharge = bpy.data.objects[atom.neg_charge] if (
            options["charges"] != "none") else ""

        for (iframe, isnap) in enumerate(atom.trajectory):
            atom_obj.location = isnap.position
            atom_obj.keyframe_insert(
                data_path='location', frame=frame_map(iframe))

            if options["charges"] != "none":
                if options["charges"] == "scale":
                    pc, nc = max(0.0, isnap.charge), -min(0.0, isnap.charge)
                    pluscharge.scale = molecule.scale(pc)
                    negcharge.scale = molecule.scale(nc)
                    for x in [pluscharge, negcharge]:
                        x.keyframe_insert(data_path='scale',
                                          frame=frame_map(iframe))

        if interpolation:
            animated = [(atom_obj, 'location')]
            if options['charges'] == 'scale':
                animated.extend(((pluscharge, 'scale'), (negcharge, 'scale')))
            for obj, path in animated:
                for curve in obj.animation_data.action.fcurves:
                    if curve.data_path == path:
                        for point in curve.keyframe_points:
                            point.interpolation = interpolation

    if options["draw_bond_orders"] and options["hook_atoms"]:
        _animate_strand_tilts(molecule, frame_map, interpolation)

    if options["animate_bonds"] == "dynamic" and molecule.explicit_bonds is None:
        bondmask = molecule.bond_mask(options)
        for bond in molecule.bonds:
            i, j = bond.iatom.index, bond.jatom.index
            for bond_obj in [ bpy.data.objects[b] for b in bond.name ]:
                for (iframe, mask) in enumerate(bondmask[(i, j)]):
                    bond_obj.hide_viewport = not mask
                    bond_obj.keyframe_insert(
                        data_path="hide_viewport", frame=frame_map(iframe))
                    bond_obj.hide_render = not mask
                    bond_obj.keyframe_insert(
                        data_path="hide_render", frame=frame_map(iframe))
    return


def _animate_strand_tilts(molecule, frame_map, interpolation):
    bonds = [bond for bond in molecule.bonds if bond.order > 1 and bond.name]
    if not bonds:
        return
    neighbors = {atom: [] for atom in molecule.atoms}
    for bond in molecule.bonds:
        neighbors[bond.iatom].append(bond.jatom)
        neighbors[bond.jatom].append(bond.iatom)
    for bond in bonds:
        curves = [bpy.data.objects[name].data for name in bond.name]
        for i, tilt in enumerate(trajectory_tilts(bond, neighbors)):
            for curve in curves:
                for point in curve.splines[0].bezier_points:
                    point.tilt = tilt
                    point.keyframe_insert(data_path='tilt', frame=frame_map(i))
        for curve in curves:
            curve['mb_trajectory_tilt'] = True
            if interpolation:
                for fcurve in curve.animation_data.action.fcurves:
                    for key in fcurve.keyframe_points:
                        key.interpolation = interpolation
