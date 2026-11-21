# Copyright (C) 2014-2025 Shane Parker, Joshua Szekely
# SPDX-License-Identifier: GPL-3.0-or-later

"""Material bindings for one molecule and opt-in reuse of matching generated shaders."""
import hashlib
import json

import bpy

from .periodictable import elements


def _appearance(mat):
    """Fingerprint the generated shader's inputs and links before reusing it."""
    if not mat.use_nodes or mat.animation_data or mat.node_tree.animation_data:
        return None
    nodes = list(mat.node_tree.nodes)
    indices = {node: i for i, node in enumerate(nodes)}
    settings = []
    for node in nodes:
        if node.type not in {'OUTPUT_MATERIAL', 'BSDF_PRINCIPLED', 'BSDF_DIFFUSE', 'MIX_SHADER',
                             'LAYER_WEIGHT', 'MATH', 'BSDF_TRANSPARENT', 'EMISSION'}:
            return None
        inputs = []
        for socket in node.inputs:
            value = getattr(socket, 'default_value', None)
            if value is not None and not isinstance(value, (str, int, float)):
                value = tuple(value)
            inputs.append(value)
        settings.append((node.bl_idname, node.mute, getattr(node, 'operation', None),
                         getattr(node, 'distribution', None), getattr(node, 'use_clamp', None),
                         getattr(node, 'subsurface_method', None), inputs))
    links = [(indices[link.from_node], link.from_socket.identifier,
              indices[link.to_node], link.to_socket.identifier) for link in mat.node_tree.links]
    state = (tuple(mat.diffuse_color), mat.blend_method, settings, sorted(links))
    return hashlib.sha256(repr(state).encode()).hexdigest()


def material_for(styler, options, role, label, name, create):
    key = json.dumps([1, type(styler).__name__, role, label])
    if options['recycle_materials']:
        for material in bpy.data.materials:
            if (material.library is None and material.get('mb_material_key') == key
                    and not material.get('mb_style_material')
                    and material.get('mb_material_signature') == _appearance(material)):
                return material
    material = create(name)
    material['mb_generated'] = True
    if options['recycle_materials']:
        signature = _appearance(material)
        if signature is not None:
            material['mb_material_key'] = key
            material['mb_material_signature'] = signature
    return material


def element_material(styler, options, symbol):
    return material_for(styler, options, 'atom', symbol, symbol + '_mat',
                        lambda name: styler.atom_material(name, elements[symbol]))


def bond_material(styler, options, style):
    return material_for(styler, options, 'bond', style, style + '_bond_mat',
                        lambda name: styler.bond_material(name, style))


def make_materials(molecule, styler, options):
    make_atom_materials(molecule, styler, options)
    if molecule.bonds:
        make_bond_materials(molecule, styler, options)
    if molecule.rings:
        make_ring_materials(molecule, styler, options)


def make_atom_materials(molecule, styler, options):
    for symbol in sorted({atom.element.symbol for atom in molecule.atoms}):
        molecule.materials[symbol] = element_material(styler, options, symbol).name
    if options['charges'] != 'none':
        positive, negative = styler.charge_material('pluscharge_mat', 'negcharge_mat', None)
        molecule.pluscharge_mat, molecule.negcharge_mat = positive.name, negative.name


def make_bond_materials(molecule, styler, options):
    for style in sorted({bond.style for bond in molecule.bonds}):
        molecule.bond_materials[style] = bond_material(styler, options, style).name


def make_ring_materials(molecule, styler, options):
    mat = material_for(styler, options, 'ring', '', 'ring_mat', styler.ring_material)
    molecule.materials['ring'] = mat.name


def make_iso_materials(molecule, styler, vertex_sets, options):
    groups = {}
    for surface in vertex_sets:
        phase = 'plus' if surface['isovalue'] >= 0 else 'negative'
        groups.setdefault((surface['name'], phase), []).append(surface)

    materials = {}
    for (field, phase), surfaces in groups.items():
        ordered = sorted(surfaces, key=lambda s: abs(s['isovalue']), reverse=True)
        for index, surface in enumerate(ordered):
            solid = index == 0 and len(ordered) > 1
            layer = 'inner' if index == 0 else 'outer'
            name = f'{field}_orbital_mat_{phase}_{layer}'
            material = material_for(styler, options, 'solid' if solid else 'outer', [field, phase], name,
                                    styler.isosurface_material if solid else styler.outer_isosurface_material)
            material['mb_surface_material'] = True
            materials[material.name] = surface['material'] = material
    return materials
