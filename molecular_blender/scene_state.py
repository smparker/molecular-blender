"""Versioned scene records and molecule ownership independent of object names.

Atom/bond/surface IDs are local to a molecule. Parent ancestry is authoritative:
Blender duplicates custom properties verbatim, but remaps copied parent links.
A duplication handler allocates a fresh root UUID without rebinding children to
an unrelated original molecule. Records contain no pickles or executable data.
"""
import json
import math
import os
import uuid

import bpy
from bpy.app.handlers import persistent

from .options import validated_options

RECORD = 'mb_molecule'
ROLE = 'mb_role'
ID = 'mb_id'
MOLECULE_ID = 'mb_molecule_id'
VERSION = 1
_normalizing = False


def new_id():
    return uuid.uuid4().hex


def resolve_owner(obj):
    """Find the nearest molecule root, including for untagged child decorations."""
    while obj is not None:
        if RECORD in obj:
            return obj
        obj = obj.parent
    return None


def owned_objects(owner, role=None):
    return [obj for obj in owner.children_recursive
            if resolve_owner(obj) == owner and (role is None or obj.get(ROLE) == role)]


def read_record(owner):
    record = json.loads(owner[RECORD])
    if record.get('version') != VERSION:
        raise ValueError('Unsupported molecule settings version; reimport with this addon version')
    record['settings'] = validated_options(record['settings'])
    return record


def write_record(owner, record):
    if record.get('version') != VERSION:
        raise ValueError('Unsupported molecule settings version')
    record = dict(record, settings=validated_options(record['settings']))
    owner[RECORD] = json.dumps(record, allow_nan=False, separators=(',', ':'))
    owner[MOLECULE_ID] = record['id']
    owner[ROLE] = 'molecule'


def tag(obj, owner, role, identifier):
    obj[ROLE] = role
    obj[ID] = identifier
    obj[MOLECULE_ID] = owner[MOLECULE_ID]


@persistent
def normalize_duplicate_ids(*_args):
    global _normalizing
    if _normalizing:
        return
    _normalizing = True
    try:
        seen = set()
        for owner in bpy.data.objects:
            if RECORD not in owner:
                continue
            identifier = owner.get(MOLECULE_ID)
            if identifier in seen:
                try:
                    record = read_record(owner)
                except (ValueError, KeyError):
                    continue  # Never rewrite an unknown future schema.
                record['id'] = new_id()
                write_record(owner, record)
                for obj in owned_objects(owner):
                    obj[MOLECULE_ID] = record['id']
                identifier = record['id']
            seen.add(identifier)
    finally:
        _normalizing = False


def orbital_metadata(model):
    if model is None:
        return []
    return [dict(index=i+1, energy=float(e) if e is not None and math.isfinite(e) else None,
                 symmetry=sym, occupation=float(occ), spin=model.spins[i] if model.spins is not None else None)
            for i, (e, sym, occ) in enumerate(zip(model.energies, model.symmetries, model.occupations))]


def persist_molecule(molecule, source, settings):
    owner = molecule.scene_owner
    atoms = []
    atom_ids = {}
    for atom in molecule.atoms:
        identifier = new_id()
        atom_ids[atom.index] = identifier
        atoms.append(dict(id=identifier, index=int(atom.index), element=atom.element.symbol,
                          position=list(atom.position), charge=float(atom.charge),
                          formal_charge=int(atom.formal_charge), hidden=bool(atom.hidden)))
    provenance = ('explicit' if molecule.explicit_bonds is not None else
                  'perceived' if settings['perceive_bond_orders'] else 'distance')
    bonds = [dict(id=new_id(), atoms=[atom_ids[b.iatom.index], atom_ids[b.jatom.index]],
                  order=int(b.order), provenance=provenance, locked=False,
                  direction=list(b.direction) if b.direction is not None else None)
             for b in molecule.bonds]
    record = dict(version=VERSION, id=new_id(), source=os.path.abspath(source), units='Angstrom',
                  settings=settings, atoms=atoms, bonds=bonds, surfaces=[], orbitals=orbital_metadata(molecule.orbitals),
                  fields=(dict(datasets=list(molecule.volume.labels), field_type=molecule.volume.field_type)
                          if molecule.volume is not None else None))
    if molecule.volume is not None:
        for key in ("units", "linear_combination"):
            value = getattr(molecule.volume, key)
            if value is not None:
                record["fields"][key] = value
    write_record(owner, record)
    for atom, entry in zip(molecule.atoms, atoms):
        obj = bpy.data.objects.get(atom.name)
        if obj is not None and resolve_owner(obj) == owner:
            tag(obj, owner, 'atom', entry['id'])
            obj['mb_element'] = atom.element.symbol
            obj['mb_base_scale'] = list(obj.scale)
    for bond, entry in zip(molecule.bonds, bonds):
        for name in bond.name if isinstance(bond.name, list) else []:
            obj = bpy.data.objects.get(name)
            if obj is not None:
                tag(obj, owner, 'bond', entry['id'])
                obj['mb_endpoints'] = entry['atoms']
                obj['mb_order'] = entry['order'] if settings['draw_bond_orders'] else 1
                obj['mb_radius'] = obj.data.bevel_depth
                obj['mb_nominal_radius'] = settings['bond_thickness']
                if settings['plot_style'] == 'sticks' and obj.active_material is not None:
                    for symbol, material_name in molecule.materials.items():
                        if obj.active_material.name == material_name:
                            obj['mb_element'] = symbol
                            break
    return owner


def persist_surfaces(owner, definitions, options):
    """Each calculation of a field has one switch for all signs, levels and frames."""
    record = read_record(owner)
    groups = {}
    for definition in definitions:
        field = definition.get('series', definition['name'])
        if field not in groups:
            groups[field] = new_id()
        identifier = new_id()
        tag(definition['mesh'], owner, 'surface', identifier)
        record['surfaces'].append(dict(id=identifier, field=field, group=groups[field],
                                       orbital_spin=options.get('orbital_spin', 'all'),
                                       isovalue=float(definition['isovalue'])))
    write_record(owner, record)


def register():
    # Blender restricts scene access during add-on registration.
    for handlers in (bpy.app.handlers.load_post, bpy.app.handlers.depsgraph_update_post):
        if normalize_duplicate_ids not in handlers:
            handlers.append(normalize_duplicate_ids)


def unregister():
    for handlers in (bpy.app.handlers.load_post, bpy.app.handlers.depsgraph_update_post):
        if normalize_duplicate_ids in handlers:
            handlers.remove(normalize_duplicate_ids)
