"""Inspect the molecule owning the current Blender selection."""
import os

import bpy
from bpy.props import BoolProperty, IntProperty

from .scene_state import ID, ROLE, read_record, resolve_owner, owned_objects

REPRESENTATIONS = {'vdw': 'Van der Waals', 'bs': 'Ball and stick',
                   'fixedbs': 'Fixed ball and stick', 'sticks': 'Sticks'}
ORBITALS_PER_PAGE = 5


def molecule_summary(obj):
    """Return saved metadata without loading the source or evaluating fields."""
    owner = resolve_owner(obj)
    if owner is None:
        return None
    record = read_record(owner)
    settings = record['settings']
    atom = None
    if obj.get(ROLE) == 'atom':
        atom_id = obj.get(ID)
        atom = next((a for a in record['atoms'] if a['id'] == atom_id), None)
    return dict(owner=owner, source=record['source'], units=record['units'],
                selected_atom=atom,
                atoms=len(record['atoms']), visible_atoms=len(owned_objects(owner, 'atom')),
                bonds=len(record['bonds']), surfaces=len(owned_objects(owner, 'surface')),
                representation=REPRESENTATIONS[settings['plot_style']],
                object_type=settings['object_type'], charge=settings['molecular_charge'],
                orbitals=record['orbitals'], fields=record['fields'])


def _page_count(orbitals):
    return max(1, (len(orbitals) + ORBITALS_PER_PAGE - 1) // ORBITALS_PER_PAGE)


def _get_orbital_page(owner):
    return owner.get('mb_orbital_page', 1)


def _set_orbital_page(owner, page):
    orbitals = read_record(owner)['orbitals']
    owner['mb_orbital_page'] = max(1, min(page, _page_count(orbitals)))


def _draw_orbitals(box, owner, orbitals):
    box.prop(owner, 'mb_show_orbitals', text='Orbitals ({})'.format(len(orbitals)),
             icon='TRIA_DOWN' if owner.mb_show_orbitals else 'TRIA_RIGHT', emboss=False)
    if not owner.mb_show_orbitals:
        return
    pages = _page_count(orbitals)
    if pages > 1:
        box.prop(owner, 'mb_orbital_page', text='Page (1–{})'.format(pages))
    start = (min(owner.mb_orbital_page, pages) - 1) * ORBITALS_PER_PAGE
    for orbital in orbitals[start:start + ORBITALS_PER_PAGE]:
        entry = box.column(align=True)
        entry.label(text='Orbital {}'.format(orbital['index']))
        row = entry.row()
        energy = orbital['energy']
        row.label(text='Energy: {}'.format('—' if energy is None else format(energy, '.7g')))
        row.label(text='Occ: {:g}'.format(orbital['occupation']))
        row = entry.row()
        row.label(text='Spin: {}'.format(orbital['spin'] or '—'))
        row.label(text='Sym: {}'.format(orbital['symmetry'] or '—'))


class MB_PT_molecule(bpy.types.Panel):
    bl_label = 'Molecule'
    bl_idname = 'MB_PT_molecule'
    bl_space_type = 'VIEW_3D'
    bl_region_type = 'UI'
    bl_category = 'Molecule'

    def draw(self, context):
        layout = self.layout
        try:
            summary = molecule_summary(context.object)
        except (ValueError, KeyError, TypeError) as exc:
            layout.label(text=str(exc), icon='ERROR')
            return
        if summary is None:
            layout.label(text='Select an imported atom, bond, or surface')
            layout.label(text='Older molecules need to be reimported')
            return
        owner = summary['owner']
        layout.label(text=owner.name, icon='OUTLINER_OB_EMPTY')
        box = layout.box()
        box.label(text='Source: ' + os.path.basename(summary['source']))
        box.label(text=summary['source'])
        if not os.path.isfile(summary['source']):
            box.label(text='Source unavailable; saved scene remains editable', icon='INFO')
        box.label(text='Atoms: {}    Bonds: {}'.format(summary['atoms'], summary['bonds']))
        if summary['object_type'] != 'wireframe' and summary['visible_atoms'] != summary['atoms']:
            box.label(text='Displayed atoms: {}'.format(summary['visible_atoms']))
        box.label(text='Coordinates: ' + summary['units'])
        box.label(text='Molecular charge setting: {:+d}'.format(summary['charge']))
        box.label(text='{} ({})'.format(summary['representation'], summary['object_type']))

        atom = summary['selected_atom']
        if atom is not None:
            box = layout.box()
            box.label(text='Selected atom: {}'.format(atom['element'].capitalize()))
            box.label(text='Atom index: {}'.format(atom['index']))
            # A cube's atom-line charge can carry a nuclear charge.
            charge_label = 'Partial charge' if summary['fields'] is None else 'Cube atom charge'
            box.label(text='{}: {:g} e'.format(charge_label, atom['charge']))
            box.label(text='Formal charge: {:+d} e'.format(atom['formal_charge']))

        box = layout.box()
        box.label(text='Scientific data')
        orbitals = summary['orbitals']
        if orbitals:
            _draw_orbitals(box, owner, orbitals)
            labels = sorted({o['spin'] for o in orbitals if o['spin']})
            box.label(text='Spin labels: ' + (', '.join(labels) if labels else 'not supplied'))
            box.label(text='Density: total, alpha, beta, spin')
        fields = summary['fields']
        if fields is not None:
            box.label(text='Cube datasets: {}'.format(len(fields['datasets'])))
            box.label(text='Dataset IDs: ' + ', '.join(str(x) for x in fields['datasets']))
            box.label(text='Field type: ' + fields['field_type'])
        if not orbitals and fields is None:
            box.label(text='No orbital or grid data')
        box.label(text='Displayed surfaces: {}'.format(summary['surfaces']))
        if summary['object_type'] != 'wireframe':
            layout.operator('molecular_blender.restyle', text='Restyle molecule…')
            box = layout.box()
            box.label(text='Select two atoms or a bond')
            row = box.row(align=True)
            for order, label in (('1', 'Single'), ('2', 'Double'), ('3', 'Triple')):
                row.operator('molecular_blender.set_bond', text=label).order = order
            box.operator('molecular_blender.set_bond_tilt', text='Set strand tilt…')
            box.operator('molecular_blender.delete_bond')
            box.label(text='Manual orders are locked', icon='LOCKED')
            box.operator('molecular_blender.reperceive')


def register():
    bpy.types.Object.mb_show_orbitals = BoolProperty(name='Show saved orbitals', default=False)
    bpy.types.Object.mb_orbital_page = IntProperty(name='Orbital page', min=1,
                                                  get=_get_orbital_page, set=_set_orbital_page)
    bpy.utils.register_class(MB_PT_molecule)


def unregister():
    bpy.utils.unregister_class(MB_PT_molecule)
    del bpy.types.Object.mb_orbital_page
    del bpy.types.Object.mb_show_orbitals
