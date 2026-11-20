"""Inspect the molecule owning the current Blender selection."""
import os

import bpy
from bpy.props import BoolProperty, EnumProperty, FloatProperty, FloatVectorProperty, IntProperty

from .importers import file_format
from .orbitals import frontier_index
from .scene_state import ID, RECORD, ROLE, read_record, resolve_owner, owned_objects
from .scene_surfaces import (MB_OT_add_surface, MB_OT_remove_surface, MB_OT_surface_visibility,
                             set_surface_color, set_surface_opacity, surface_color, surface_groups,
                             surface_opacity, surface_visible)
from .scene_view import MB_OT_frame_molecule
from .surface_operator import DENSITIES, ORBITAL_SPINS

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
    groups = surface_groups(owner, record)
    return dict(owner=owner, source=record['source'], units=record['units'],
                selected_atom=atom,
                atoms=len(record['atoms']), visible_atoms=len(owned_objects(owner, 'atom')),
                bonds=len(record['bonds']), surfaces=sum(len(g['objects']) for g in groups.values()),
                groups=groups,
                representation=REPRESENTATIONS[settings['plot_style']],
                object_type=settings['object_type'], charge=settings['molecular_charge'],
                orbitals=record['orbitals'], fields=record['fields'])


def _page_count(orbitals):
    return max(1, (len(orbitals) + ORBITALS_PER_PAGE - 1) // ORBITALS_PER_PAGE)


def _level_labels(levels):
    labels = []
    for magnitude in sorted({abs(value) for value in levels}):
        if magnitude and magnitude in levels and -magnitude in levels:
            labels.append('±{:.4g}'.format(magnitude))
        else:
            value = magnitude if magnitude in levels else -magnitude
            labels.append(format(value, '.4g'))
    return labels


def _get_orbital_page(owner):
    return owner.get('mb_orbital_page', 1)


def _set_orbital_page(owner, page):
    orbitals = read_record(owner)['orbitals']
    owner['mb_orbital_page'] = max(1, min(page, _page_count(orbitals)))


def _get_browser_spin(owner):
    if 'mb_browser_spin' in owner:
        return owner['mb_browser_spin']
    spin = read_record(owner)['settings']['orbital_spin'] if RECORD in owner else 'all'
    return next(i for i, item in enumerate(ORBITAL_SPINS) if item[0] == spin)


def _set_browser_spin(owner, value):
    owner['mb_browser_spin'] = value
    owner['mb_frontier_index'] = 0


class MB_OT_jump_orbital(bpy.types.Operator):
    bl_idname = 'molecular_blender.jump_orbital'
    bl_label = 'Jump to Orbital'
    bl_description = 'Show the orbital for the chosen spin in the saved orbital list'

    kind: EnumProperty(items=(('homo', 'HOMO', ''), ('lumo', 'LUMO', '')))

    def execute(self, context):
        owner = resolve_owner(context.object)
        if owner is None:
            return {'CANCELLED'}
        orbitals = read_record(owner)['orbitals']
        spins = [o['spin'] for o in orbitals]
        channel = owner.mb_browser_spin
        try:
            index = frontier_index([o['occupation'] for o in orbitals], spins if any(spins) else None,
                                   [o['energy'] for o in orbitals], self.kind,
                                   None if channel == 'all' else channel)
        except ValueError as exc:
            self.report({'WARNING'}, str(exc))
            return {'CANCELLED'}
        owner.mb_show_orbitals = True
        owner.mb_orbital_page = (index - 1) // ORBITALS_PER_PAGE + 1
        owner['mb_frontier_index'] = index
        owner['mb_frontier_kind'] = self.kind.upper()
        return {'FINISHED'}


def _draw_orbitals(box, owner, orbitals, molden, available):
    box.prop(owner, 'mb_show_orbitals', text='Orbitals ({})'.format(len(orbitals)),
             icon='TRIA_DOWN' if owner.mb_show_orbitals else 'TRIA_RIGHT', emboss=False)
    if not owner.mb_show_orbitals:
        return
    box.prop(owner, 'mb_browser_spin')
    jump = box.row(align=True)
    for kind in ('homo', 'lumo'):
        jump.operator('molecular_blender.jump_orbital', text='Go to ' + kind.upper()).kind = kind
    pages = _page_count(orbitals)
    if pages > 1:
        box.prop(owner, 'mb_orbital_page', text='Page (1–{})'.format(pages))
    start = (min(owner.mb_orbital_page, pages) - 1) * ORBITALS_PER_PAGE
    for orbital in orbitals[start:start + ORBITALS_PER_PAGE]:
        entry = box.column(align=True)
        row = entry.row(align=True)
        label = 'Orbital {}'.format(orbital['index'])
        if orbital['index'] == owner.get('mb_frontier_index'):
            row.label(text=label + ' · ' + owner['mb_frontier_kind'], icon='BOOKMARKS')
        else:
            row.label(text=label)
        if molden:
            action = row.row(align=True)
            action.enabled = available
            op = action.operator('molecular_blender.add_surface', text='Add Surface…')
            op.volume = 'orbital'
            op.orbital = str(orbital['index'])
        row = entry.row()
        energy = orbital['energy']
        row.label(text='Energy: {}'.format('—' if energy is None else format(energy, '.7g')))
        row.label(text='Occ: {:g}'.format(orbital['occupation']))
        row = entry.row()
        row.label(text='Spin: {}'.format(orbital['spin'] or '—'))
        row.label(text='Sym: {}'.format(orbital['symmetry'] or '—'))


class MB_OT_select_atoms(bpy.types.Operator):
    bl_idname = 'molecular_blender.select_atoms'
    bl_label = 'Select atoms'
    bl_description = 'Select matching visible atoms on the active molecule, including the starting atom'
    bl_options = {'UNDO'}

    match: EnumProperty(items=(('element', 'Same element', ''), ('neighbors', 'Bonded neighbors', ''),
                               ('fragment', 'Bonded fragment', 'All atoms connected through bonds')))

    @classmethod
    def poll(cls, context):
        obj = context.active_object
        return context.mode == 'OBJECT' and obj is not None and obj.get(ROLE) == 'atom'

    def execute(self, context):
        atom = context.active_object
        owner = resolve_owner(atom)
        if owner is None:
            return {'CANCELLED'}
        record = read_record(owner)
        if self.match == 'element':
            ids = {a['id'] for a in record['atoms'] if a['element'] == atom['mb_element']}
        else:
            adjacency = {}
            for bond in record['bonds']:
                left, right = bond['atoms']
                adjacency.setdefault(left, set()).add(right)
                adjacency.setdefault(right, set()).add(left)
            ids = {atom[ID]}
            pending = [atom[ID]]
            while pending:
                new = adjacency.get(pending.pop(), set()) - ids
                ids.update(new)
                if self.match == 'fragment':
                    pending.extend(new)
        selectable = set(context.selectable_objects)
        matches = {obj for obj in owned_objects(owner, 'atom') if obj in selectable and obj[ID] in ids}
        for obj in context.selected_objects:
            if obj not in matches:
                obj.select_set(False)
        for obj in matches:
            obj.select_set(True)
        return {'FINISHED'}


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
        layout.operator('molecular_blender.frame_molecule', icon='VIEWZOOM')
        box = layout.box()
        box.label(text='Source: ' + os.path.basename(summary['source']))
        box.label(text=summary['source'])
        available = os.path.isfile(summary['source'])
        if not available:
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
            box.operator('molecular_blender.select_atoms', text='Select Same Element').match = 'element'
            box.operator('molecular_blender.select_atoms', text='Select Bonded Neighbors').match = 'neighbors'
            box.operator('molecular_blender.select_atoms', text='Select Bonded Fragment').match = 'fragment'

        box = layout.box()
        box.label(text='Scientific data')
        orbitals = summary['orbitals']
        if orbitals:
            molden = file_format(summary['source']) == 'molden'
            _draw_orbitals(box, owner, orbitals, molden, available)
            labels = sorted({o['spin'] for o in orbitals if o['spin']})
            box.label(text='Spin labels: ' + (', '.join(labels) if labels else 'not supplied'))
            if molden:
                row = box.row(align=True)
                row.enabled = available
                row.prop(owner, 'mb_density_field', text='')
                row.operator('molecular_blender.add_surface', text='Add Surface…').volume = owner.mb_density_field
        fields = summary['fields']
        if fields is not None:
            box.label(text='Cube datasets: {}'.format(len(fields['datasets'])))
            box.label(text='Dataset IDs: ' + ', '.join(str(x) for x in fields['datasets']))
            box.label(text='Field type: ' + fields['field_type'])
            if 'units' in fields:
                box.label(text='Field units: ' + fields['units'])
        if not orbitals and fields is None:
            box.label(text='No orbital or grid data')
        if summary['groups']:
            box = layout.box()
            box.label(text='Surfaces')
            for key, group in summary['groups'].items():
                visible = any(surface_visible(obj) for obj in group['objects'])
                entry = box.column(align=True)
                row = entry.row(align=True)
                row.operator('molecular_blender.surface_visibility', text=group['label'],
                             icon='HIDE_OFF' if visible else 'HIDE_ON').group = key
                row.operator('molecular_blender.remove_surface', text='', icon='TRASH').group = key
                colors = entry.row(align=True)
                for sign, label in ((True, 'Positive'), (False, 'Negative')):
                    phase = group['phases'].get(sign)
                    if phase:
                        colors.prop(phase[0], 'mb_surface_color', text=label)
                entry.prop(group['objects'][0], 'mb_surface_opacity', text='Opacity', slider=True)
                levels = _level_labels(group['levels'])
                for i in range(0, len(levels), 3):
                    prefix = 'Isovalues  ' if i == 0 else ''
                    entry.label(text=prefix + ', '.join(levels[i:i+3]))
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


_CLASSES = (MB_OT_select_atoms, MB_OT_surface_visibility, MB_OT_add_surface, MB_OT_remove_surface,
            MB_OT_jump_orbital, MB_OT_frame_molecule, MB_PT_molecule)


def register():
    bpy.types.Object.mb_surface_color = FloatVectorProperty(name='Surface color', subtype='COLOR',
        size=3, min=0., max=1., get=surface_color, set=set_surface_color, options=set())
    bpy.types.Object.mb_surface_opacity = FloatProperty(name='Surface opacity', min=0., max=1.,
        get=surface_opacity, set=set_surface_opacity, options=set())
    bpy.types.Object.mb_browser_spin = EnumProperty(name='Orbital Spin', items=ORBITAL_SPINS,
        get=_get_browser_spin, set=_set_browser_spin)
    bpy.types.Object.mb_density_field = EnumProperty(name='Density', items=DENSITIES, default='density')
    bpy.types.Object.mb_show_orbitals = BoolProperty(name='Show saved orbitals', default=False)
    bpy.types.Object.mb_orbital_page = IntProperty(name='Orbital page', min=1,
                                                  get=_get_orbital_page, set=_set_orbital_page)
    for cls in _CLASSES:
        bpy.utils.register_class(cls)


def unregister():
    for cls in reversed(_CLASSES):
        bpy.utils.unregister_class(cls)
    del bpy.types.Object.mb_browser_spin
    del bpy.types.Object.mb_surface_color
    del bpy.types.Object.mb_surface_opacity
    del bpy.types.Object.mb_density_field
    del bpy.types.Object.mb_orbital_page
    del bpy.types.Object.mb_show_orbitals
