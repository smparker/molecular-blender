"""Molecule sidebar edits use the saved owner, even across renames and selection changes."""
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

try:
    import bpy
except ImportError:
    bpy = None


@unittest.skipUnless(bpy, 'requires Blender')
class TestSidebarActions(unittest.TestCase):
    def setUp(self):
        from molecular_blender import scene_panel
        self.objects = set(bpy.data.objects)
        self.collections = set(bpy.data.collections)
        scene_panel.register()
        self.addCleanup(scene_panel.unregister)
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.path = Path(self.tmp.name) / 'water.molden'
        self.path.write_bytes(Path('examples/water-sto-3g.molden').read_bytes())

    def tearDown(self):
        for obj in set(bpy.data.objects) - self.objects:
            bpy.data.objects.remove(obj, do_unlink=True)
        for coll in set(bpy.data.collections) - self.collections:
            bpy.data.collections.remove(coll)

    def molecule(self):
        from molecular_blender.plotter import BlendMoleculeSetup
        model, _, _ = BlendMoleculeSetup(bpy.context, str(self.path), isovalues='', ignore_hydrogen=False)
        return model

    def select(self, obj):
        bpy.ops.object.select_all(action='DESELECT')
        obj.select_set(True)
        bpy.context.view_layer.objects.active = obj

    def test_indexed_field_mesh_keeps_winding_and_welds_only_sub_grid_seams(self):
        from molecular_blender.plotter import create_surface_mesh_data, index_field_triangles
        a, b, c, d = (0, 0, 0), (1, 0, 0), (0, 1, 0), (1, 1, 0)
        seam = (1.000001, 0, 0)
        triangles = [(a, b, c), (seam, d, c), (c, b, a), (b, seam, d),
                     [(x, y, .0001) for x, y, z in (a, c, b)]]
        verts, faces = index_field_triangles(triangles, .05)
        mesh = create_surface_mesh_data('Indexed field', verts, faces, weld_distance=None)
        try:
            self.assertEqual((len(mesh.vertices), len(mesh.polygons)), (7, 3))
            self.assertEqual([round(p.normal.z) for p in mesh.polygons], [1, 1, -1])
            self.assertAlmostEqual(sum(p.area for p in mesh.polygons), 1.5, places=5)
            self.assertFalse(mesh.validate())
        finally:
            bpy.data.meshes.remove(mesh)

    def test_legacy_surface_convention_survives_recomputing(self):
        from molecular_blender.scene_state import read_record
        from molecular_blender.scene_surfaces import recompute_settings
        from molecular_blender.options import validated_options
        model = self.molecule()
        self.assertEqual(read_record(model.scene_owner)['molden_convention'], 'standard')
        for legacy in (False, True):
            record = {'settings': validated_options({'orca_fix': not legacy})}
            group = {'field': 'homo', 'spin': 'all', 'sampling': {'orca_fix': legacy,
                     'resolution': .1, 'cumulative': False, 'isovalues': '.05'}}
            opts = validated_options(dict(record['settings'], **recompute_settings(record, group)))
            self.assertEqual(opts['molden_convention'], 'orca' if legacy else 'standard')

    def test_sidebar_add_and_toggle_without_touching_selection_or_other_molecules(self):
        from mathutils import Matrix, Vector
        from molecular_blender.plotter import compute_surface_vertex_sets, create_surface_meshes
        from molecular_blender.scene_state import read_record, owned_objects
        from molecular_blender.scene_surfaces import MB_OT_add_surface, surface_groups
        model = self.molecule()
        owner = model.scene_owner
        peer = self.molecule().scene_owner
        original_center = model.center_of_mass()
        owner.name = 'Moved water'
        owner.location += Vector((3, 2, 1))
        owner.rotation_euler.z = 0.7
        owner.scale = (2, 3, 1)
        bpy.context.view_layer.update()
        transform = owner.matrix_world @ Matrix.Translation(-original_center)
        atom = owned_objects(owner, 'atom')[0]
        self.select(atom)
        peer.select_set(True)
        saved_settings = read_record(owner)['settings']

        def finish_now(op, context, source, styler, settings, *, select_surfaces):
            self.assertFalse(select_surfaces)
            self.assertEqual(source.scene_owner, owner)
            sets, frames = compute_surface_vertex_sets(source.surface_input.load(settings), settings)
            create_surface_meshes(source, styler, context, settings, sets, frames, select_surfaces=False)
            for definition in sets:
                triangles = definition['triangles']
                if len(triangles):
                    mesh = definition['mesh']
                    context.view_layer.update()
                    point = transform @ Vector(triangles[0][0])
                    distance = min((point - mesh.matrix_world @ v.co).length for v in mesh.data.vertices)
                    self.assertLess(distance, .04)  # Mesh welding tolerance, scaled by the parent.
            return {'FINISHED'}

        with patch.object(MB_OT_add_surface, 'start_surface_job', finish_now):
            for field in ('orbital', 'density', 'orbital'):
                status = bpy.ops.molecular_blender.add_surface(volume=field, orbital='5', isovalues='0.05',
                                                              cumulative=False, resolution=.35, remesh=0)
                self.assertEqual(status, {'FINISHED'})
        self.assertEqual(set(bpy.context.selected_objects), {peer, atom})
        self.assertIsNone(peer.parent)
        groups = surface_groups(owner, read_record(owner))
        self.assertEqual([g['label'] for g in groups.values()], ['Orbital 5', 'Total density', 'Orbital 5'])
        self.assertEqual(read_record(owner)['settings'], saved_settings)
        self.assertEqual(read_record(peer)['surfaces'], [])
        self.path.unlink()
        key, first = next(iter(groups.items()))
        self.assertEqual(first['levels'], {-.05, .05})
        self.select(first['objects'][0])
        bpy.ops.molecular_blender.surface_visibility(group=key)
        self.assertEqual(bpy.context.active_object, owner)
        self.assertTrue(all(o.hide_render and o.hide_viewport for o in first['objects']))
        self.assertTrue(all(not o.hide_render for g in list(groups.values())[1:] for o in g['objects']))
        bpy.ops.molecular_blender.surface_visibility(group=key)
        self.assertTrue(all(not (o.hide_render or o.hide_viewport) for o in first['objects']))

    def test_missing_and_changed_sources_leave_the_scene_intact(self):
        from molecular_blender.scene_state import read_record
        from molecular_blender.scene_surfaces import surface_model
        owner = self.molecule().scene_owner
        self.select(owner)
        before = set(bpy.data.objects)
        original = self.path.read_text()
        for content, message in ((original.replace('1.857299866800e-01', '2.0'), 'coordinates'),
                                 (original.replace('-20.670655', '-21.0'), 'orbital metadata'),
                                 (None, 'source unavailable')):
            with self.subTest(content=message):
                if content is None:
                    self.path.unlink()
                else:
                    self.path.write_text(content)
                if content is None:
                    with self.assertRaisesRegex(RuntimeError, message):
                        bpy.ops.molecular_blender.add_surface()
                else:
                    record = read_record(owner)
                    request = surface_model(owner, record).surface_input
                    with self.assertRaisesRegex(ValueError, message):
                        request.load(record['settings'])
                self.assertEqual(set(bpy.data.objects), before)
                self.assertEqual(read_record(owner)['surfaces'], [])

    def test_row_style_and_removal_isolate_phase_group_and_linked_copies(self):
        from molecular_blender.plotter import BlendMoleculeSetup, draw_surfaces
        from molecular_blender.scene_edit import restyle_molecule
        from molecular_blender.scene_state import read_record
        from molecular_blender.scene_surfaces import surface_groups
        model, styler, options = BlendMoleculeSetup(bpy.context, str(self.path), orbital='5',
            isovalues='.05,.1', cumulative=False, resolution=.35, remesh=0, recycle_materials=True)
        owner = model.scene_owner
        for _ in range(2):
            draw_surfaces(model, styler, bpy.context, options)
        self.path.unlink()
        (key, group), (_, other) = surface_groups(owner, read_record(owner)).items()
        plus, minus = group['phases'][True], group['phases'][False]
        meshes = [obj.data for obj in group['objects']]
        copy = plus[0].copy()
        copy.parent = None
        model.scene_collection.objects.link(copy)
        color = (.2, .4, .6)
        original = copy.active_material
        shared_mesh = copy.data
        minus_before = tuple(minus[0].mb_surface_color)
        plus[0].mb_surface_color = color
        plus[0].mb_surface_opacity = .25
        plus[0].mb_surface_color = (.1, .3, .5)
        self.assertEqual(original, copy.active_material)
        self.assertAlmostEqual(copy.mb_surface_opacity, 1)
        self.assertEqual(tuple(minus[0].mb_surface_color), minus_before)
        for obj in plus:
            mat = obj.active_material
            shader = next(n for n in mat.node_tree.nodes if n.type in {'BSDF_PRINCIPLED', 'EMISSION'})
            socket = shader.inputs['Color' if shader.type == 'EMISSION' else 'Base Color']
            for value, expected in zip(socket.default_value, (.1, .3, .5)):
                self.assertAlmostEqual(value, expected)
        restyle_molecule(owner, atom_scale=1.2)
        for obj in group['objects']:
            self.assertAlmostEqual(obj.mb_surface_opacity, .25)
        self.assertEqual([o.data for o in group['objects']], meshes)
        for obj in other['objects']:
            self.assertAlmostEqual(obj.mb_surface_opacity, 1)
        names = [obj.name for obj in group['objects']]
        collections = {coll.name for obj in group['objects'] for coll in obj.users_collection}
        self.select(plus[0])
        self.assertEqual(bpy.ops.molecular_blender.remove_surface(group=key), {'FINISHED'})
        self.assertEqual(bpy.context.active_object, owner)
        self.assertTrue(all(name not in bpy.data.objects for name in names))
        self.assertTrue(all(name not in bpy.data.collections for name in collections))
        record = read_record(owner)
        self.assertEqual({s['id'] for s in record['surfaces']}, other['ids'])
        self.assertEqual(copy.data, shared_mesh)

    def test_fine_orbital_mesh_welds_only_repeated_intersections(self):
        from molecular_blender.plotter import BlendMoleculeSetup, create_surface_meshes
        model, styler, options = BlendMoleculeSetup(bpy.context, str(self.path), orbital='5',
            isovalues='.05', cumulative=False, resolution=.05, remesh=0)
        a, b, c, d = (0., 0., 0.), (.004, 0., 0.), (0., .005, 0.), (0., 0., .006)
        sets = [dict(name='5', isovalue=.05, triangles=[(a,c,b), (a,b,d), (b,c,d), (c,a,d)])]
        create_surface_meshes(model, styler, bpy.context, options, sets, {}, select_surfaces=False)
        mesh = sets[0]['mesh'].data
        self.assertEqual((len(mesh.vertices), len(mesh.edges), len(mesh.polygons)), (4, 6, 4))

    def test_recompute_keeps_identity_presentation_and_other_groups(self):
        from molecular_blender.plotter import BlendMoleculeSetup, compute_surface_vertex_sets, draw_surfaces
        from molecular_blender.scene_state import read_record
        from molecular_blender.scene_surfaces import MB_OT_recompute_surface, surface_groups, set_surface_visibility
        model, styler, options = BlendMoleculeSetup(bpy.context, str(self.path), orbital='5',
            isovalues='.05', cumulative=False, resolution=.32, remesh=0)
        owner = model.scene_owner
        self.select(owner)
        for _ in range(2):
            draw_surfaces(model, styler, bpy.context, options)
        (key, group), (other_key, other) = surface_groups(owner, read_record(owner)).items()
        minus, plus = group['phases'][False][0], group['phases'][True][0]
        plus.name = 'custom lobe'
        plus.mb_surface_color = (.1, .8, .3)
        minus.mb_surface_color = (.2, .3, .9)
        plus.mb_surface_opacity = .37
        set_surface_visibility(group['objects'], False)
        before = {o: o.data.as_pointer() for o in group['objects']}
        other_meshes = [o.data for o in other['objects']]
        def finish(op, ctx, model, styler, options, **kwargs):
            op._owner, op._styler, op._options = owner, styler, options
            sets, frames = compute_surface_vertex_sets(model.surface_input.load(options), options)
            op.apply_surface_result(ctx, [v for v in sets if len(v['triangles'])], frames)
            return {'FINISHED'}
        with patch.object(MB_OT_recompute_surface, 'start_surface_job', finish):
            self.assertEqual(bpy.ops.molecular_blender.recompute_surface(
                group=key, isovalues='.05,.13', resolution=.24), {'FINISHED'})
            first, second = surface_groups(owner, read_record(owner)).items()
            self.assertEqual((first[0], second[0]), (key, other_key))
            updated = first[1]
            self.assertEqual(updated['levels'], {-.13, -.05, .05, .13})
            self.assertEqual(updated['sampling']['isovalues'], '0.05, 0.13')
            self.assertAlmostEqual(updated['sampling']['resolution'], .24)
            for obj in updated['objects']:
                self.assertAlmostEqual(obj.mb_surface_opacity, .37)
                self.assertTrue(obj.hide_viewport and obj.hide_render)
                self.assertAlmostEqual(obj.mb_surface_color[1], .8 if obj in updated['phases'][True] else .3)
            self.assertTrue(all(o.data.as_pointer() != mesh for o, mesh in before.items()))
            self.assertIn(bpy.data.objects['custom lobe'], updated['objects'])
            self.assertEqual(bpy.ops.molecular_blender.recompute_surface(group=key, isovalues='.05'), {'FINISHED'})
            self.assertEqual(set(surface_groups(owner, read_record(owner))[key]['objects']), set(before))
            self.assertEqual([o.data for o in other['objects']], other_meshes)

    def test_failed_recompute_retains_original_meshes_and_saved_settings(self):
        from molecular_blender.plotter import BlendMoleculeSetup, draw_surfaces, create_surface_mesh_data
        from molecular_blender.scene_state import read_record
        from molecular_blender.scene_surfaces import replace_surface_geometry, surface_groups, recompute_settings
        model, styler, options = BlendMoleculeSetup(bpy.context, str(self.path), orbital='5',
            isovalues='.05', cumulative=False, resolution=.32, remesh=0)
        draw_surfaces(model, styler, bpy.context, options)
        owner = model.scene_owner
        record = read_record(owner)
        for entry in record['surfaces']:
            entry.pop('sampling')
        key, group = next(iter(surface_groups(owner, record).items()))
        self.assertEqual(recompute_settings(record, group)['isovalues'], '0.05')
        self.assertFalse(recompute_settings(record, group)['cumulative'])
        baseline = (owner['mb_molecule'], set(bpy.data.meshes), set(bpy.data.objects))
        triangles = [[(0,0,0), (1,0,0), (1,1,0)]]
        defs = [dict(isovalue=x, triangles=triangles) for x in (.02, -.02)]
        def fail_second(name, verts, faces, weld_distance):
            if set(bpy.data.meshes) != baseline[1]:
                raise ValueError('second mesh failed')
            return create_surface_mesh_data(name, verts, faces, weld_distance)
        with patch('molecular_blender.scene_surfaces.create_surface_mesh_data', fail_second):
            with self.assertRaisesRegex(ValueError, 'second mesh failed'):
                replace_surface_geometry(bpy.context, owner, record, key, group, defs, options, styler)
        self.assertEqual((owner['mb_molecule'], set(bpy.data.meshes), set(bpy.data.objects)), baseline)

    def test_frontier_jumps_use_saved_populations_spin_and_import_order(self):
        from molecular_blender.scene_state import read_record, write_record
        from test_scene_panel import panel_labels
        owner = self.molecule().scene_owner
        self.path.unlink()
        self.select(owner)
        self.assertEqual(bpy.ops.molecular_blender.jump_orbital(kind='lumo'), {'FINISHED'})
        self.assertTrue(owner.mb_show_orbitals)
        self.assertEqual(owner.mb_orbital_page, 2)
        self.assertIn('Orbital 6 · LUMO', panel_labels(owner))
        record = read_record(owner)
        for o in record['orbitals']:
            o['spin'] = 'Beta' if o['index'] in {1, 6, 7} else 'Alpha'
            o['occupation'] = 1 if o['index'] <= 5 else 0
        record['orbitals'][5]['energy'], record['orbitals'][6]['energy'] = 5., 1.
        record['settings']['orbital_spin'] = 'beta'
        write_record(owner, record)
        self.assertEqual(owner.mb_browser_spin, 'beta')
        bpy.ops.molecular_blender.jump_orbital(kind='lumo')
        self.assertIn('Orbital 7 · LUMO', panel_labels(owner))
        bpy.ops.molecular_blender.jump_orbital(kind='homo')
        self.assertEqual(owner.mb_orbital_page, 1)
        self.assertIn('Orbital 1 · HOMO', panel_labels(owner))
        owner.mb_browser_spin = 'alpha'
        self.assertNotIn('Orbital 1 · HOMO', panel_labels(owner))
        self.assertEqual(bpy.ops.molecular_blender.jump_orbital(kind='lumo'), {'CANCELLED'})
        bpy.ops.molecular_blender.jump_orbital(kind='homo')
        self.assertIn('Orbital 5 · HOMO', panel_labels(owner))

    def test_element_and_neighbor_selection_read_current_bonds(self):
        from molecular_blender.scene_state import owned_objects, read_record, write_record
        owner = self.molecule().scene_owner
        peer = self.molecule().scene_owner
        atoms = owned_objects(owner, 'atom')
        hydrogens = [o for o in atoms if o['mb_element'] == 'h']
        oxygen, = [o for o in atoms if o['mb_element'] == 'o']
        self.select(hydrogens[0])
        peer.select_set(True)
        bpy.ops.molecular_blender.select_atoms(match='element')
        self.assertEqual(set(bpy.context.selected_objects), set(hydrogens))
        bpy.ops.molecular_blender.select_atoms(match='neighbors')
        self.assertEqual(set(bpy.context.selected_objects), {hydrogens[0], oxygen})
        self.assertEqual(bpy.context.active_object, hydrogens[0])
        oxygen.hide_set(True)
        bpy.ops.molecular_blender.select_atoms(match='fragment')
        self.assertEqual(set(bpy.context.selected_objects), set(hydrogens))
        self.assertEqual(bpy.context.active_object, hydrogens[0])
        oxygen.hide_set(False)
        self.select(oxygen)
        bpy.ops.molecular_blender.select_atoms(match='neighbors')
        self.assertEqual(set(bpy.context.selected_objects), set(atoms))
        record = read_record(owner)
        record['bonds'] = []
        write_record(owner, record)
        bpy.ops.molecular_blender.select_atoms(match='neighbors')
        self.assertEqual(set(bpy.context.selected_objects), {oxygen})
        bpy.ops.molecular_blender.select_atoms(match='fragment')
        self.assertEqual(set(bpy.context.selected_objects), {oxygen})

    def test_index_ranges_validate_before_changing_selection(self):
        from molecular_blender.scene_state import ID, owned_objects, read_record
        owner = self.molecule().scene_owner
        peer = self.molecule().scene_owner
        objects = {o[ID]: o for o in owned_objects(owner, 'atom')}
        atoms = [objects[a['id']] for a in read_record(owner)['atoms']]
        self.select(owner)
        peer.select_set(True)
        ops = bpy.ops.molecular_blender
        self.assertEqual(ops.select_atom_indices(indices='3, 1 - 2, 1'), {'FINISHED'})
        self.assertEqual(set(bpy.context.selected_objects), set(atoms))
        self.assertEqual(bpy.context.active_object, atoms[2])
        for invalid in ('1,4', '3-1', '1-999999999999', '1,,2', ''):
            self.assertEqual(ops.select_atom_indices(indices=invalid), {'CANCELLED'})
            self.assertEqual(set(bpy.context.selected_objects), set(atoms))
        atoms[0].hide_set(True)
        self.assertEqual(ops.select_atom_indices(indices='1,2'), {'FINISHED'})
        self.assertEqual(bpy.context.selected_objects, [atoms[1]])
        self.assertEqual(ops.select_atom_indices(indices='1'), {'CANCELLED'})
        self.assertEqual(bpy.context.active_object, atoms[1])

    def test_manifest_positions_and_orbitals_use_the_same_blender_frames(self):
        import json
        from molecular_blender.plotter import BlendMoleculeSetup, compute_surface_vertex_sets
        from molecular_blender.scene_state import owned_objects
        manifest = self.path.with_suffix('.json')
        manifest.write_text(json.dumps({'molecules': [
            {'filename': 'skipped.molden', 'frame': 0}, {'filename': 'water.molden', 'frame': 1},
            {'filename': 'skipped.molden', 'frame': 2}, {'filename': 'water.molden', 'frame': 4}]}))
        model, _, options = BlendMoleculeSetup(bpy.context, str(manifest), orbital='5', keystride=2, frame_start=1, frame_step=2,
                                               isovalues='.05', cumulative=False, resolution=.4, remesh=0)
        atom = owned_objects(model.scene_owner, 'atom')[0]
        times = {p.co.x for curve in atom.animation_data.action.fcurves for p in curve.keyframe_points}
        self.assertEqual(times, {3, 9})
        definitions, frames = compute_surface_vertex_sets(model, options)
        self.assertEqual(set(frames), times)
        self.assertTrue(all(frames.values()))
        self.assertNotEqual(frames[3], frames[9])
        self.assertTrue(all(len(d['triangles']) for d in definitions))

    def test_relink_validates_and_changes_only_the_root_source(self):
        from molecular_blender.scene_state import read_record
        from molecular_blender.scene_surfaces import surface_model
        owner = self.molecule().scene_owner
        peer = self.molecule().scene_owner
        self.select(owner)
        before = read_record(owner)
        moved = self.path.with_name('moved.molden')
        self.path.rename(moved)
        bad = self.path.with_name('other.molden')
        original = moved.read_text()
        for content, message in [(original.replace('1.857299866800e-01', '2.0'), 'coordinates'),
                                 (original.replace('-20.670655', '-21.0'), 'orbital metadata'),
                                 ('invalid file', 'Molden Format')]:
            bad.write_text(content)
            with self.assertRaisesRegex(RuntimeError, message):
                bpy.ops.molecular_blender.locate_source(filepath=str(bad))
            self.assertEqual(read_record(owner), before)
        self.assertEqual(bpy.ops.molecular_blender.locate_source(filepath=str(moved)), {'FINISHED'})
        after = read_record(owner)
        self.assertEqual(after, dict(before, source=str(moved)))
        self.assertEqual(read_record(peer)['source'], before['source'])
        self.assertEqual(len(surface_model(owner, after).surface_input.load(after['settings']).atoms), 3)

    def test_solo_switches_and_restores_individual_states_even_after_adding_a_surface(self):
        from molecular_blender.plotter import BlendMoleculeSetup, draw_surfaces
        from molecular_blender.scene_state import read_record
        from molecular_blender.scene_surfaces import SOLO_GROUP, surface_groups, visibility_state
        model, styler, options = BlendMoleculeSetup(bpy.context, str(self.path), orbital='5,6',
            isovalues='.05', cumulative=False, resolution=.35, remesh=0)
        owner = model.scene_owner
        draw_surfaces(model, styler, bpy.context, options)
        groups = surface_groups(owner, read_record(owner))
        (a, first), (b, second) = groups.items()
        first['objects'][0].hide_render = True
        second['objects'][0].hide_viewport = True
        states = {o: visibility_state(o) for g in groups.values() for o in g['objects']}
        self.select(first['objects'][0])
        ops = bpy.ops.molecular_blender
        self.assertEqual(ops.solo_surface(group=b), {'FINISHED'})
        self.assertEqual(bpy.context.active_object, owner)
        self.assertTrue(all(visibility_state(o) == 3 for o in first['objects']))
        self.assertTrue(all(visibility_state(o) == 0 for o in second['objects']))
        draw_surfaces(model, styler, bpy.context, options)
        extra = [o for k, g in surface_groups(owner, read_record(owner)).items() if k not in groups for o in g['objects']]
        self.assertTrue(extra and all(visibility_state(o) == 3 for o in extra))
        self.assertEqual(ops.solo_surface(group=a), {'FINISHED'})
        self.assertTrue(all(visibility_state(o) == 0 for o in first['objects']))
        self.assertTrue(all(visibility_state(o) == 3 for o in second['objects']))
        self.assertEqual(ops.solo_surface(group=a), {'FINISHED'})
        self.assertEqual(states, {o: visibility_state(o) for o in states})
        self.assertTrue(all(visibility_state(o) == 0 for o in extra))
        ops.solo_surface(group=b)
        ops.remove_surface(group=b)
        self.assertNotIn(SOLO_GROUP, owner)
        self.assertTrue(all(visibility_state(o) == states[o] for o in first['objects']))

    def test_png_series_uses_camera_and_restores_settings_and_solo_after_failure(self):
        from molecular_blender import surface_export
        from molecular_blender.plotter import BlendMoleculeSetup, draw_surfaces
        from molecular_blender.scene_state import read_record
        from molecular_blender.scene_surfaces import SOLO_GROUP, surface_groups, visibility_state
        model, styler, options = BlendMoleculeSetup(bpy.context, str(self.path), orbital='5,6',
            isovalues='.05', cumulative=False, resolution=.35, remesh=0)
        draw_surfaces(model, styler, bpy.context, options)
        owner = model.scene_owner
        self.select(owner)
        scene = bpy.context.scene
        camera = scene.camera
        if camera is None:
            camera = bpy.data.objects.new('Test camera', bpy.data.cameras.new('Test camera'))
            scene.collection.objects.link(camera)
            scene.camera = camera
        groups = surface_groups(owner, read_record(owner))
        keys = list(groups)
        keys.reverse()
        bpy.ops.molecular_blender.solo_surface(group=keys[1])
        states = {o: visibility_state(o) for g in groups.values() for o in g['objects']}
        directory = Path(self.tmp.name)
        stem = bpy.path.clean_name(owner.name + '-' + groups[keys[0]]['label'], replace='-')
        occupied = directory / (stem + '.png')
        occupied.write_bytes(b'keep me')
        render = scene.render
        def settings():
            return (render.filepath, render.film_transparent, render.use_multiview, render.use_sequencer,
                    render.image_settings.file_format, render.image_settings.color_depth, render.image_settings.color_mode)
        original_settings = settings()
        calls = []
        def simulated_render(**kwargs):
            key = keys[len(calls)]
            self.assertEqual(scene.camera, camera)
            self.assertTrue(render.film_transparent)
            self.assertEqual(render.image_settings.color_mode, 'RGBA')
            for other, group in groups.items():
                self.assertTrue(all(o.hide_render == (other != key) for o in group['objects']))
            calls.append(Path(render.filepath))
            calls[-1].write_bytes(b'png')
            return {'FINISHED'} if len(calls) == 1 else {'CANCELLED'}
        with patch.object(surface_export, 'bpy', wraps=bpy) as api:
            api.path.clean_name = bpy.path.clean_name
            api.ops.render.render.side_effect = simulated_render
            with self.assertRaisesRegex(RuntimeError, 'did not finish'):
                surface_export.export_surface_pngs(bpy.context, owner, groups, keys, directory)
        self.assertEqual(occupied.read_bytes(), b'keep me')
        self.assertEqual(calls[0].read_bytes(), b'png')
        self.assertFalse(calls[1].exists())
        self.assertEqual(settings(), original_settings)
        self.assertEqual(owner[SOLO_GROUP], keys[1])
        self.assertTrue(all(visibility_state(o) == state for o, state in states.items()))
        self.assertEqual(bpy.context.selected_objects, [owner])

        from mathutils import Quaternion
        cameras = set(bpy.data.cameras)
        with patch.object(surface_export, 'bpy', wraps=bpy) as api:
            api.data = bpy.data
            api.path.clean_name = bpy.path.clean_name
            api.ops.render.render.return_value = {'CANCELLED'}
            with self.assertRaisesRegex(RuntimeError, 'did not finish'):
                surface_export.export_surface_pngs(bpy.context, owner, groups, keys, directory, fit_rotation=Quaternion())
        self.assertEqual(scene.camera, camera)
        self.assertEqual(set(bpy.data.cameras), cameras)
        self.assertTrue(all(visibility_state(o) == state for o, state in states.items()))

    def test_legacy_animated_groups_preserve_keyframes_and_save_visibility(self):
        from molecular_blender.scene_state import tag, read_record, write_record
        from molecular_blender.scene_surfaces import surface_groups, surface_visible, set_surface_visibility
        model = self.molecule()
        owner = model.scene_owner
        record = read_record(owner)
        record['settings']['plot_type'] = 'animate'
        for i, name in enumerate(('homo-1', 'homo-1-   1')):
            obj = bpy.data.objects.new(name, None)
            obj.parent = owner
            model.scene_collection.objects.link(obj)
            tag(obj, owner, 'surface', str(i))
            record['surfaces'].append(dict(field=name, id=str(i), isovalue=.05))
            for frame in (1, 2):
                obj.hide_render = obj.hide_viewport = (frame != i+1)
                obj.keyframe_insert(data_path='hide_viewport', frame=frame)
                obj.keyframe_insert(data_path='hide_render', frame=frame)
        write_record(owner, record)
        group, = surface_groups(owner, record).values()
        self.assertEqual(group['label'], 'Orbital HOMO-1')
        objects = group['objects']
        set_surface_visibility(objects, False)
        old_frame = bpy.context.scene.frame_current
        self.addCleanup(bpy.context.scene.frame_set, old_frame)
        bpy.context.scene.frame_set(2)
        self.assertFalse(any(surface_visible(obj) for obj in objects))
        self.assertTrue(all(o.users_collection[0].hide_render for o in objects))
        path = str(Path(self.tmp.name) / 'hidden.blend')
        bpy.data.libraries.write(path, {model.scene_collection})
        with bpy.data.libraries.load(path, link=False) as (_, data):
            data.collections = [model.scene_collection.name]
        restored = data.collections[0]
        bpy.context.scene.collection.children.link(restored)
        copy = next(o for o in restored.objects if 'mb_molecule' in o)
        copied_group, = surface_groups(copy, read_record(copy)).values()
        copied = copied_group['objects']
        self.assertFalse(any(surface_visible(o) for o in copied))
        set_surface_visibility(copied, True)
        self.assertTrue(all(o.users_collection[0] == restored for o in copied))
        self.assertFalse(any(surface_visible(o) for o in objects))
        set_surface_visibility(objects, True)
        for frame in (1, 2):
            bpy.context.scene.frame_set(frame)
            self.assertEqual([o.hide_render for o in objects], [frame != 1, frame != 2])
        self.select(owner)
        set_surface_visibility(objects, False)
        gate_name = objects[0].users_collection[0].name
        self.assertEqual(bpy.ops.molecular_blender.remove_surface(group='homo-1'), {'FINISHED'})
        self.assertNotIn(gate_name, bpy.data.collections)
        self.assertEqual(read_record(owner)['surfaces'], [])
        self.assertEqual(len(surface_groups(copy, read_record(copy))['homo-1']['objects']), 2)
