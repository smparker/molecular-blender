"""Verify multiple bonds using Blender's evaluated extrusion geometry."""
import unittest
from pathlib import Path

try:
    import bpy
    import bmesh
    from mathutils import Vector
    from molecular_blender.plotter import (set_bond_order, PlotMolecule,
                                          default_options, plot_prep, make_materials,
                                          BlendMoleculeSetup)
    from molecular_blender.importers import molecule_from_file
    from molecular_blender.molecule import Atom, Bond, Molecule
    from molecular_blender.stylers import get_styler
    HAS_BPY = True
except ImportError:
    HAS_BPY = False


@unittest.skipUnless(HAS_BPY, 'requires Blender')
class TestBondOrders(unittest.TestCase):
    def test_spatial_neighbors_match_exhaustive_static_and_trajectory_pairs(self):
        import random
        from molecular_blender.molecule import Snapshot, bonded_pairs
        rng = random.Random(719)
        atoms = [Atom(rng.choice(('h', 'c', 'n', 's')), [rng.uniform(-7, 7) for _ in range(3)], i)
                 for i in range(120)]
        for atom in atoms:
            atom.trajectory = [Snapshot(atom.position + Vector([rng.uniform(-.7, .7) for _ in range(3)]),
                                        0, (0, 0, 0)) for _ in range(3)]
        for static in (True, False):
            expected = []
            for i, a in enumerate(atoms):
                for j, b in enumerate(atoms[:i]):
                    bond = Bond(a, b)
                    if static:
                        exists = bond.is_bonded()
                    else:
                        exists = any(bond.is_bonded((x.position-y.position).length)
                                     for x, y in zip(a.trajectory, b.trajectory))
                    if exists:
                        expected.append((i, j))
            self.assertEqual(bonded_pairs(atoms, static), expected)

    def test_animated_sdf_keeps_explicit_bonds_in_dynamic_distance_mode(self):
        from molecular_blender.scene_state import owned_objects
        model, _, _ = BlendMoleculeSetup(bpy.context, 'examples/water_bend.json', ignore_hydrogen=False,
                                        isovalues='', animate_bonds='dynamic', keystride=2)
        self.assertEqual(model.frame_numbers, [0, 4, 8])
        self.assertEqual(len(model.bonds), 2)
        curves = owned_objects(model.scene_owner, 'bond')
        self.assertTrue(curves)
        self.assertTrue(all(o.animation_data is None or not o.animation_data.action for o in curves))

    def setUp(self):
        self.objects = set(bpy.data.objects)
        self.curves = set(bpy.data.curves)

    def tearDown(self):
        for obj in set(bpy.data.objects) - self.objects:
            bpy.data.objects.remove(obj, do_unlink=True)
        for curve in set(bpy.data.curves) - self.curves:
            if curve.users == 0:
                bpy.data.curves.remove(curve)

    def bond_curve(self, end):
        curve = bpy.data.curves.new('test_bond', 'CURVE')
        curve.dimensions = '3D'
        curve.bevel_depth = 0.2
        curve.use_fill_caps = True
        spline = curve.splines.new('BEZIER')
        spline.bezier_points.add(1)
        for point, co in zip(spline.bezier_points, [(0, 0, 0), end]):
            point.co = co
            point.handle_left_type = point.handle_right_type = 'AUTO'
        obj = bpy.data.objects.new('test_bond', curve)
        bpy.context.collection.objects.link(obj)
        return obj

    def strand_centers(self, obj):
        bpy.context.view_layer.update()
        evaluated = obj.evaluated_get(bpy.context.evaluated_depsgraph_get())
        mesh = evaluated.to_mesh()
        try:
            # Blender emits cap vertices separately from the tube walls.
            welded = bmesh.new()
            welded.from_mesh(mesh)
            bmesh.ops.remove_doubles(welded, verts=list(welded.verts), dist=1e-6)
            welded.to_mesh(mesh)
            welded.free()
            neighbors = {v.index: set() for v in mesh.vertices}
            for edge in mesh.edges:
                a, b = edge.vertices
                neighbors[a].add(b)
                neighbors[b].add(a)
            remaining = set(neighbors)
            centers = []
            while remaining:
                todo = [remaining.pop()]
                component = set(todo)
                while todo:
                    for index in neighbors[todo.pop()] - component:
                        component.add(index)
                        remaining.discard(index)
                        todo.append(index)
                centers.append(sum((mesh.vertices[i].co for i in component), Vector()) / len(component))
            return centers
        finally:
            evaluated.to_mesh_clear()

    def test_extruded_strands_and_orientation(self):
        # Includes both vertical singularities, a reversed bond, and an oblique
        # separation with a component along the bond that must be projected out.
        for end, direction in [((3, 0, 0), (0, 1, 0)),
                               ((0, 3, 0), (0, 0, 1)),
                               ((0, 0, 3), (1, 1, 0)),
                               ((0, 0, -3), (1, 1, 0)),
                               ((-3, 0, 0), (0, 1, 1)),
                               ((1, 2, 3), (2, -1, 1))]:
            for order in (2, 3):
                with self.subTest(end=end, direction=direction, order=order):
                    obj = self.bond_curve(end)
                    set_bond_order(obj, order, direction)
                    centers = self.strand_centers(obj)
                    self.assertEqual(len(centers), order)
                    tangent = Vector(end).normalized()
                    desired = Vector(direction)
                    desired = (desired - tangent * desired.dot(tangent)).normalized()
                    centers.sort(key=lambda c: c.dot(desired))
                    for a, b in zip(centers, centers[1:]):
                        self.assertAlmostEqual((b - a).length, 0.34, places=4)
                        self.assertAlmostEqual((b - a).normalized().dot(desired), 1, places=5)
                    self.assertEqual(obj.type, 'CURVE')
                    self.assertEqual(len(obj.data.splines), 1)

    def test_reset_and_invalid_direction(self):
        obj = self.bond_curve((3, 0, 0))
        with self.assertRaises(ValueError):
            set_bond_order(obj, 2, (1, 0, 0))
        with self.assertRaises(ValueError):
            set_bond_order(obj, 4)
        set_bond_order(obj, 2, (0, 1, 0))
        profile = obj.data.bevel_object.name
        set_bond_order(obj, 1)
        self.assertNotIn(profile, bpy.data.objects)
        self.assertEqual(len(self.strand_centers(obj)), 1)
        self.assertAlmostEqual(obj.data.bevel_depth, 0.2)

    def test_split_profiles_have_no_collapsed_faces(self):
        # Bevel-factor trimming at a sampled midpoint used to duplicate an end
        # ring on the first half. Its zero-area quads caused dark color streaks.
        for order in (2, 3):
            for bounds in ((0, 0.5), (0.5, 1)):
                with self.subTest(order=order, bounds=bounds):
                    obj = self.bond_curve((1, 2, 0))
                    obj.data.resolution_u = 2
                    obj.data.bevel_factor_start, obj.data.bevel_factor_end = bounds
                    set_bond_order(obj, order, (2, -1, 0))
                    bpy.context.view_layer.update()
                    evaluated = obj.evaluated_get(bpy.context.evaluated_depsgraph_get())
                    mesh = evaluated.to_mesh()
                    try:
                        self.assertTrue(mesh.polygons)
                        self.assertTrue(all(p.area > 1e-10 for p in mesh.polygons))
                    finally:
                        evaluated.to_mesh_clear()

    def test_split_colors_and_hooks(self):
        atoms = [Atom('c', (0, 0, 0), 0), Atom('o', (1.3, 0, 0), 1)]
        mol = Molecule('multiple_bond_test', atoms)
        bond = Bond(*atoms, order=2, direction=(0, 1, 0))
        mol.bonds = [bond]
        options = default_options({'plot_style': 'sticks', 'hook_atoms': 'on', 'draw_bond_orders': True})
        make_materials(mol, get_styler(options), options)
        plot_prep(mol, options)
        PlotMolecule(bpy.context, mol, options)
        self.assertEqual(len(bond.name), 2)
        halves = [bpy.data.objects[name] for name in bond.name]
        self.assertNotEqual(halves[0].data.materials[0], halves[1].data.materials[0])
        for obj, bounds in zip(halves, [(0, 0.5), (0.5, 1)]):
            self.assertEqual(len(self.strand_centers(obj)), 2)
            self.assertEqual(len(obj.modifiers), 2)
            self.assertEqual((obj.data.bevel_factor_start, obj.data.bevel_factor_end), bounds)
        before = self.strand_centers(halves[1])
        bpy.data.objects[atoms[1].name].location.x += 1
        after = self.strand_centers(halves[1])
        self.assertGreater(sum(c.x for c in after) - sum(c.x for c in before), 1)

    def test_mixed_element_bond_radii(self):
        for style in ('bs', 'fixedbs', 'sticks'):
            for universal in (False, True):
                for symbols in (('c', 'o'), ('o', 'c'), ('c', 'h')):
                    with self.subTest(style=style, universal=universal, symbols=symbols):
                        atoms = [Atom(s, (1.3*i, 0, 0), i) for i, s in enumerate(symbols)]
                        mol = Molecule('bond_radius_test', atoms)
                        bond = Bond(*atoms, style='universal' if universal else '')
                        mol.bonds = [bond]
                        options = default_options({'plot_style': style,
                                                   'universal_bonds': universal,
                                                   'bond_thickness': 0.17})
                        make_materials(mol, get_styler(options), options)
                        plot_prep(mol, options)
                        PlotMolecule(bpy.context, mol, options)
                        expected = (min(a.element.vdw for a in atoms) * 0.1
                                    if style == 'bs' and not universal else 0.17)
                        for name in bond.name:
                            self.assertAlmostEqual(bpy.data.objects[name].data.bevel_depth,
                                                   expected)

    def test_charge_defaults_for_python_import(self):
        path = Path(__file__).resolve().parents[1] / 'examples/bond_orders.xyz'
        molecule, _, options = BlendMoleculeSetup(
            bpy.context, str(path), charges='scale', ignore_hydrogen=False)
        self.assertEqual(options['charge_offset'], 0.0)
        self.assertEqual(options['charge_factor'], 1.0)
        self.assertEqual(molecule.chgoff, 0.0)
        self.assertEqual(molecule.chgfac, 1.0)
        for atom in molecule.atoms:
            self.assertIn(atom.plus_charge, bpy.data.objects)
            self.assertIn(atom.neg_charge, bpy.data.objects)
        explicit = default_options({'charge_offset': 0.4, 'charge_factor': 2.0})
        self.assertEqual(explicit['charge_offset'], 0.4)
        self.assertEqual(explicit['charge_factor'], 2.0)

    def test_density_volume_selector_ignores_orbital_text(self):
        from molecular_blender.plotter import compute_surface_vertex_sets
        path = Path(__file__).resolve().parents[1] / 'examples/water-sto-3g.molden'
        molecule, _, options = BlendMoleculeSetup(
            bpy.context, str(path), volume='density', orbital='not-an-orbital,also-invalid',
            isovalues='0.05', cumulative=False, resolution=0.25)
        progress = {}
        surfaces, _ = compute_surface_vertex_sets(molecule, options, progress)
        self.assertEqual([s['name'] for s in surfaces], ['density'])
        self.assertGreater(len(surfaces[0]['triangles']), 0)
        self.assertEqual(progress['total'], 1)
        self.assertEqual(progress['completed'], 1)

    def test_sdf_import_to_extruded_geometry(self):
        path = Path(__file__).resolve().parents[1] / 'examples' / 'bond_orders.sdf'
        for style in ('sticks', 'fixedbs'):
            with self.subTest(style=style):
                molecule, _, options = BlendMoleculeSetup(
                    bpy.context, str(path), plot_style=style, ignore_hydrogen=False,
                    hook_atoms='on', bond_thickness=0.14, draw_bond_orders=True)
                self.assertEqual(options['plot_type'], 'frame')
                self.assertEqual([b.order for b in molecule.bonds], [2, 1, 2, 1, 3, 1, 1, 1])
                for bond in molecule.bonds:
                    for name in bond.name:
                        obj = bpy.data.objects[name]
                        self.assertEqual(obj.type, 'CURVE')
                        self.assertEqual(len(obj.modifiers), 2)
                        self.assertEqual(len(self.strand_centers(obj)), bond.order)
                    if bond.order == 2:
                        # Planar fixture: separation stays in XY, perpendicular
                        # to the bond, rather than sticking out of its plane.
                        self.assertIsNotNone(bond.direction)
                        self.assertAlmostEqual(bond.direction.z, 0)

    def test_explicit_graph_overrides_distances(self):
        path = Path(__file__).resolve().parents[1] / 'examples' / 'bond_orders.sdf'
        data = molecule_from_file(str(path), {})
        for i, atom in enumerate(data['atoms']):
            atom['position'] = [10*i, 0, 0]
        molecule = Molecule.from_dict('explicit_graph', data)
        options = default_options({'plot_type': 'frame'})
        molecule.determine_bonding(options)
        self.assertEqual(len(molecule.bonds), 8)  # even far beyond distance cutoffs
        molecule.determine_bonding(options)
        self.assertEqual(len(molecule.bonds), 8)  # no duplicates on a second call
        options['bonds'] = False
        molecule.determine_bonding(options)
        self.assertEqual(molecule.bonds, [])
        data['bonds'] = []
        for atom in data['atoms']:
            atom['position'] = [0, 0, 0]
        molecule = Molecule.from_dict('empty_graph', data)
        options['bonds'] = True
        molecule.determine_bonding(options)
        self.assertEqual(molecule.bonds, [])  # explicit empty table is authoritative

    def test_xyz_perception_to_extruded_geometry(self):
        root = Path(__file__).resolve().parents[1]
        data = molecule_from_file(str(root / 'examples/bond_orders.sdf'), {})
        expected = {tuple(sorted(b['atoms'])): b['order'] for b in data['bonds']}
        for style in ('sticks', 'fixedbs'):
            with self.subTest(style=style):
                molecule, _, _ = BlendMoleculeSetup(
                    bpy.context, str(root / 'examples/bond_orders.xyz'),
                    plot_style=style, perceive_bond_orders=True, draw_bond_orders=True,
                    bond_thickness=0.14, ignore_hydrogen=True)
                actual = {tuple(sorted((b.iatom.index, b.jatom.index))): b.order
                          for b in molecule.bonds}
                self.assertEqual(actual, expected)
                self.assertTrue(all(a.formal_charge == 0 for a in molecule.atoms))
                for bond in molecule.bonds:
                    if bond.iatom.hidden or bond.jatom.hidden:
                        continue
                    for name in bond.name:
                        self.assertEqual(len(self.strand_centers(bpy.data.objects[name])), bond.order)

    def test_perception_is_opt_in_and_keeps_explicit_orders(self):
        root = Path(__file__).resolve().parents[1]
        options = default_options({'plot_type': 'frame', 'perceive_bond_orders': False})
        data = molecule_from_file(str(root / 'examples/bond_orders.xyz'), options)
        molecule = Molecule.from_dict('no_perception', data)
        molecule.determine_bonding(options)
        self.assertTrue(all(b.order == 1 for b in molecule.bonds))
        # A nonsensical requested charge must not override a supplied SDF table.
        options.update(perceive_bond_orders=True, molecular_charge=99)
        data = molecule_from_file(str(root / 'examples/bond_orders.sdf'), options)
        molecule = Molecule.from_dict('explicit_orders', data)
        molecule.determine_bonding(options)
        self.assertEqual([b.order for b in molecule.bonds], [2, 1, 2, 1, 3, 1, 1, 1])
        # Union/dynamic trajectory graphs are not suitable for a single valence assignment.
        molecule = Molecule('trajectory', [Atom('c', (0, 0, 0), 0)])
        options.update(plot_type='animate', animate_bonds='dynamic', molecular_charge=0)
        with self.assertRaisesRegex(ValueError, 'static'):
            molecule.determine_bonding(options)

    def test_perceive_fused_aromatic_dimer(self):
        root = Path(__file__).resolve().parents[1]
        options = default_options({'plot_type': 'frame', 'perceive_bond_orders': True,
                                   'find_aromatic': True})
        data = molecule_from_file(str(root / 'examples/tetracene_dimer.xyz'), options)
        molecule = Molecule.from_dict('aromatic_perception', data)
        molecule.determine_bonding(options)
        self.assertEqual(len(molecule.rings), 8)
        self.assertTrue(any(b.order == 2 for b in molecule.bonds))
        for atom in molecule.atoms:
            valence = sum(b.order for b in molecule.bonds if atom in (b.iatom, b.jatom))
            self.assertEqual(valence, 4 if atom.element.symbol == 'c' else 1)

    def test_perception_total_charge_and_failure(self):
        atoms = [Atom('n', (0, 0, 0), 0, charge=-0.6)]
        for i, pos in enumerate(((.6, .6, .6), (.6, -.6, -.6),
                                  (-.6, .6, -.6), (-.6, -.6, .6)), 1):
            atoms.append(Atom('h', pos, i, charge=0.2))
        molecule = Molecule('ammonium', atoms)
        options = default_options({'plot_type': 'frame', 'perceive_bond_orders': True})
        with self.assertRaisesRegex(ValueError, 'No valid'):
            molecule.determine_bonding(options)
        self.assertTrue(all(a.formal_charge == 0 for a in atoms))
        self.assertTrue(all(b.order == 1 for b in molecule.bonds))
        options['molecular_charge'] = 1
        molecule.determine_bonding(options)
        self.assertEqual([a.formal_charge for a in atoms], [1, 0, 0, 0, 0])
        self.assertEqual(len(molecule.bonds), 4)
        self.assertEqual(atoms[0].charge, -0.6)  # partial charges are independent

    def test_bond_order_drawing_defaults_off(self):
        path = Path(__file__).resolve().parents[1] / 'examples/bond_orders.sdf'
        molecule, _, options = BlendMoleculeSetup(bpy.context, str(path), plot_style='sticks')
        self.assertFalse(options['draw_bond_orders'])
        self.assertTrue(any(b.order > 1 for b in molecule.bonds))
        for bond in molecule.bonds:
            if bond.iatom.hidden or bond.jatom.hidden:
                continue
            for name in bond.name:
                obj = bpy.data.objects[name]
                self.assertEqual(obj.data.bevel_mode, 'ROUND')
                self.assertIsNone(obj.data.bevel_object)
                self.assertEqual(len(self.strand_centers(obj)), 1)
        for atom in molecule.atoms:
            if not atom.hidden:
                self.assertEqual(tuple(bpy.data.objects[atom.name].scale), (1, 1, 1))

    def test_symbolic_orbital_requests_cannot_wrap(self):
        from molecular_blender.plotter import compute_surface_vertex_sets
        root = Path(__file__).resolve().parents[1]
        data = molecule_from_file(str(root / 'examples/water-sto-3g.molden'), {})
        molecule = Molecule.from_dict('orbital_validation', data)
        for expression in (f'homo-{molecule.orbitals.nocc}', 'lumo+999', 'homo-1oops'):
            options = default_options({'orbital': expression, 'isosurfaces': True,
                                       'isovalues': [.05], 'cumulative': False})
            with self.subTest(expression=expression), self.assertRaises(ValueError):
                compute_surface_vertex_sets(molecule, options)

    def test_explicit_density_modes_and_zero_spin(self):
        from molecular_blender.plotter import compute_surface_vertex_sets
        path = Path(__file__).resolve().parents[1] / 'examples/water-sto-3g.molden'
        for mode in ('density', 'alpha', 'beta', 'spin'):
            molecule, _, options = BlendMoleculeSetup(
                bpy.context, str(path), volume=mode, orbital='ignored',
                isovalues='0.05', cumulative=True, resolution=0.3)
            progress = {}
            surfaces, _ = compute_surface_vertex_sets(molecule, options, progress)
            if mode == 'spin':
                self.assertEqual(surfaces, [])  # closed-shell zero spin is valid
                self.assertEqual(set(options['isovalues']), {-.05, .05})
            else:
                self.assertTrue(surfaces)
                self.assertEqual(surfaces[0]['name'], mode)
            self.assertEqual(progress, {'total': 1, 'completed': 1})

    def test_spin_frontier_control_selects_the_matching_input_orbital(self):
        from unittest.mock import patch
        from molecular_blender.orbitals import MOData
        from molecular_blender.plotter import compute_surface_vertex_sets
        root = Path(__file__).resolve().parents[1]
        data = molecule_from_file(str(root / 'examples/water-sto-3g.molden'), {})
        molecule = Molecule.from_dict('spin_frontier', data)
        base = molecule.orbitals
        molecule.orbitals = MOData(base.shells, base.coeff[:4], 2, [1.,0.,1.,0.],
                                   ['Alpha','Alpha','Beta','Beta'], energies=[-1.,1.,-2.,2.])
        options = default_options({'orbital': 'homo', 'orbital_spin': 'beta',
                                   'isosurfaces': True, 'isovalues': [.1], 'cumulative': False})
        with patch('molecular_blender.plotter.molden_isosurface', return_value=[]) as generate:
            compute_surface_vertex_sets(molecule, options)
        selected = generate.call_args[0][0]
        expected = molecule.orbitals.get_orbital(3)
        self.assertAlmostEqual(selected.value(0., .1, .2), expected.value(0., .1, .2))

    def test_importer_bond_order_toggle(self):
        from types import SimpleNamespace
        from molecular_blender import MolecularBlender
        registered = MolecularBlender.is_registered
        if not registered:
            bpy.utils.register_class(MolecularBlender)
        try:
            # Blender does not allow direct construction of Operator instances.
            # Read the registered UI defaults and exercise its options builder.
            rna = getattr(bpy.ops, 'import').molecular_blender.get_rna_type()
            operator = SimpleNamespace(**{p.identifier: p.default for p in rna.properties
                                          if p.identifier != 'dataset_choice' and hasattr(p, 'default')})
            self.assertFalse(operator.draw_bond_orders)
            options = MolecularBlender._build_options(operator)
            self.assertFalse(options['draw_bond_orders'])
            self.assertFalse(options['perceive_bond_orders'])
            operator.draw_bond_orders = True
            options = MolecularBlender._build_options(operator)
            self.assertTrue(options['draw_bond_orders'])
            self.assertTrue(options['perceive_bond_orders'])
        finally:
            if not registered:
                bpy.utils.unregister_class(MolecularBlender)
