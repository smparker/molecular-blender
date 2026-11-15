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
        options = default_options({'plot_style': 'sticks', 'hook_atoms': 'on'})
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

    def test_sdf_import_to_extruded_geometry(self):
        path = Path(__file__).resolve().parents[1] / 'examples' / 'bond_orders.sdf'
        for style in ('sticks', 'fixedbs'):
            with self.subTest(style=style):
                molecule, _, options = BlendMoleculeSetup(
                    bpy.context, str(path), plot_style=style, ignore_hydrogen=False,
                    hook_atoms='on', bond_thickness=0.14)
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
