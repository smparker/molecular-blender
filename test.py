#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Unit testing for molecular blender"""

import unittest
import sys
import re

import molecular_blender as mb
import numpy as np

import bpy

class TestBlendSticksRings(unittest.TestCase):
    """Test Suite for plotting with sticks and drawing rings"""
    def setUp(self):
        """Setup function"""
        self.file = 'examples/tetracene_dimer.xyz'
        self.atom_re = re.compile(r"^tetracene_dimer_c\d+$")
        self.bond_re = re.compile(r"^tetracene_dimer_c\d+-c\d+$")
        self.context = bpy.context

    def test_sticks_and_rings_import(self):
        """Test importer for sticks/rings"""
        mb.BlendMolecule(self.context, self.file,
                         plot_style="sticks",
                         find_aromatic=True)

        nrings = sum(
            ["tetracene_dimer_ring" in x.name for x in bpy.data.objects])
        self.assertEqual(nrings, 8, 'Incorrect number of rings.')

        ncarbon = sum([self.atom_re.match(x.name)
                       is not None for x in bpy.data.objects])
        self.assertEqual(ncarbon, 36, 'Incorrect number of Carbon atoms.')

        ncc = sum([self.bond_re.match(x.name)
                   is not None for x in bpy.data.objects])
        self.assertEqual(ncc, 42, 'Incorrect number of Carbon-Carbon bonds.')

        nobj = sum(["tetracene_dimer" in x.name for x in bpy.data.objects])
        self.assertEqual(nobj, 87, 'Incorrect number of objects created.')


class TestAnimateBallAndSticks(unittest.TestCase):
    """Test Suite for animation"""
    def setUp(self):
        """Setup function"""
        self.file = 'examples/benzopentyl_twist.xyz'
        self.context = bpy.context

    def test_balls_and_sticks(self):
        """Test importer for animation involving balls/sticks"""
        mb.BlendMolecule(self.context, self.file,
                         plot_style="fixedbs",
                         plot_type="animate",
                         keystride=3)

        nact = sum(
            ["benzopentyl_twist" in action.name for action in bpy.data.actions])
        self.assertEqual(nact, 11, 'Incorrect number of actions registered.')

        action = bpy.data.actions["benzopentyl_twist_c5Action"]
        nfcurves = len(action.fcurves)
        self.assertEqual(
            nfcurves, 3, 'Incorrect number of fcurves registered.')

        nkeyframe_points = len(action.fcurves[0].keyframe_points)
        self.assertEqual(nkeyframe_points, 4,
                         'Incorrect number of keyframe points registered.')

class TestAnimateSticksAndCharges(unittest.TestCase):
    """Test suite for animation including charges"""
    def setUp(self):
        """Setup function"""
        self.file = 'examples/tio2_excited_state.xyz'
        self.context = bpy.context

    def test_sticks_and_charges(self):
        """Test importer for animation involving sticks/charges"""
        mb.BlendMolecule(self.context, self.file,
                         plot_style="sticks",
                         plot_type="animate",
                         animate_bonds="dynamic",
                         ignore_hydrogen=False,
                         keystride=2,
                         charges="scale",
                         charge_offset=0.95,
                         charge_factor=4.0)

        nneg = sum(["_neg" in x.name for x in bpy.data.objects if "tio2_excited_state" in x.name])
        self.assertEqual(nneg, 42, 'Incorrect number of negative spheres found')

        npos = sum(["_plus" in x.name for x in bpy.data.objects if "tio2_excited_state" in x.name])
        self.assertEqual(npos, 42, 'Incorrect number of positive spheres found')

        # pick some objects that have disappearing bonds and count the number of hidden frames
        objs =  [ bpy.data.objects[x] for x in ["tio2_excited_state_h18-o4_a",
            "tio2_excited_state_h18-o14_a"]]

        nhidden = [ 0 for o in objs ]
        for frame in range(120):
            self.context.scene.frame_set(frame)

            for i, o in enumerate(objs):
                nhidden[i] += o.hide_viewport
        self.assertEqual(nhidden, [71,47], 'Incorrect number of hidden frames found')


class TestXYZRead(unittest.TestCase):
    """Test Suite for reading XYZ"""
    def setUp(self):
        """Setup function"""
        self.geo = mb.importers.molecule_from_file(
            'examples/tetracene_dimer.xyz', {})

    def test_read_xyz(self):
        """Test to read XYZ"""
        self.assertEqual(len(self.geo["atoms"]), 60)


class TestCubeRead(unittest.TestCase):
    """Test Suite for reading cube files"""
    def setUp(self):
        """Setup function"""
        self.geo = mb.importers.molecule_from_file('examples/water.cube', {})
        dxyz = [np.linalg.norm(self.geo["volume"]["axes"][:, i])
                for i in range(3)]
        self.dvol = np.prod(dxyz)

    def test_total_density(self):
        """Test integrate density"""
        vol = self.geo["volume"]["data"]

        total_density = np.sum(vol * vol) * self.dvol
        self.assertAlmostEqual(total_density, 0.141058531145, places=4)


class TestMoldenRead(unittest.TestCase):
    """Test Suite for reading molden files"""
    def setUp(self):
        """Setup function"""
        self.geo = mb.importers.molecule_from_file(
            'examples/water-sto-3g.molden', {})
        self.orbitals = mb.orbitals.MOData.from_dict(self.geo)
        self.xyz = np.linspace(-10, 10, 50)
        self.dvol = (self.xyz[1] - self.xyz[0])**3

    def test_orbital_norm(self):
        """Test to verify orbital norm"""
        orb = self.orbitals.get_orbital(2)

        integration = 0.0
        for x in self.xyz:
            for y in self.xyz:
                for z in self.xyz:
                    phi = orb.value(x, y, z)
                    integration += phi * phi
        integration *= self.dvol

        self.assertAlmostEqual(integration, 0.98273143395905316, places=6)

    def test_orbital_overlap(self):
        """Test to verify orbital overlaps"""
        iorb = self.orbitals.get_orbital(3)
        jorb = self.orbitals.get_orbital(4)

        integration = 0.0
        for x in self.xyz:
            for y in self.xyz:
                for z in self.xyz:
                    phi_i = iorb.value(x, y, z)
                    phi_j = jorb.value(x, y, z)
                    integration += phi_i * phi_j
        integration *= self.dvol

        self.assertAlmostEqual(integration, -0.0010284808822276, places=6)

class TestMoldenFunctions(unittest.TestCase):
    """Test Suite for reading molden files"""
    def setUp(self):
        """Setup function"""
        self.geo = mb.importers.molecule_from_file(
            'examples/fake-atom.molden', {})
        self.orbitals = mb.orbitals.MOData.from_dict(self.geo)
        self.xyz = np.linspace(-30, 30, 100, dtype=np.float32)
        self.dvol = (self.xyz[1] - self.xyz[0])**3

    def test_orbital_norms(self):
        """Test to verify orbital norm"""

        for i in range(1,25+1):
            orb = self.orbitals.get_orbital(i)
            integration = np.sum(orb.box_values(self.xyz, self.xyz, self.xyz)**2) * self.dvol
            self.assertAlmostEqual(integration, 1.0, places=1)

    def test_cumulative_orbital_isovalue(self):
        """Test to verify cumulative isovalue"""
        orb = self.orbitals.get_orbital(1)

        isovals = orb.isovalue_containing_proportion([0.5, 0.9], resolution=0.25, box=[[-3, -3, -3], [3, 3, 3]])
        self.assertAlmostEqual(isovals[0], 0.12, places=2)
        self.assertAlmostEqual(isovals[1], 0.048, places=2)

    def test_density_norm(self):
        """Test to verify density norm"""
        density = self.orbitals.get_density()
        integration = np.sum(density.box_values(self.xyz, self.xyz, self.xyz)) * self.dvol
        self.assertAlmostEqual(integration, 4.0, places=1)

    def test_cumulative_density_isovalue(self):
        """Test to verify cumulative isovalue"""
        orb = self.orbitals.get_density()

        isovals = orb.isovalue_containing_proportion([0.5, 0.9], resolution=0.25, box=([-3, -3, -3], [3, 3, 3]))
        self.assertAlmostEqual(isovals[0], 0.10, places=2)
        self.assertAlmostEqual(isovals[1], 0.010, places=2)

class TestTurbomoleFunctions(unittest.TestCase):
    """Testing turbomole related functions"""
    def setUp(self):
        """Setup function"""
        self.geo = mb.importers.molecule_from_file(
            'examples/phenanth_polar.json', {})
        self.orbitals = mb.orbitals.MOData.from_dict(self.geo)
        self.xyz = np.linspace(-30, 30, 100, dtype=np.float32)
        self.dvol = (self.xyz[1] - self.xyz[0])**3

    def test_polarizability(self):
        """Test to verify polarizability"""
        polar = self.orbitals.get_density("polar 0")

        pxyz = polar.box_values(self.xyz, self.xyz, self.xyz)
        # the polarization density should integrate to zero
        integrate = np.sum(pxyz) * self.dvol
        self.assertLess(np.abs(integrate), 0.5)

        # integrating polarization * r should give the polarizability
        alphax = np.zeros([3], dtype=np.float32)
        for i in range(3):
            # first sum along the other two axes
            other_axes = tuple([ x for x in range(3) if x != i ])
            partial = np.sum(pxyz, axis=other_axes)
            # then integrate along the remaining axis
            integrate = np.sum(partial * self.xyz) * self.dvol
            alphax[i] = integrate

        self.assertAlmostEqual(alphax[0], 206.0, delta=2.0)
        self.assertAlmostEqual(alphax[1], 33.0, delta=1.0)
        self.assertAlmostEqual(alphax[2], -0.00012, delta=0.05)

def blender_argv(argv):
    """Processes argv to process the Blender specific parts"""
    out = argv
    if "blender" in sys.argv[0] and "--python" in sys.argv:
        # probably called through blender
        out = [argv[argv.index('--python') + 1]]
        if '--' in argv:
            out.extend(argv[argv.index('--') + 1:])

    return out


# When called through Blender, sys.argv includes the blender command. The
# code block here is constructing argv to match what unittest expects.
if __name__ == "__main__":
    import sys
    print(sys.version)
    unittest.main(argv=blender_argv(sys.argv))
