#!/usr/bin/env python

import molecular_blender as mb
import numpy as np

import bpy

import unittest
import sys
import re

class TestBlendSticksRings(unittest.TestCase):
    def setUp(self):
        self.file = 'examples/tetracene_dimer.xyz'
        self.atom_re = re.compile(r"^tetracene_dimer_Carbon\d+$")
        self.bond_re = re.compile(r"^tetracene_dimer_Carbon\d+-Carbon\d+$")
        self.context = bpy.context

    def test_sticks_and_rings_import(self):
        mb.BlendMolecule(self.context, self.file,
            plot_style = "sticks",
            find_aromatic = True)

        nrings = sum([ "tetracene_dimer_ring" in x.name for x in bpy.data.objects ])
        self.assertEqual(nrings, 8, 'Incorrect number of rings.')

        nCarbon = sum([ self.atom_re.match(x.name) is not None for x in bpy.data.objects ])
        self.assertEqual(nCarbon, 36, 'Incorrect number of Carbon atoms.')

        nCC = sum([ self.bond_re.match(x.name) is not None for x in bpy.data.objects ])
        self.assertEqual(nCC, 42, 'Incorrect number of Carbon-Carbon bonds.')

        nobj = sum([ "tetracene_dimer" in x.name for x in bpy.data.objects ])
        self.assertEqual(nobj, 88, 'Incorrect number of objects created.')

class TestAnimateBallAndSticks(unittest.TestCase):
    def setUp(self):
        self.file = 'examples/benzopentyl_twist.xyz'
        self.context = bpy.context

    def test_balls_and_sticks(self):
        mb.BlendMolecule(self.context, self.file,
            plot_style = "fixedbs",
            plot_type = "animate",
            keystride = 3)

        nact = sum(["benzopentyl_twist" in action.name for action in bpy.data.actions])
        self.assertEqual(nact, 11, 'Incorrect number of actions registered.')

        action = bpy.data.actions["benzopentyl_twist_Carbon5Action"]
        nfcurves = len(action.fcurves)
        self.assertEqual(nfcurves, 3, 'Incorrect number of fcurves registered.')

        nkeyframe_points = len(action.fcurves[0].keyframe_points)
        self.assertEqual(nkeyframe_points, 4, 'Incorrect number of keyframe points registered.')

class TestXYZRead(unittest.TestCase):
    def setUp(self):
        self.geo = mb.importers.molecule_from_file('examples/tetracene_dimer.xyz', {})

    def test_read_xyz(self):
        self.assertEqual(len(self.geo["atoms"]), 60)

class TestCubeRead(unittest.TestCase):
    def setUp(self):
        self.geo = mb.importers.molecule_from_file('examples/water.cube', {})
        dxyz = [ np.linalg.norm(self.geo["volume"]["axes"][:,i]) for i in range(3) ]
        self.dV = np.prod(dxyz)

    def test_total_density(self):
        vol = self.geo["volume"]["data"]

        total_density = np.sum(vol*vol) * self.dV
        self.assertAlmostEqual(total_density, 0.141058531145, places=4)

class TestMoldenRead(unittest.TestCase):
    def setUp(self):
        self.geo = mb.importers.molecule_from_file('examples/water-sto-3g.molden', {})
        self.orbitals = mb.orbitals.MOData.from_dict(self.geo)
        self.xyz = np.linspace(-10, 10, 50)
        self.dV = (self.xyz[1] - self.xyz[0])**3

    def test_orbital_norm(self):
        orb = self.orbitals.get_orbital(2)

        integration = 0.0
        for x in self.xyz:
            for y in self.xyz:
                for z in self.xyz:
                    phi = orb.value(x,y,z)
                    integration += phi*phi
        integration *= self.dV

        self.assertAlmostEqual(integration, 0.98273143395905316, places=6)

    def test_orbital_overlap(self):
        iorb = self.orbitals.get_orbital(3)
        jorb = self.orbitals.get_orbital(4)

        integration = 0.0
        for x in self.xyz:
            for y in self.xyz:
                for z in self.xyz:
                    phi_i = iorb.value(x,y,z)
                    phi_j = jorb.value(x,y,z)
                    integration += phi_i * phi_j
        integration *= self.dV

        self.assertAlmostEqual(integration, -0.0010284808822276, places=6)

# When called through Blender, sys.argv includes the blender command. The
# code block here is constructing argv to match what unittest expects.
if __name__ == "__main__":
    if "blender" in sys.argv[0] and "--python" in sys.argv:
        # probably called through blender
        argv = [ sys.argv[sys.argv.index('--python')+1] ]
        if '--' in sys.argv:
            argv.extend( sys.argv[sys.argv.index('--')+1:] )
    else:
        argv = sys.argv

    unittest.main(argv=argv)
