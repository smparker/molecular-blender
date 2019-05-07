#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Unit testing for molecular blender"""

import unittest
import sys
import re

import numpy as np

import cProfile

import molecular_blender as mb

def run():
    geo = mb.importers.molecule_from_file('zine.molden', {})
    orbitals = mb.orbitals.MOData.from_dict(geo)
    xyz = np.linspace(-10, 10, 50)

    mb.isosurfaces.molden_isosurface(orbitals.get_orbital(12), [0.01, -0.01], 0.01)

cProfile.run('run()', sort="time")
#run()
