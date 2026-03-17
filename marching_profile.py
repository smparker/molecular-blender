#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Unit testing for molecular blender"""

import unittest
import sys
import re

import numpy as np

import cProfile

import molecular_blender as mb


def orbital_isosurface():
    geo = mb.importers.molecule_from_file('examples/zine.molden', {})
    orbitals = mb.orbitals.MOData.from_dict(geo)

    mb.isosurfaces.molden_isosurface(orbitals.get_orbital(12), [0.01, -0.01], 0.1)

def density_isosurface():
    geo = mb.importers.molecule_from_file('examples/zine.molden', {})
    orbitals = mb.orbitals.MOData.from_dict(geo)

    mb.isosurfaces.molden_isosurface(orbitals.get_density("density"), [0.01], 0.1)

def run():
    density_isosurface()

cProfile.run('run()', sort="time")
