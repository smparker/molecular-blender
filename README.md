# Molecular Blender [![Build](https://github.com/smparker/molecular-blender/actions/workflows/blender-plugin.yml/badge.svg?branch=master)](https://github.com/smparker/molecular-blender/actions/workflows/blender-plugin.yml)
An import plugin specialized for .xyz files of molecules used in, for example,
quantum chemistry

## Capabilities

- imports using linked objects making rapid changes to the aesthetics feasible
- basic styling (stick model, VDW, ball and stick, wireframe) and sensible
  defaults
- support for animations (input as multiframe .xyz files) including dynamically
  drawing bonds
- find and fill in aromatic rings
- draw spheres sitting on atoms to represent atomic charges and dynamically
  scale them during an animation
- draw molecular orbital and density isosurfaces with .cube files or .molden files

## Installation
There are two basic ways to install Molecular Blender, depending on whether you
want to use the optional cython enhanced isosurface routines (highly recommended,
if you're drawing isosurfaces at all).

### Using setup.py (cython compatible)
The most robust way to install is to use the included `setup.py`.

    cd /path/to/molecular/blender/repo
    python3 setup.py build
    python3 setup.py link

This will build the cython extensions to the code and then create symlinks to all of
the versions of blender found in the normal blender install directories (only tested on macOS).
For this to work, you'll need to use python 3.7 (the same version used in Blender 2.80),
and have cython installed for it. The same setup can be used to export a build:

    python3 setup.py build
    python3 setup.py export

This will create a zip file that should be compatible with the add-on installation
machinery of blender itself.

### Direct symlink (no cython)
If you are okay with skipping sython, you can directly symlink from the python
package to the blender scripts directory. The non-cython version has all the same
functionality as the cython version, but may be significantly slower for certain
tasks (cython version is roughly 100 times faster at drawing orbitals from molden
files).

    ln -s /path/to/molecular/blender/repo/molecular_blender /path/to/blender/scripts/addons/

and then it should appear in the list of Import-Export addons that can be
activated like any other addon.

On Mac OS X, the path for a user supplied addon is

    /Users/<username>/Library/Application\ Support/Blender/<version>/scripts/addons

where `<username>` and `<version>` should be replaced with your username and the
Blender version you are using.

## Isosurfaces

To draw isosurfaces, you need either a cube file or a molden file. If using a molden file,
you also need to specify which orbital or density you want to draw. In the orbital selection
box, you can specify a comma separated list of orbitals to draw. molecular_blender also
understands the keywords homo, lumo, and density, and they can be combined with additional math.
For example:

  - `homo - 1` would compute the second highest occupied molecular orbital
  - `lumo + 2` would compute the third highest unoccupied molecular orbital
  - `4` would compute the 4th molecular orbital
  - `density` would compute the density

## Design

The only guiding design principle so far has been to try and separate the molecular
data structures (Molecule, Bond, Atom) from the Blender manipulations. This is
subject to change.
