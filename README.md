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

## Explicit double and triple bonds

Multiple bonds use the same extruded Bézier centerlines as single bonds, with a
bevel object containing two or three circular profiles. Each strand has 70% of
the single-bond radius, with a gap of 30% of that radius between strands. Stick
model junction spheres widen at multiply bonded atoms to cover the strand ends.
Single bonds retain their existing drawing behavior.
Bond orders are currently explicit; importing a molecule still only perceives
connectivity.

From Blender Python, an existing bond curve can be changed with:

```python
from molecular_blender.plotter import set_bond_order

set_bond_order(bond_object, 2, direction=(0, 1, 0))
```

The direction specifies the separation between strands in the curve's local
coordinates, projected perpendicular to the bond. For generated bonds these are
the molecule's coordinate axes. For a bond in the XY molecular plane, supply an
in-plane direction perpendicular to it, rather than the plane normal. A direction
parallel to the bond is rejected. Omit the direction for a deterministic default,
or use order `1` to restore the round profile. Apply this to both halves of a
split-color bond. Alternatively set `Bond.order` and `Bond.direction` before
calling `PlotMolecule`.

Endpoint hooks and split colors continue to work. The initial orientation is
encoded as curve-point tilt; it does not track a changing molecular plane during
animation. Wireframe rendering does not use these profiles.

Run the manual example from the repository root (this clears the current scene):

```bash
blender --background --factory-startup --python examples/multiple_bonds.py -- \
  --save /tmp/multiple_bonds.blend --render /tmp/multiple_bonds.png
```

## Installation

### For Blender Users (Recommended)

The easiest way to install Molecular Blender in Blender is using the pre-built addon zip:

1. **Download** `molecular_blender.zip` from the [GitHub Releases](https://github.com/smparker/molecular-blender/releases) page
2. **In Blender**: Edit → Preferences → Add-ons → Install
3. **Select** the downloaded zip file (do not extract it!)
4. **Enable** the "Molecular Blender" addon

The zip includes pre-compiled Cython extensions for all platforms (Linux, macOS, Windows). Blender automatically uses the correct ones for your system.

### For Python Developers (Advanced)

If you want to use molecular_blender as a Python package outside of Blender, or contribute to development:

**Install from PyPI** (once published):
```bash
pip install molecular_blender
```

**Or download wheels** from [GitHub Releases](https://github.com/smparker/molecular-blender/releases):
```bash
pip install molecular_blender-*.whl
```

### Building from Source

There are two basic ways to install Molecular Blender from source, depending on whether you
want to use the optional cython enhanced isosurface routines (highly recommended,
if you're drawing isosurfaces at all).

#### Using blender_setup.py (cython compatible)
The most robust way to install from source is to use the included `blender_setup.py`.

**To create a Blender-ready zip file:**
```bash
cd /path/to/molecular/blender/repo
python3 blender_setup.py build
python3 blender_setup.py export
```

This creates a zip file that can be installed directly in Blender (Edit → Preferences → Add-ons → Install).

**To install directly to Blender's addon directory:**
```bash
python3 blender_setup.py build
python3 blender_setup.py install  # or use --mode=link for a symlink
```

This will build the Cython extensions and install/link to all Blender versions found in the standard directories (tested on macOS).

**Requirements:**
- Python version matching your Blender version (3.9 for Blender 3.0-3.1, 3.10 for 3.2-3.6, 3.11 for 4.0-4.1, 3.12 for 4.2+)
- Cython installed: `pip install cython`
- NumPy installed: `pip install numpy`
- C++ compiler (gcc, clang, or MSVC)

#### Direct symlink (no cython - not recommended)
If you are okay with skipping Cython compilation, you can directly symlink from the python
package to the Blender scripts directory. **Note:** The non-cython version is roughly 100x slower
at drawing orbitals from molden files, so this is not recommended if you're using isosurfaces.

```bash
ln -s /path/to/molecular/blender/repo/molecular_blender /path/to/blender/scripts/addons/
```

On macOS, the user addon path is:
```
/Users/<username>/Library/Application Support/Blender/<version>/scripts/addons
```

On Linux:
```
~/.config/blender/<version>/scripts/addons
```

On Windows:
```
%APPDATA%\Blender Foundation\Blender\<version>\scripts\addons
```

After symlinking, activate the addon in Blender's preferences.

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

## Development and Contributing

### For Blender Addon Development

The standard development workflow for working on the addon:

```bash
python blender_setup.py build      # Build Cython extensions
python blender_setup.py install    # Install to Blender (or --mode=link for symlink)
python blender_setup.py export     # Export Blender-ready addon zip
```

This creates a zip file that can be distributed to Blender users or uploaded to GitHub releases.
