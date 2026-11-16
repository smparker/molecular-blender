# Molecular Blender [![Build](https://github.com/smparker/molecular-blender/actions/workflows/blender-plugin.yml/badge.svg?branch=master)](https://github.com/smparker/molecular-blender/actions/workflows/blender-plugin.yml)
An import plugin for molecular XYZ, SDF/MOL, Molden, cube, and JSON data used in
quantum chemistry.

The import dialog's **Import Details** panel reports the active orbital evaluator
and surface mesher separately as **Cython** or **Python**. Python fallbacks work
without a compiler. If an installed compiled module fails to load, the original
error is reported so an incompatible binary or missing dependency can be diagnosed.

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
SDF/MOL imports retain their explicit bond tables. For coordinate-only formats
such as XYZ, optional perception assigns bond orders from connectivity and valence.

### Perceiving bond orders from coordinates

Enable **Draw bond orders** under Import Details, or pass
`draw_bond_orders=True, perceive_bond_orders=True` to `BlendMoleculeSetup`. Supply `molecular_charge`
for ions (default: `0`). This is an integer total charge, not a sum of partial
charges from the input. The checkbox is off by default: bonds draw as single
strands, including SDF/MOL bonds. Enabling it draws explicit orders or perceives
them when no bond table is supplied. Stored SDF/MOL orders are never replaced.
The Python API exposes drawing and perception separately, both off by default;
use only `draw_bond_orders=True` when orders are already supplied.

The approach follows the connectivity/valence/unsaturation/charge separation in
[RDKit DetermineBonds](https://github.com/rdkit/rdkit/blob/master/Code/GraphMol/DetermineBonds/DetermineBonds.cpp)
and [xyz2mol](https://github.com/jensengroup/xyz2mol). The existing
`1.2 * (covalent_radius_i + covalent_radius_j)` distance cutoff supplies the
single-bond graph. A dependency-free constrained search then assigns orders 1–3
so every atom reaches an allowed valence and formal charges add up to the supplied
total. It prunes impossible assignments, tries shorter normalized distances first
among equally constrained edges, and backtracks if a guess blocks completion.
It minimizes total absolute formal charge before accepting charge-separated forms.

Unlike RDKit's maximum-matching implementation, this search checks alternative
edge orders directly. Allowed main-group states cover H, B, C, N, O, F, Si, P, S,
Cl, Ge, Br and I, including common B/C electron-deficient states and P/S expanded
valences. For example, it prefers neutral `S=O` for sulfoxides; RDKit's
DetermineBonds currently selects the charge-separated `S+–O−` representation.
An optional test compares both implementations across 24 molecular graphs.
No RDKit installation is needed to run the importer.

All atoms, including hydrogens, must be present in the file. Hiding hydrogens is
fine; they still participate in the assignment. Missing hydrogens, incorrect
connectivity, unsupported elements, or an inconsistent total charge can prevent
an assignment. Failures raise an error without applying tentative orders or
charges; the search also stops at `bond_order_max_iterations` nodes (default:
100,000). This is a limited closed-shell valence model, not full RDKit sanitization:
radicals, coordination chemistry, stereochemistry, and aromaticity classification
are not inferred. Resonance rings receive one valid integer Kekule assignment;
equivalent or ambiguous forms need not have a unique result. Bond lengths guide
the search rather than proving chemical identity.

For trajectories, only `animate_bonds="staticfirst"` is supported with perception.
Orders and initial strand orientation then stay fixed; they are not recalculated
for each frame. Dynamic and union connectivity are rejected for perception.

To render the same example from XYZ coordinates, with no bond table:

```bash
blender --background --factory-startup --python examples/sdf_bonds.py -- \
  --perceive --save /tmp/perceived_bonds.blend --render /tmp/perceived_bonds.png
```

### Importing bond orders from SDF/MOL

The native reader supports `.sdf`, `.sd`, and `.mol` files with V2000 connection
tables ([CTfile specification](https://www.daylight.com/meetings/mug05/Kappler/ctfile.pdf)).
It reads coordinates in angstroms and explicit single, double, and triple bonds.
The supplied bond table is authoritative, including when it is empty; distance
cutoffs do not add or remove bonds. Neighboring atoms in that graph set the
initial strand separation direction in the local molecular plane. Collinear or
isolated bonds use the deterministic default orientation.

Import through File > Import > Molecular Blender, or from Blender Python:

```python
from molecular_blender.plotter import BlendMoleculeSetup

molecule, styler, options = BlendMoleculeSetup(
    bpy.context, "examples/bond_orders.sdf",
    plot_style="fixedbs", bond_thickness=0.14, ignore_hydrogen=False,
    draw_bond_orders=True,
)
```

An SDF can contain multiple records. The importer loads the first by default;
pass `sdf_record=1` to load the second (zero-based indexing). Each selected record
is a static model, not an animation. Hydrogen hiding preserves indices and
connectivity. Formal charges are retained as `Atom.formal_charge`, separately
from partial charges used by charge visualization.

This first reader rejects V3000 and aromatic/query bond codes; aromatic molecules
can be supplied with explicit alternating single/double bonds. Stereo, isotope,
radical, S-group and SD data annotations are not interpreted. Implicit hydrogens
are not generated.

The example contains `O=CH-CH=CH-C#N` with all three bond orders. Render the same
file as sticks and ball-and-stick models with:

```bash
blender --background --factory-startup --python examples/sdf_bonds.py -- \
  --save /tmp/sdf_bonds.blend --render /tmp/sdf_bonds.png
```

### Setting bond orders in Python

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

Use `--link-target` with either install mode to select discovered Blender versions,
and `--dry-run` to preview destinations before replacing existing installations:

```bash
python3 blender_setup.py install --mode=link --link-target=4.3 --dry-run
python3 blender_setup.py install --mode=copy --link-target=4.3,4.10
```

`blender_setup.py` builds, installs, and exports the Blender addon. `setup.py` and
`pyproject.toml` define the Python wheel and source distributions used by packaging CI.

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

### Selecting cube datasets

Cube imports retain every dataset and its file label. **Cube dataset index** in
the import dialog chooses a zero-based dataset; Python callers can pass
`cube_dataset=1` to select the second. The default remains the first dataset.
`VolumeData.datasets` stores all fields, `labels` retains their identifiers,
and `select_dataset(index)` changes the active three-dimensional `data` array.
Comment lines are preserved without guessing physical field types or units.

Cumulative fractions must lie within `[-1, 1]`; signs choose surface signs.
Zero containment excludes all samples, full containment includes all nonzero
samples, and an all-zero field has no defined containment threshold. Finite
grids can only approximate intermediate fractions. The numerical helper's
`square=True` integrates orbital squared amplitude; `square=False` integrates
absolute density or signed-field magnitude.

Orbital requests use positive one-based indices or HOMO-relative offsets
(`0` means HOMO, `-1` means the preceding orbital); invalid indices raise an
error. Occupations must be finite and within `[0, 2]`. Missing spin labels allow
a visualization alpha default. For a single spatial-orbital set (maximum
occupation above 1), alpha fills first up to 1 and the remainder goes to beta;
A Beta-labelled spatial orbital reverses that majority convention. Thus doubly occupied
orbitals contribute equally to both spins. Fractional occupations are retained
in densities. These defaults are visualization heuristics, not a determination
of the wavefunction type.

### Density and frontier selection

The **Volume** selector offers total, alpha, beta, and alpha-minus-beta spin
density. Closed-shell zero spin produces no surface. Missing spin labels use
the visualization defaults above; no additional input is required.

For frontier names, the occupation scale is 1 when the maximum occupation is
at most 1, otherwise 2. Occupations strictly above half that scale are occupied;
exactly half and lower are empty. Actual density weights remain fractional.
HOMO is the highest-energy occupied orbital and LUMO the lowest-energy empty
orbital. If energy metadata is incomplete in a channel, selection uses input
order for that channel. `homo-N` and `lumo+N` traverse occupied/empty candidates
in that order. **Orbital spin** chooses alpha/beta, or use `alpha:homo` and
`beta:lumo+1`. Positive numeric selections always use the original one-based
file index regardless of spin selection.

**Cube field type** controls cumulative integration: orbital amplitude uses
squared values; density and signed scalar fields use absolute magnitude. Auto
uses field metadata if available and otherwise defaults to scalar magnitude,
so an untyped cube is not silently squared. Absolute-isovalue surfaces are
unaffected by this setting. Python callers use `cube_field_type='orbital'`,
`'density'`, `'signed'`, or `'auto'`.

Surface calculations started from the import dialog run in a separate background
Blender process. Escape stops that process, clears progress, and discards its
results; the already imported atoms and bonds remain available. Deleting the
molecule, opening another file, or disabling the addon also stops pending work.
The synchronous `BlendMolecule` Python API continues to finish before returning.

New imports store a versioned molecule record in the `.blend`: source path,
Angstrom coordinate units, import settings, atom identities and charges, bond
endpoints and order provenance, orbital metadata, and generated surface levels.
The source calculation is not embedded; the saved geometry and metadata remain
available when the source file is missing. Rename objects freely. To duplicate a
whole molecule, duplicate its root and children together; the copy receives a
new molecule identity while preserving atom identities within that molecule.
Older imports without this record must be reimported to use molecule editing.

In the 3D Viewport, open the sidebar (`N`) and choose **Molecule**. Select an
imported atom, bond, surface, or molecule root to inspect its source, counts,
coordinate units, charge setting, representation, and available scientific data.
Inspection uses the saved record and does not recalculate surfaces or reread the
source. The charge setting is the requested molecular charge for perception;
it is distinct from each atom's partial charge.

Choose **Restyle molecule…** in that sidebar to change atom size, bond thickness,
palette, and surface opacity. Apply acts on the whole owning molecule and supports
Blender Undo/Redo. Linked duplicates and recycled materials are isolated on edit;
multiple-bond profiles and stick junctions follow the new thickness. Surface
opacity wraps the existing shader without recomputing the field or mesh.
Individual atom/bond restyling is available for mesh, NURBS, and metaball imports;
wireframe imports remain inspectable. Use positive bond thickness when restyling.

For manual topology edits, select exactly two atoms in one molecule and click
**Single**, **Double**, or **Triple**. These buttons add a missing bond or update
an existing one. You can also select one bond (either colored half, or both) to
change its order or use **Delete bond**. Both halves and neighboring stick
junctions update together. New manual bonds follow their atoms through hooks.
All these actions support Undo/Redo and survive saving and reopening the file.

Manual orders are saved as locked assignments. **Perceive unlocked bonds** works
on the current graph, retaining manual additions/deletions and protecting both
manual orders and explicit source orders. Conflicting constraints report an
error without applying a partial assignment. This post-import perception action
requires a static molecule with all atoms, including hydrogens, displayed.
Manual atom moves change the displayed geometry; they do not recalculate orbital
or density surfaces. Bond perception itself remains an explicit, opt-in action.
