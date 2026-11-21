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
- toggle a molecule’s charge spheres together in its **-charges** collection in the
  Outliner. Each sphere scales independently during playback
- render XYZ gradients as arrows with **Plot gradient**. Supply rows in the form
  `element x y z charge gx gy gz`. Each arrow starts at its atom and points to
  `(x + gx, y + gy, z + gz)`. A zero gradient produces no visible arrow
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
Orders stay fixed. With endpoint hooks enabled, double and triple bonds follow
the molecular plane at each sampled pose. For a ring bond, the smallest ring
containing it defines the plane. For an open chain, a neighboring atom defines it.
When the atoms become collinear, the previous orientation follows the bond axis.

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

Endpoint hooks and split colors continue to work. Imported trajectories animate
the curve-point tilt on both halves of a multiple bond. The sampling and
interpolation settings control the timing of those keys.

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
- A build targeting Blender's bundled Python ABI. Check `sys.version` in Blender's Python Console. This project's local Blender 4.3.2 uses Python 3.11.9 and NumPy 1.24.3.
- Cython 3 for generating the native extension sources
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

The import dialog adapts to the highlighted filename. XYZ and SDF/MOL show
structure controls. Molden and cube files also show **Surfaces**. For JSON
trajectories the first referenced frame determines the surface controls.
The filename selects the reader, so give the file an extension matching its contents.
If building the molecular geometry fails, Blender removes the partial import and
restores your active object, selection and editing mode. You can retry in the same scene.

Cube imports retain every dataset and its file label. **Dataset** appears for
multi-dataset cubes and shows file IDs when available. Very large sets use a
zero-based index field. Python callers can pass `cube_dataset=1` to select the
second dataset. The default remains the first dataset.
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

For Molden files, choose **Molecular orbital** under **Field** to select HOMO,
LUMO, or a **Custom…** expression. The same field menu offers total, alpha,
beta, and alpha-minus-beta spin density. Choosing a density hides the orbital
controls. Closed-shell zero spin produces no surface. Missing spin labels use
the visualization defaults above.

For a JSON entry that pairs a Molden file with a Turbomole DIPL file through
`polar`, choose **Field → Response density**, then select an **Eigenpair** by
its DIPL ID and eigenvalue. Both signs of the induced density share an absolute
contour level. Try `examples/phenanth_polar.json`, which contains six eigenpairs.
For the Python setup API, use `volume="response", response_eigenpair=3` to
select eigenpair 3. Zero selects the first pair listed in the file.

For frontier names, the occupation scale is 1 when the maximum occupation is
at most 1, otherwise 2. Occupations strictly above half that scale are occupied;
exactly half and lower are empty. Actual density weights remain fractional.
HOMO is the highest-energy occupied orbital and LUMO the lowest-energy empty
orbital. If energy metadata is incomplete in a channel, selection uses input
order for that channel. `homo-N` and `lumo+N` traverse occupied/empty candidates
in that order. **Orbital Spin** appears for unqualified frontier expressions
and HOMO-relative indices. It chooses alpha/beta, or use `alpha:homo` and
`beta:lumo+1`. Positive numeric selections always use the original one-based
file index regardless of spin selection.

**Level mode** selects **Fraction enclosed** or **Absolute isovalue**. Clear the
levels to import just the structure. For cubes, **Advanced → Include both signs**
defaults on. Turn it off to plot only the signs entered.

Cube voxel counts set one length unit for the grid axes, origin, and atoms.
Positive counts use Bohr and negative counts use Å. The reader accepts rotated
orthogonal grids with unequal spacings. It reports mixed signs and sheared axes.

**Cube field type** appears for enclosed fractions and controls integration: orbital amplitude uses
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

If a Molden file moves, select its molecule and click **Locate Source…**
beside the missing path. Choose the relocated file. The add-on checks the saved
atom indices, elements, coordinates, and orbital metadata before updating the
path. You can then add or recompute surfaces in the same scene.

In the 3D Viewport, open the sidebar (`N`) and choose **Molecule**. Select an
imported atom, bond, surface, or molecule root to inspect its source, counts,
coordinate units, charge setting, representation, and available scientific data.
Inspection uses the saved record and does not recalculate surfaces or reread the
source. The charge setting is the requested molecular charge for perception;
it is distinct from each atom's partial charge.

Click **Frame Molecule** to fit the molecule and its visible surfaces in the
current viewport. It preserves the selection. In camera view it switches to a
perspective view and keeps the render camera in place. Select three noncollinear
atoms and click **View Face-On (3 Atoms)** to look perpendicular to their plane.
**View Principal Axes** uses the displayed atom positions at the current frame.
It looks along the molecule's shortest dimension and puts its longest dimension
across the screen. For a nearly circular planar molecule, the current view
sets the in-plane angle. These axes give each displayed atom equal weight.

Select an atom to see its element, stored atom index, partial charge, and formal
charge. The index follows the importer's numbering. Charge values come from the
saved atom record, which defaults missing charges to zero. Cube imports label
the charge from the atom line as **Cube atom charge**.

**Select Atoms by Index…** accepts a comma-separated list such as
`4, 7, 12-16`. It uses the indices shown in the sidebar, includes both ends
of a range, and selects visible, selectable atoms on the current molecule.
You can then use **Style Selected…** or, with three atoms, **View Face-On**.
Turn on **Show Atom Indices** to see these numbers beside the visible atoms of
the active molecule. The labels follow the atoms as you move the view or play
an animation. Blender's viewport Overlays switch controls their display.

**Select Same Element** selects the visible atoms of that element on this
molecule. **Select Bonded Neighbors** selects the active atom and its directly
bonded visible neighbors, including bonds added in the sidebar. **Select Bonded
Fragment** follows those bonds to select every visible, selectable atom in the
same connected fragment. A hidden atom can connect two visible atoms.

**Style Selected…** applies an accent color to the selected atoms and to the
bonds between them. You can also select either half of a bond to color the
whole bond. **Restore Element Colors** returns the selection to the current
molecule palette, including the usual neutral bonds for ball-and-stick models.

**Reset Atom and Bond Style** restores the selected molecule's imported atom
size multiplier, bond thickness and palette, and clears all fragment accents.
You can use it after editing bond orders or styling a surface. Older blend files
use the atom and bond settings saved before the first restyle or reset with this
version.

Expand **Orbitals** to browse the saved index, energy, occupation, spin, and
symmetry in input order. Use the **Page** arrows or enter a page number.
**Go to HOMO** and **Go to LUMO** jump to and mark the matching row for the
chosen **Orbital Spin**, even when the source file is unavailable.
A dash marks missing metadata. Energies are shown as supplied by the source.

For a single-frame Molden import, click **Add Surface…** next to an orbital,
or choose a density and click **Add Surface…** beside it. Adjust the levels
and grid spacing in the dialog. The calculation reads the original Molden file
and runs in the background. Escape cancels it. New surfaces follow the molecule
root's transforms and use the coordinates from the source calculation.

Under **Surfaces**, each row shows a field and its calculated isovalues. Click
the eye to show or hide all lobes and levels of that calculation in both the
viewport and render. Visibility is saved in the blend file and works offline.
The color swatches change the positive and negative lobes independently. The
opacity slider affects every lobe, level, and frame in that row. Both controls
edit the materials directly, so you can see the result in Material Preview or
Rendered shading as you adjust them. The trash button removes that calculation
and supports Blender Undo.

Click **Solo** to show just that calculation on its molecule. Click another
row's **Solo** to switch, or click the active **Solo** again to restore the
previous visibility of each lobe. Using a row's eye also exits Solo.

Click **Export Surface PNGs…** to render the existing surface groups using the
active scene camera at the current frame. Choose a destination folder and check
the surfaces to export in the file browser's options sidebar (press **N** to
show it). Each PNG has a transparent background and shows one chosen surface
with the rest of your scene. The exporter restores visibility and render
settings afterward and adds a numeric suffix if a filename already exists.
Enable **Include field legends** to save an SVG beside each mapped PNG.
The pair shares a filename stem, including when an existing PNG or SVG
requires a new suffix. Ordinary orbital panels produce just a PNG.
Enable **Fit series camera from view** to use the current viewport orientation
for all checked panels. The exporter frames the molecule and the selected
surfaces together at the current frame, with a margin for the output aspect
ratio. It makes a reusable orthographic camera active so the framing is saved
with the scene. You can return to your original camera in Scene Properties.
For an ordinary still, **Output Properties → Output → Use Transparent PNG**
sets transparent film and RGBA PNG output for the scene.

To regenerate a series from a saved figure, run Blender with the add-on enabled.
Use the molecule root's Outliner name when the scene contains several molecules.
List the active camera, molecule names, surface labels and group IDs from the
saved scene with the first command. Pass a label or group ID to `--surface`
to render a particular panel.

```sh
blender -b "figure.blend" -c molecular_render --list
blender -b "figure.blend" -c molecular_render --molecule "water" --output "panels"
```

Add `--legends` to the render command to include the SVG field legends.

The command creates the output directory and renders one PNG per saved surface
group with the scene's camera, lighting, resolution and current frame. To render
a subset, repeat `--surface` with exact sidebar labels such as
`--surface "Orbital HOMO" --surface "Orbital LUMO"`. A failed render or unknown
name returns a nonzero exit code. Run `blender -c molecular_render --help` for
the arguments. Blender 4.3 supports this command.

**Make Molecule Independent** detaches geometry and materials shared with
objects outside the molecule, including multiple-bond bevel shapes. It keeps
internal sharing, so identical atoms can still reuse a mesh. The action
preserves placement, selection, and surface definitions.

For a single-frame Molden field, the refresh button (**Recompute Surface**)
changes its levels or grid spacing using the original file. The current
surface stays visible or hidden as you set it, including during calculation.
The replacement keeps its lobe colors and opacity. For surfaces saved by older
versions, the dialog starts from their stored absolute isovalues.
On a field-mapped density it resamples the original scalar field on the new
vertices and keeps the edited ramp and limits. If the grid is missing or too
small for the requested contour, the previous surface and map stay in place.

On a single-frame cube surface, **Recompute Surface** changes the contour
levels on the original dataset and grid. **Include both signs** mirrors each
level. Turn it off to enter signed levels individually. For a field-mapped
cube, Blender resamples the map on the new mesh and keeps the edited color ramp.

Choose **Restyle molecule…** in that sidebar to change atom size, bond thickness,
palette, and surface opacity. Apply acts on the whole owning molecule and supports
Blender Undo/Redo. Linked duplicates and recycled materials are isolated on edit;
multiple-bond profiles and stick junctions follow the new thickness. Surface
opacity wraps the existing shader without recomputing the field or mesh.
Individual atom/bond restyling is available for mesh, NURBS, and metaball imports;
wireframe imports remain inspectable. Use positive bond thickness when restyling.

On single-frame molecules with individually selectable atoms, choose **Change Representation…**
below the current representation in the Molecule sidebar. Switching to sticks
splits a bond between different elements into two colored halves. Switching to
ball-and-stick restores element-sized spheres and single-color bonds. The action
uses the displayed atoms and saved bond graph, including manual orders and tilt.
Existing orbital and density surfaces keep their placement and materials.
Choose **Space-filling (van der Waals)** to see the molecular envelope. Returning
to the previous representation restores each bond's viewport and render visibility.

For manual topology edits, select exactly two atoms in one molecule and click
**Single**, **Double**, or **Triple**. These buttons add a missing bond or update
an existing one. You can also select one bond (either colored half, or both) to
change its order or use **Delete bond**. Both halves and neighboring stick
junctions update together. New manual bonds follow their atoms through hooks.
All these actions support Undo/Redo and survive saving and reopening the file.

Select a displayed double or triple bond, or its two atoms, and choose
**Set strand tilt…** to rotate the strands around the bond axis. The dialog
sets the angle on both colored halves. On a trajectory with animated strand
orientation, the angle is an offset from the moving plane. Restyling preserves it.

Manual orders are saved as locked assignments. **Perceive unlocked bonds** works
on the current graph, retaining manual additions/deletions and protecting both
manual orders and explicit source orders. Conflicting constraints report an
error without applying a partial assignment. This post-import perception action
requires a static molecule with all atoms, including hydrogens, displayed.
Manual atom moves change the displayed geometry; they do not recalculate orbital
or density surfaces. Bond perception itself remains an explicit, opt-in action.

### Figure and animation examples

To color an atom selection and render the water HOMO and LUMO with the startup
camera, run the included example from the repository root. The output folder
will contain one transparent PNG per orbital.

```sh
blender --background --factory-startup --python-exit-code 1 --python examples/edit_orbitals.py -- /tmp/water-orbitals
```

To try a trajectory, import `examples/water_bend.json` with **Ignore Hydrogens**
off and scrub frames 1, 5 and 9. The manifest selects two illustrative poses
from a multi-record SDF and returns to the first pose.

In a multi-entry JSON manifest, each entry reads one structure. Paths are
relative to the manifest and `sdf_record` selects a zero-based SDF record.
Optional `frame` values are strictly increasing nonnegative integers. The
first defaults to zero and each omitted value advances one step. Blender
places the key at `1 + frame × keystride`. Atom identities and explicit SDF
topology must match across frames. Molden frames retain their own coefficients
and occupations and share a compatible basis and orbital channels. Cube frames
can have different grid origins and spacings and must use the same dataset IDs.

### Combine cube fields in JSON

A manifest entry can supply a `linear_combination` of cube datasets. For example,
this plots the polarized density minus the reference density.

```json
{
  "molecules": [{
    "units": "e/bohr^3",
    "resample": true,
    "linear_combination": [
      {"filename": "polarized.cube", "coefficient": 1},
      {"filename": "reference.cube", "coefficient": -1}
    ]
  }]
}
```

Add more terms or use fractional coefficients for weighted sums and averages.
Each term accepts a zero-based `cube_dataset` (default 0). The first cube
supplies the displayed atoms and output grid. By default, all grids must match
after coordinate conversion to angstroms. With `"resample": true`, the importer
uses trilinear interpolation for terms with a different spacing, origin, or
axis orientation. Each cube must cover the entire first grid. Supply all fields
in the same molecular coordinate frame before import. The `units` string declares the
scalar units shared by every term. A term can also declare its units, which
must match. Supply values in those units before import. Cube comments are
free text. **Include both signs** draws accumulation and depletion at equal
absolute levels. A combination can also be one frame of a JSON trajectory.

### Preview a Molden vibrational mode

Import a Molden file with normal coordinates, select one of its atoms, and open
**N sidebar → Molecule → Vibrational modes**. Choose a frequency in cm⁻¹, set
**Peak displacement (Å)**, and click **Preview Mode**. Press **Play** or scrub
the timeline. **Frames per cycle** controls the viewing speed. **Starting
phase** sets the position within the cycle at frame 1. Click **Update Preview**
after changing these settings or selecting another mode. **Restore Geometry**
removes the oscillation and restores the previous atom offsets.

Molecular Blender reads files containing only the vibration sections as well
as files with orbitals. It converts the [Molden equilibrium coordinates and
displacements](https://www.theochem.ru.nl/molden/molden_format.html) from Bohr.
Peak displacement sets the largest atomic motion in the mode and preserves the
relative lengths and directions of the other vectors. The preview and its mode
data are saved with the scene. For example, import `examples/zine.molden` and
select mode 7 (14.55 cm⁻¹).

The preview uses the atom objects of a single geometry. Mesh, NURBS and metaball
atoms support it. Leave **Hook bonds** on **Automatic** to let bonds follow the
motion. Turn off **Ignore Hydrogens** at import to see their displacements.
Electronic surfaces describe the source geometry and stay fixed during playback.

### Sample trajectory frames

For XYZ and JSON imports, **First input frame** and **Last input frame** are
inclusive, zero-based positions. A last frame of `-1` reads through the end.
**Read every Nth frame** counts from the first selected position. For example,
`100` through `500` with a step of `10` loads 41 frames. Explicit JSON frame
times survive sampling, and **Keystride** still controls the spacing of Blender
keyframes. The Python options are `frame_start`, `frame_end`, and `frame_step`.
A single JSON entry wrapping an XYZ trajectory samples the XYZ frames.

Under **Animation → Playback**, **Hold** keeps each sampled pose and charge
size until the next key. **Linear** moves atoms and changes charge sizes at a
constant rate between keys. **Blender Default** uses your keyframe preference.
The Python option is `trajectory_interpolation="hold"`, `"linear"`, or `"blender"`.
Bond and surface visibility switch on snapshot frames. Gradient arrows follow the
sampled positions and vectors.

Cube and Molden trajectory surfaces are built one snapshot at a time. Press
**Esc** during calculation to discard the partial surface series and keep the
imported molecular geometry. The finished meshes are stored in the blend file
for scrubbing and rendering.

For a JSON sequence of Molden files, choose **Molecular orbital** and enable
**Align orbital phase across frames** to keep lobe colors consistent. In Python,
use `align_orbital_phase=True`. For each selected orbital, the importer estimates
its signed spatial overlap with the preceding sampled frame in the supplied
coordinate system. A normalized overlap of −0.5 or lower reverses the sign for
display. If the overlap magnitude is below 0.5 or the spin channel changes,
the importer starts from the file's sign and flags that Blender frame in the
Molecule sidebar. Orbital indices and frontier selections still resolve separately
on each frame. An orbital crossing can require a different selection.

With dynamic inferred bonds, **Break margin (%)** keeps an existing bond
visible until its length exceeds the formation cutoff by the chosen percentage.
The default is 5%. Set `bond_break_margin=0` for a single cutoff. The rule is
evaluated on the sampled snapshots and applies only to distance-inferred bonds.

### Color a density surface with a cube or a JSON field combination

Import a density cube or add a **Total**, **Alpha**, or **Beta Density** surface
to a Molden molecule. In **Molecule → Surfaces**, choose **Map Field…** on
the single-frame density surface. Select the cube, dataset and color limits.
Use **Electrostatic Potential** for an ESP cube in Hartree/e or volts. For other
quantities, enter a field name and units (`1` for dimensionless values).
**Signed Field** and ESP use red for negative values, white at zero and blue
for positive values. **Unsigned Field** runs from white at the lower limit to
blue at the upper limit. Use the same units and limits across figures.
Use Material Preview or Rendered shading to see the colors.

To map induced density, select a JSON file containing one `linear_combination`
entry, such as the [polarized-minus-reference example](#combine-cube-fields-in-json).
Enter a field name and color limits spanning zero. The map uses the term
datasets, weights, resampling setting and units from the manifest. Its first
cube supplies the atom coordinates and sampling grid. **Export Legend…** uses
the field name and manifest units.

The field cube must list the same atoms in the same order and source
positions (within 0.0001 Å). Its grid can have different spacing, orientation,
or origin, and must cover the whole density surface. The add-on interpolates
the supplied values onto the displayed vertices after remeshing and saves
them as a point attribute (`mb_esp` for ESP or `mb_field` for a custom field).
The material's Color Ramp is editable in the Shader Editor.
Use **Edit Color Scale…** on a mapped surface to adjust its limits, colors
and ramp stops. To match another panel, choose it under **Copy scale from**
and confirm. The picker includes maps in the same scene with matching units
and signed or unsigned ranges. The edit uses the saved mesh values and works
after the source files move. A stop at zero follows zero when the limits change.
**Change Field Map…** reopens the file, field and range settings.
**Remove Map** restores the original geometry, material, and Remesh settings.

**Export Legend…** saves an SVG for compositing with an editable field name,
units and numeric labels. It reads the connected material's Color Ramp and
Map Range, including edited colors and limits. A signed scale marks zero at
its position between the two limits. The SVG uses sRGB swatches for the
material colors before scene lighting and the render view transform.

### Render a water ESP figure

Run this command from the checkout to create a transparent PNG and a matching
SVG legend. Blender reads the supplied Molden file and ESP cube directly.

```sh
blender -b --factory-startup --python-exit-code 1 --python examples/render_water_esp.py -- water-figure
```

The example uses an RHF water calculation from PySCF. It plots total density
at 0.001 e/Bohr³ and colors it by electrostatic potential from −0.035 to
+0.035 Hartree/e. Red marks negative potential. Blue marks positive potential.
To use your own calculation, change the two input paths, color limits and camera
framing in the [figure script](examples/render_water_esp.py).

### Molden normalization by writer

New imports use **Auto** under **Surfaces → Advanced → Molden convention**.
A producer title from `orca_2mkl` selects **ORCA**. Other files use **Standard
Molden**. The Molecule sidebar shows the convention used for the import. For
files with missing or edited titles, select either convention explicitly. In
Python use `molden_convention="auto"`, `"standard"`, or `"orca"`. Existing
scripts and saved scenes using `orca_fix=True` or `False` retain their explicit
ORCA or standard choice.

The ORCA convention uses the normalization already present in the contraction
coefficients, including the spherical d-shell factor. For spherical f and g
shells it also converts ORCA’s radial factors and the phases of components with
|m| ≥ 3. The [ethene calculation](test/fixtures/orca_ethene/ethene.inp) checks a
π orbital on both lobes and near the node. The [ORCA 6.1.1 water
calculation](test/fixtures/orca_water_fg/README.md) checks signed amplitudes for
all 23 orbitals.

The [PySCF 2.11 water references](test/fixtures/pyscf_water/generate.py) check
standard Molden with contracted d, f and g functions on oxygen. They cover all
28 spherical and 38 Cartesian orbitals, plus total density, at 64 off-axis
points in each calculation. The generator also recreates the ESP figure inputs.
