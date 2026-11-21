# Molecular Blender roadmap

An editable menu of possible work, based on the 2026-11-15 code review. Status was checked against the working tree on 2026-11-21. Checked items are implemented. Open items describe remaining work.

Keep the current curve-extrusion approach for bonds as the baseline. Keep bond-order drawing and perception opt-in until their behavior has broader validation.

## Suggested starting points

- **Extend post-import editing.** The sidebar already supports restyling, bond edits and surface management. Reloading coordinates and orbitals from a changed source is open.
- **Extend field maps.** External scalar cubes can color density surfaces. The sidebar exports field legends for compositing. Direct wavefunction ESP evaluation remains open.
- **Extend bond perception.** Diagnostics and manual overrides are available. Solver scaling and per-frame bond-order assignment remain open.

## 1. Post-import editing and scene interaction

Imports retain their settings and object identities in the blend file. The sidebar resolves the owning molecule through its root. See [scene records](molecular_blender/scene_state.py), [sidebar controls](molecular_blender/scene_panel.py) and [surface actions](molecular_blender/scene_surfaces.py).

- [x] **Persist molecule identity and settings in the blend file.** Versioned records store the source, settings, atom IDs, bond endpoints and order provenance, orbital metadata, and surface definitions. IDs survive renames and save/reopen. Duplicating a hierarchy gives its root an independent identity. [Scene record tests](test/test_scene_state.py)
- [x] **Add a selected-molecule sidebar panel.** Selecting the root, an atom, a bond, or a surface shows its source, counts, units, charge setting, representation, and available fields. [Panel tests](test/test_scene_panel.py)
- [x] **Restyle a molecule and individual surfaces.** Adjust atom size, bond thickness, palette, and surface opacity. Multiple-bond profiles and junctions resize together. Each surface calculation has independent lobe colors and opacity, including for linked duplicates. [Restyle tests](test/test_scene_edit.py) and [surface tests](test/test_scene_surfaces.py)
- [x] **Style selected atoms and bonds.** Apply an accent to a fragment or individual bonds and restore the current element palette.
- [x] **Add an orbital and surface browser.** Browse saved orbital metadata, jump to HOMO/LUMO by spin, and show, hide, style or remove a surface calculation. Single-frame Molden imports can add orbitals and densities after import. [Sidebar action tests](test/test_scene_surfaces.py)
- [x] **Replace a surface calculation in place.** Single-frame Molden surfaces can change levels and grid spacing. Cube contours reuse the selected dataset and grid. Both preserve lobe styling and placement. Mapped cubes resample their field colors.
- [x] **Recontour a mapped Molden density.** Changing levels or grid spacing resamples the supplied field on the replacement mesh with the edited ramp and limits.
- [x] **Browse additional fields on an imported molecule.** Choose a DIPL eigenpair by ID and value or a dataset from a multi-field cube. Each saved group can be styled, recontoured, soloed and exported independently.
- [x] **Expose bond editing tools.** Set an order, add/delete a bond, lock a manual assignment against reperception, and adjust strand tilt. Split bonds update both halves. Order changes update junction sizes, and profile rebuilds preserve manual tilt.
- [x] **Inspect and select atoms.** Inspect element, index, and charges. Select visible atoms by element, direct neighbor, or connected fragment using the current bond graph. [Panel tests](test/test_scene_panel.py) and [selection tests](test/test_scene_surfaces.py)
- [x] **Select atoms by index.** Select individual saved indices and inclusive ranges on the active molecule.
- [x] **Show optional atom labels.** The active molecule can display its saved atom indices as viewport overlays.
- [x] **Solo a surface.** Switch between calculations and restore the molecule's previous viewport and render visibility.
- [x] **Relink a moved Molden source.** Check saved atom identities, coordinates, and orbital metadata before updating the source path.
- [ ] **Refresh source data while preserving presentation.** Reload coordinates or orbitals and retain user colors, visibility, transforms, and locked edits where atom identities still match. Report incompatible atom lists rather than silently rebinding objects.
- [x] **Support undo and redo for sidebar edits.** Molecule restyling, manual bond edits, and surface removal support Blender Undo/Redo. Moving atoms changes the displayed geometry. Wavefunction surfaces use the original calculation coordinates.
- [x] **Reset atom and bond style.** The sidebar restores the original atom size, bond thickness and palette, and clears accents for the selected molecule. Older scene records use their saved settings as the baseline.
- [x] **Switch sticks, ball-and-stick and space-filling after import.** The sidebar converts single-frame molecules using their current geometry, preserving atom accents, manual bond edits, and existing surfaces.
- [ ] **Add a general Rebuild Geometry action.** Scope cleanup to the owning molecule. The sidebar already covers single-frame representation changes.
- [x] **Make a molecule independent.** Detach shared geometry, multiple-bond bevel sources and materials while retaining sharing within the molecule.

## 2. Electrostatic-potential maps

External potential cubes can color single-frame density surfaces. [Sampling](molecular_blender/esp.py) and [Blender mapping](molecular_blender/scene_esp.py).

- [x] **Map supplied ESP grids.** Sample a selected cube dataset on the final displayed density mesh, with atom alignment and coverage checks. Save per-vertex values, declared potential units, editable diverging colors and explicit limits. Removing the map restores the original surface.
- [x] **Generalize field-mapped materials.** Single-frame density surfaces support signed and unsigned scalar cubes with named fields, explicit units and editable color ramps.
- [x] **Map JSON field combinations.** Map a single weighted combination of cube datasets directly onto a density mesh using the manifest’s units and a signed color ramp.
- [x] **Export a field legend for compositing.** Save an SVG with the mapped material’s field name, units, edited Color Ramp, limits and signed zero tick.
- [x] **Edit and reuse field scales.** The mapped surface dialog edits the ramp and limits and copies compatible scales across molecules using the saved scene.
- [x] **Export legends with a surface series.** The sidebar and headless renderer can write a matching SVG beside each mapped PNG.

- [ ] **Implement a wavefunction ESP evaluator.** In atomic units, evaluate `V(r) = Σ_A Z_A / |r-R_A| - ∫ ρ(r') / |r-r'| dr'`. Reuse the AO density matrix, but add the electron–point Coulomb integrals; existing pointwise density evaluation does not supply this integral. Compare a native Gaussian-integral backend with an optional external calculation path. [PySCF's MEP example](https://github.com/pyscf/pyscf/blob/master/examples/1-advanced/031-MEP.py) is a useful implementation/reference calculation.
- [ ] **Define the scientific input contract before claiming general Molden ESP support.** Check occupations, orbital completeness, normalization, nuclear charges, ghost centers, and effective-core-potential treatment. Preserve needed atom metadata; require missing core-charge information or explicitly restrict the first implementation to validated all-electron inputs. See the [Molden format](https://www.theochem.ru.nl/molden/molden_format.html).
- [ ] **Keep exact and approximate electrostatics distinct.** A quick map from supplied partial charges could be useful, but label it as a point-charge approximation and retain its charge source. Do not substitute it silently for a wavefunction-derived ESP.
- [ ] **Validate ESP numerically.** Compare selected points and surface values for neutral and ionic reference molecules against an independent code; test the nuclear/electronic sign convention, unit conversion, distance from nuclear singularities, and far-field total-charge behavior.
- [ ] **Make ESP computation practical.** Evaluate surface points in cancellable batches, screen negligible shell pairs where justified, cache the density matrix and evaluated values, and allow coarse preview followed by a final calculation. Benchmark both time and peak memory before optimizing kernels.

## 3. Other scientific visualization features

Entry points: [orbitals.py](molecular_blender/orbitals.py), [compute_ao_py.py](molecular_blender/compute_ao_py.py), [compute_ao_cy.pyx](molecular_blender/compute_ao_cy.pyx), [isosurfaces.py](molecular_blender/isosurfaces.py), and [stylers.py](molecular_blender/stylers.py).

- [x] **Expose existing density modes properly.** The field menu offers total, alpha, beta, and spin density. Restricted, unrestricted and fractional populations follow the documented visualization defaults. Commit `ba0d216`.
- [x] **Show orbital metadata and support spin-specific frontier orbitals.** Retain energies and symmetry labels in `MOData`; select HOMO/LUMO separately for alpha/beta sets and handle fractional occupations and unsorted orbitals deliberately. Completed with commit `0a287ce` and visualization defaults approved 2026-11-16.
- [x] **Combine fields on matching grids.** JSON manifest entries support weighted sums of cube datasets with declared common scalar units. Signed surfaces share an absolute threshold.
- [x] **Resample field combinations across different grids.** JSON entries can request trilinear resampling onto the first term's grid. Each term must cover that grid in the supplied molecular coordinate frame.
- [x] **Import and animate normal modes.** The Molecule sidebar previews a saved Molden mode by frequency, peak displacement, phase and cycle length. Restore Geometry removes the preview. [Playback tests](test/test_vibrations.py)
- [x] **Frame the visible molecule in the viewport.** Frame Molecule fits the current geometry and visible surfaces while preserving the active selection and render camera.
- [x] **View a selected plane face-on.** Three atoms define the viewing direction for the framed molecule.
- [x] **Export a surface series.** Render selected surface groups to individual transparent PNGs with the active scene camera.
- [x] **Fit the camera to a surface series.** Export selected groups at a common scale from the chosen viewport orientation.
- [x] **Choose a Turbomole response density.** A Molden manifest with a DIPL file offers an eigenpair picker for rendering signed induced densities.
- [x] **Set ordinary renders to transparent PNG.** Output Properties includes a one-click film and RGBA output preset.
- [x] **Align the viewport to molecular principal axes.** Frame the current geometry with its longest dimension horizontal and best-fit plane face-on, using the current view to resolve ambiguous axes.
- [ ] **Improve publication output.** Standalone field legends are available for compositing. Scale bars remain a compositing task.

## 4. Bond perception and chemical representation

Current implementation includes covalent-radius connectivity, bounded valence/charge search, explicit SDF precedence, and animated strand orientation. See [bond perception](molecular_blender/bond_perception.py) and [strand orientation](molecular_blender/bond_orientation.py).

- [x] **Make failures explainable.** Report the offending atom IDs/elements, degrees, candidate valences, requested charge, and search budget. Distinguish bad connectivity, unsupported chemistry, missing hydrogens, incompatible charge, and exhausted search. Commit `a4c6961`.
- [x] **Allow manually specifying bonds and bond orders after import.** Select two atoms to create or change a single, double, or triple bond. Both halves of a split-color bond update together. [Bond edit tests](test/test_bond_edit.py)
- [ ] **Scale the solver and expose ambiguity.** Partition independent graph work where possible, propagate charge constraints across fragments, add timing/node diagnostics, and benchmark symmetric/fused systems. Preserve search limits and atomic application of successful results.
- [x] **Keep multiple-bond orientation coherent across frames.** Hooked trajectory bonds follow the smallest ring containing the bond or a local neighbor. Degenerate frames carry the preceding direction forward. Sidebar tilt offsets apply to both halves across the animation.

## 5. Import, export, and data fidelity

Entry points: [importers.py](molecular_blender/importers.py), [transform.py](molecular_blender/transform.py), and [turbomole_reader.py](molecular_blender/turbomole_reader.py).

- [x] **Harden Molden basis parsing.** Expands `sp` shells, registers standalone `[7F]`, maps GTO blocks by atom reference, reads Fortran D/d numbers, and validates indexed coefficient blocks and section endings. Existing molecular fixtures and focused synthetic cases pass. Commit `f1547e0`.
- [x] **Detect the Molden normalization convention.** Auto recognizes the ORCA producer title. Import and sidebar controls support explicit standard and ORCA choices, including settings from older scenes and scripts.
- [x] **Validate Cartesian d, f and g orbitals.** The [PySCF fixture](test/fixtures/pyscf_water/generate.py) checks 38 signed Cartesian orbitals and total density by point and by grid.
- [x] **Validate spherical d orbitals against ORCA.** An RHF/def2-SVP ethene fixture checks both lobes and the nodal plane of the occupied π orbital. Corrected the normalized d-harmonic expansion and ORCA spherical d contraction factor. [Reference test](test/test_transform.py)
- [x] **Validate spherical f and g orbitals.** The coefficient tables reproduce all 16 real harmonics in Molden order. ORCA imports handle its radial and phase conventions. Point and grid tests include weak tails from an independent [ORCA water calculation](test/fixtures/orca_water_fg/README.md).
- [x] **Validate spherical Molden import from a second writer.** A PySCF RHF water calculation supplies independent values for all 28 signed orbitals and total density on a 64-point grid with the ORCA fix off. The density check exposed premature screening of a contracted oxygen s shell. [Reference generator](test/fixtures/pyscf_water/generate.py)
- [x] **Preserve all cube datasets.** All cube fields and dataset IDs are retained, with a zero-based importer/API selector and first-field compatibility. Optional header value counts and malformed lengths are tested. Free-form comments are retained without guessing physical field types or units.
- [x] **Audit cube coordinates and grid conventions.** Tests cover both voxel-count unit conventions, nonzero origins, atom positions on a known surface, rotated orthogonal grids with unequal spacings, and surface normals on left-handed grids. Dataset tests cover optional header values and malformed data lengths. The reader reports mixed signs, zero counts, and degenerate or sheared axes. General sheared-grid support remains a possible extension.
- [ ] **Consider PDB/mmCIF, MOL2, and periodic inputs.** Choose formats based on actual use cases: biomolecular hierarchy, typed bonds, or unit cells. The disabled Babel fallback has been removed; unsupported formats now report the supported formats.
- [x] **Check atom counts and element order across JSON trajectory frames.** The reader resolves paths relative to the manifest and checks each subsequent frame against the first. [JSON reader](molecular_blender/importers.py)
- [x] **Validate JSON trajectories.** Check frame indices, references, atom identities, SDF topology, Molden basis and orbital channels, and cube dataset IDs. Each cube frame uses its own grid. Atom and surface keyframes honor manifest timing and SDF entries can select records.

## 6. Trajectories and interactive performance

Entry points: `AnimateMolecule`, `compute_surface_vertex_sets`, and `create_surface_meshes` in [plotter.py](molecular_blender/plotter.py).

- [x] **Expose trajectory sampling.** XYZ and JSON imports accept inclusive input ranges and a sampling stride. Only selected snapshots are allocated and original key times are preserved.
- [x] **Choose trajectory interpolation.** Hold or linearly interpolate sampled positions and charge sizes, or use Blender's keyframe preference.
- [x] **Reduce surface visibility keyframes.** Each lobe records its initial viewport and render visibility, then keys only transitions, including transitions through empty frames.
- [ ] **Expose physical playback timing.** Add time units and frame rate options.
- [x] **Stabilize dynamic inferred bonds.** Separate formation and breaking distances with an adjustable break margin.
- [ ] **Design per-frame bond-order perception.** Handle strand appearance/disappearance and charge conservation.
- [x] **Align orbital phases across frames.** Molden trajectories can use sampled real-space overlap to stabilize lobe colors. The sidebar identifies frames with ambiguous phase.
- [x] **Stream trajectory surface construction.** Cube and Molden trajectories build Blender meshes a frame at a time. The interactive worker queues at most two completed frames on disk. Escape, worker failure and add-on shutdown discard staged surfaces. A 24-frame signed cube benchmark produced the same 48 meshes and 1,498,176 faces, with peak memory falling from 333 to 307 MiB in Blender and from 341 to 302 MiB in the worker.
- [ ] **Load trajectory meshes on demand during playback.** Completed scenes currently retain all final meshes for scrubbing, rendering and saving.

## 7. Architecture and data ownership

The existing split between field computation and main-thread mesh creation is a useful starting point; extend it instead of discarding it.

- [x] **Validate and persist import settings.** The importer, Python setup API and scene records use the same validation function. Runtime options are derived from a separate copy. [Settings model](molecular_blender/options.py) and [settings tests](test/test_scene_state.py)
- [ ] **Unify option declarations across entry points.** Share UI and API property definitions while preserving their intentional default differences.
- [ ] **Separate scientific data from Blender objects.** Keep parsing, topology, perception, scalar fields, and numerical calculations usable without `bpy`/`mathutils`; let a scene adapter own object/material names and handles. Avoid a missing Blender-only import being indistinguishable from a broken package import.
- [ ] **Split `plotter.py` by responsibility.** Extract atom/bond geometry, surface generation, animation, scene ownership, and driver orchestration behind small interfaces. Factor the duplicated split/unsplit curve setup while preserving the tested extrusion and hook behavior.
- [ ] **Introduce a scalar-field interface.** Standardize point/batch/grid evaluation, bounds, units, sign semantics, provenance, and optional derivatives for cube fields, orbitals, densities, and future ESP. Keep numerical values independent of color/material choices.
- [x] **Locate malformed Molden basis and orbital data.** The Blender import report gives the source path, physical line number and GTO atom or MO orbital for bad or incomplete blocks.
- [ ] **Use structured errors and diagnostics across the remaining readers.** Carry filename/record/section/atom context from readers and chemistry into Blender reports; retain detailed tracebacks in logs. Replace broad exception swallowing and misleading fallback messages with precise causes.
- [x] **Locate malformed XYZ and cube input.** Report physical lines and atom, frame, axis or field context. Truncated cubes report the missing value count. JSON errors retain the failing term and cube location.
- [x] **Locate malformed SDF and Turbomole response data.** Report the physical line and selected SDF record, atom, bond, or charge block. Polarization imports identify the DIPL section and eigenpair and retain their JSON manifest context.
- [ ] **Make object/material reuse explicit.** Cache by molecule/style identity and validated parameters, not just global names. Provide safe cleanup for molecule-owned meshes, materials, actions, bevel profiles, and collections, including fake-user profiles.
- [x] **Repair repeat imports of rings and charges.** Recycled ring materials are bound on each molecule and charge spheres use the materials returned for that import.
- [x] **Document the orbital evaluator interface.** The public module documents Bohr and grid axis order and exposes `OrbitalCalculator` and `DensityCalculator` alongside the original spellings.
- [ ] **Clarify remaining numerical APIs.** Add type hints and shape/unit contracts at other module boundaries and distinguish atom IDs from positional indices.

## 8. Correctness and lifecycle fixes to investigate first

These are concrete targets spotted in the current code, not a claim that every path has been reproduced in a running Blender session.

- [x] **Fix element-dependent bond-radius lookup.** Uses the selected bond-style element consistently; mixed-element Blender regressions cover radius selection and unaffected styles. Commit `dbbafb3`.
- [x] **Audit isosurface method forwarding and resolution units.** Forwards the requested method, converts both adaptive spacings from angstrom to Bohr, and rejects invalid methods/spacings. Analytic Gaussian-sphere regressions check displayed geometry and final sampling spacing with both backends. Commit `537bf83`.
- [x] **Harden cumulative-isovalue calculations.** Validate fractions, signs, all-zero fields, non-finite inputs, tiny grids, and termination tolerances. Define whether a field integrates `ψ²`, density, or absolute signed density; an untyped cube should not silently inherit orbital-square semantics. Completed with commit `500a232` and visualization defaults approved 2026-11-16.
- [x] **Audit occupations and spin densities.** Validate absent/inconsistent spin labels, fractional occupations, electron counts, and alpha/beta frontier indexing. Check the Cython density constructor's integer electron-count argument against fractional-occupation inputs.
- [x] **Implement real cancellation and worker cleanup.** Escape, a deleted owner, loading another file, and disabling the add-on cancel pending work. Workers are reaped, temporary files removed, and cancelled results discarded. Progress includes orbital trajectory frames. [Worker tests](test/test_surface_jobs.py)
- [x] **Move surface calculations into a cancellable subprocess.** Import and sidebar actions launch a separate Blender process. The main thread handles progress and creates meshes. Commit `f306fb7`. [Worker implementation](molecular_blender/surface_jobs.py)
- [x] **Benchmark interactive responsiveness.** A [repeatable Blender benchmark](benchmarks/surface_responsiveness.py) requests six phenanthroline orbitals at 0.025 Å and counts actual viewport redraw callbacks while the worker computes. In three local Blender 4.3.2 trials per backend, the 95th-percentile redraw interval was 10.3 ms for both with redraws requested every 10 ms. From simulated Escape to observed worker cleanup took 10–12 ms in all six trials.
- [x] **Roll back failed geometry imports.** Discard geometry, collections, materials and animation data allocated by the failed import and restore selection, active object and mode. [Rollback tests](test/test_geometry_import.py)

## 9. Performance and numerical implementation

Entry points: [molecule.py](molecular_blender/molecule.py), [compute_ao_cy.pyx](molecular_blender/compute_ao_cy.pyx), [marching_cube_cy.pyx](molecular_blender/marching_cube_cy.pyx), and [marching_profile.py](marching_profile.py).

- [x] **Replace all-pairs connectivity with spatial neighbor search.** Cell lists find static and trajectory neighbors using the original element cutoffs and deterministic ordering. On a 20 × 20 × 3 carbon lattice at 1.45 Å spacing, both methods found 3,080 bonds. Median connectivity time over three local Blender 4.3.2 runs was 6.5 ms with cells and 407 ms with exhaustive search. Tests also compare trajectories against exhaustive search.
- [ ] **Benchmark alternative scene representations.** Offer instanced/Geometry Nodes atom rendering for large systems and compare it with linked objects. Preserve editable curve bonds as a supported mode; measure whether hook creation or object count is the actual bottleneck first.
- [x] **Remove repeated scene-wide deselection from trajectory keyframing.** On the 20-atom, four-frame benzopentyl example in Blender 4.3.2, median keyframing time over three imports fell from 0.79 to 0.27 ms in an empty scene and from 10.57 to 0.29 ms with 1,000 other objects. [Benchmark](benchmarks/trajectory_import.py)
- [x] **Bind imported bonds directly to their atoms.** Geometry imports and sidebar bond edits share direct endpoint hooks for full and split curves. Importing all 52 atoms of `examples/zine.molden` as sticks in Blender 4.3.2 took 0.66 s versus 1.44 s previously (median of three runs each). All 93 curves met their expected endpoints at peak displacement of mode 7.
- [x] **Remove repeated selection and parenting operators from geometry imports.** Atoms, charge spheres, gradients and rings are parented directly. The 20-atom trajectory imported in 13.8 ms versus 36.4 ms in an empty scene and 111.5 ms versus 455.0 ms with 1,000 extra objects (medians of three local Blender 4.3.2 runs). [Benchmark](benchmarks/trajectory_import.py)
- [x] **Avoid repeated atom name searches.** Preserve the first free suffix per element. Importing 1,000 carbon atoms in Blender 4.3.2 went from 1.30 s to 0.09 s in an empty scene and from 1.53 s to 0.11 s with another molecule present (medians of three).
- [x] **Share charge overlay geometry.** Each molecule has one charge collection. Spheres of the same element and sign share geometry and retain their own transforms and animation.
- [ ] **Reduce remaining geometry construction overhead.** Profile primitive creation, material assignment and dependency-graph updates during larger imports.
- [x] **Profile a bonded ball-and-stick import.** The [384-atom, 537-bond benchmark](benchmarks/bonded_import.py) measures atoms, curves, hooks and scene evaluation. In local Blender 4.3.2, skipping a redundant update for a static import with hooks off reduced total time from 93.9 to 70.2 ms (medians of three). With hooks enabled, total time was 124.0 ms, of which 55.0 ms was dependency-graph evaluation.
- [x] **Build indexed orbital and density meshes.** Molden calculations share vertices before passing surfaces to Blender. The initial import and sidebar recomputation write indexed faces directly, using the grid spacing to resolve duplicate intersections. [Mesh benchmark](benchmarks/orbital_mesh.py)
- [x] **Extend indexed meshes to cube data.** Cube and JSON field surfaces share vertices using each frame’s smallest voxel spacing. Tests cover rotated anisotropic grids, both signs of a difference field and field mapping on a remeshed density.
- [ ] **Audit Cython memory and floating-point policy.** Check allocation failures, pointer lifetime, array layout/dtype assumptions, and Python/Cython screening differences. Measure a float64 validation mode and the accuracy impact of `-ffast-math` before changing production flags.

## 10. Tests and scientific validation

There is already useful coverage for readers, numerical routines, geometry, bond orders, and an optional RDKit comparison. Extend those foundations rather than adding tests that only mirror implementation details.

- [x] **Test persisted scene edits and worker cancellation.** Tests cover save/reopen, undo/redo, renames, duplicated molecules, shared materials, scoped bond and surface deletion, and discarded worker results. [Bond lifecycle](test/test_bond_edit.py), [surface lifecycle](test/test_scene_surfaces.py) and [worker lifecycle](test/test_surface_jobs.py)
- [ ] **Test source refresh.** Cover compatible and incompatible changes when reload is implemented.
- [x] **Check transparency and transformed surface geometry in actual renders.** A fixed CPU Cycles scene checks the location, coverage and alpha of each exported surface in its PNG.
- [x] **Check mapped scalar fields in a render.** A CPU Cycles regression checks the red and blue sides of a remeshed density surface colored by an analytic potential.
- [x] **Check split-color multiple bonds in a render.** A CPU Cycles regression checks two distinct carbonyl strands, the color change at the midpoint, and coverage to the projected atom positions. [Render test](test/test_surface_render.py)
- [x] **Render and test gradient arrows.** A transparent Cycles render checks the shaft, head and projected tip. Timeline checks cover displacement, a change of direction and a zero vector.
- [x] **Check shell bounds and fractional density across backends.** Tests compare Python and Cython bounds for s–g shells and check fractional electron counts and density values. CI exercises the installed compiled wheels and an isolated Python fallback. [Backend parity tests](test/test_compute_ao.py)
- [x] **Check diffuse open-shell evaluator parity.** Compare point, plane and grid values and cumulative levels for a diffuse orbital and signed spin density across both backends.
- [x] **Compare final Python/Cython mesh geometry.** Both evaluator/mesher pairs produce closed signed lobes with matching component and face counts and vertices within 0.00002 Å on diffuse orbital and open-shell spin-density fixtures. The check exposed an absolute interpolation tolerance that distorted weak fields. Both meshers now interpolate independently of field scale.
- [x] **Fix the empty orbital surface in the pure-Python fallback.** Discovered on 2026-11-15: with compiled `.so`/`.pyd` extensions excluded from an isolated checkout, `TestMoldenIsosurface.test_small_orbital` in `test/test_isosurfaces.py` fails because orbital 1 of `examples/fake-atom.molden`, at isovalue `0.05` and resolution `0.3`, produces no triangles. Reproduced with Python 3.14 and NumPy 2.5.3; the full Blender 4.3.2 suite with compiled extensions passes. Fixed the shell-bound formula (`l` → `l + 1`), which collapsed s-shell bounds to zero. Added analytic and s–g backend-parity bounds checks; the complete isolated fallback suite passes, and fallback CI now runs full test discovery. Commit `614be0b`.
- [x] **Validate against analytic surfaces and electron counts.** Tests check the radius of a Gaussian isosurface, the location of a rotated cube plane, and integrated densities with fractional occupations. [Gaussian surface](test/test_isosurfaces.py), [cube plane](test/test_cube_datasets.py) and [density tests](test/test_orbitals.py)
- [ ] **Grow independent scientific references.** Add ESP reference values and expand fixtures for distorted geometry. The bond-perception suite already covers ions, resonance and invalid chemistry. Add a dedicated environment for its optional RDKit comparison.
- [x] **Test the installed artifacts, not just the checkout.** Standalone wheel tests require both compiled extensions from the installed package and exercise known Gaussian/surface results independently of test ordering. Missing installations/binaries fail; a separate CI job verifies Python fallback. Commit `9f085b1`.
- [ ] **Finish test isolation and performance baselines.** The test runner handles external working directories and paths with spaces. Reader tests exercise malformed Molden, SDF and cube input. Verify scene/datablock cleanup across the suite and add peak-memory and runtime baselines.

## 11. Packaging, compatibility, and user workflows

Entry points: [pyproject.toml](pyproject.toml), [setup.py](setup.py), [blender_setup.py](blender_setup.py), [build_all.sh](build_all.sh), and [CI workflows](.github/workflows/).

- [x] **Repair and verify wheel-test wiring.** Corrected the command to `test/test_wheel.py`, verified that it ships in a built sdist, and removed CI failure swallowing. The standalone test excludes checkout import precedence. Commit `300b77b`.
- [ ] **Unify build/install/export entry points.** Clarify the roles of `setup.py` and `blender_setup.py`, share extension/compiler configuration, and make Blender target discovery explicit. `--link-target` selection and destination preview are now implemented; sharing compiler configuration remains open.
- [ ] **Maintain an explicit Blender/Python/NumPy compatibility matrix.** CI currently pins Blender 3.6.7 while material code uses newer Principled socket names and transparency settings. Test the versions actually claimed, gate API differences by capability, and verify compiled modules inside each Blender runtime rather than trusting documentation's Python-version mapping.
- [x] **Test the release ZIP in Blender 4.3.2.** Linux x64 and Apple Silicon release jobs install the same assembled ZIP in an empty profile, check both native backends, and render an orbital. The ZIP includes Intel and Apple Silicon slices in each Mac extension.
- [ ] **Validate remaining supported release artifacts.** Install the addon ZIP in a clean Blender configuration, exercise a tiny import/render, and verify binary extension selection on each supported OS/architecture. Review pinned build tools and runner availability against the chosen matrix.
- [x] **Re-render a saved surface series from the command line.** Blender’s `molecular_render` command uses the saved camera, accepts a molecule and optional surface groups, and returns a failure exit code for invalid input or render errors.
- [x] **List saved surface groups.** `molecular_render --list` prints the active camera and each molecule’s surface labels and group IDs.
- [ ] **Inspect and validate inputs from the command line.** Report file metadata and validation failures with useful exit codes.
- [x] **Adapt import controls to the selected file and field.** The highlighted file determines the available controls. Cube files expose the appropriate dataset picker, Molden densities hide orbital controls, and JSON manifests preview the first frame. Import Details reports the active numerical backends. [Dialog tests](test/test_import_dialog.py)
- [x] **Add post-import and trajectory examples.** [Edit and export water orbitals](examples/edit_orbitals.py) and [animate an SDF record sequence](examples/water_bend.json).
- [x] **Add a runnable ESP figure workflow.** The [water example](examples/render_water_esp.py) maps a supplied PySCF potential cube onto its Molden density surface and exports a transparent PNG and matching SVG legend.
- [ ] **Add workflow examples for additional quantum-chemistry writers.**

## Bugs found during implementation

- [x] **Fix point-density evaluation and inconsistent tail screening.** The compiled single-point density call used an incompatible shell API. Point, plane and grid queries now agree on the diffuse open-shell fixture, including screened values.

- [x] **Close clipped orbital lobes.** Orbital meshes sample the full requested grid in bounded tiles to capture disconnected lobes. Adaptive density meshing follows fine crossings into adjacent cells. Resolution-scaled vertex welding keeps the raw Blender mesh closed. Verified on both signs of phenanthroline orbitals 46–49. [Signed-lobe regression](test/test_isosurfaces.py)

- [x] **Use the grid determinant for cube voxel volume.** Cube containment work found that column-norm products overestimate rotated anisotropic voxel volumes. VolumeData now uses the absolute determinant, with an independent rotated-grid regression.


- [x] **Count orbital trajectory tasks in surface progress.** Density-control work found missing trajectory ticks and incomplete totals; field counts now include orbital snapshots.


- [x] **Prevent symbolic orbital expressions from wrapping into relative indexing.** Validation review found that `homo-N` resolving to zero was interpreted again as HOMO, and partial expressions were accepted by substring matching. Symbolic requests now require full matches and valid absolute bounds, with a Blender regression. Commit `7c3fbf1`.


- [x] **Implement restricted spatial-orbital spin densities.** Population validation reproduced a scientific bug: a closed-shell Molden input with Alpha labels and occupation 2 was treated as entirely alpha density, producing nonzero spin density. Restricted and fractional spatial inputs now use the approved visualization defaults: majority-spin filling up to 1, with the remainder in the other spin; doubly occupied orbitals split equally. Completed with commit `f5b6429` and visualization defaults approved 2026-11-16.


- [x] **Choose an Ubuntu image that supports the wheel-test Python matrix.** Unblocking macOS builds exposed a later failure in the `test_wheels` matrix: `actions/setup-python` reports `The version '3.9' with architecture 'x64' was not found for Ubuntu 26.04.` The floating `ubuntu-latest` runner now selects Ubuntu 26.04. Choose a supported image for older Python tests or explicitly revise the support policy; keep matrix failures independent so one setup failure does not cancel other platforms. Reproduced on final commit `a9b590e`: https://github.com/smparker/molecular-blender/actions/runs/43757795670/job/132948634610 Commit `466cd05`. Final Blender and wheel CI passed on `7c3fbf1`, including Python 3.9 on Ubuntu 22.04.

- [x] **Update unsupported macOS CI runner labels.** The pushed `451f479` wheel run on 2026-11-15 cancelled its `macos-14` job before executing any steps: GitHub reports `The configuration 'macos-14-arm64-us-default' is not supported`. Its `macos-13` job remains queued, including on the prior run. Updated the matrix and architecture expression to `macos-15-intel` and `macos-15`; both native wheel builds passed again on final commit `a9b590e` in run 43757795670. Run: https://github.com/smparker/molecular-blender/actions/runs/43753231680

- [x] **Fix the Python marching-cubes Z offset.** Discovered during analytic surface checks: `marching_cube_box` reads voxel planes `k-1`/`k` but places them at `k`/`k+1`, translating fallback surfaces by one Z voxel. Reproduce with a 2×2×2 field equal to its Z index and isovalue 0.25; the surface must lie at Z=0.25. Commit `ff53db4`.

- [x] **Correct the duplicated AO index in `fake-atom.molden`.** Stricter Molden coefficient validation exposed a fixture typo: line 384 labelled coefficient 20 as 10. Corrected the label without changing any numerical coefficient. Commit `476f7d1`.

- [x] **Honor the density volume selector.** Inspection while preparing importer diagnostics found that `compute_surface_vertex_sets` chooses the field only from `options["orbital"]`; selecting `volume="density"` still samples the requested orbital. Density selection now samples total density and ignores orbital text, with matching progress counts and an end-to-end Blender regression. Section 3 records the subsequent spin controls and population validation. Commit `6431742`.

- [x] **Harden Blender configuration discovery.** Installer review found that Linux requires `XDG_CONFIG_HOME` even when unset, absent configuration directories raise instead of being skipped, and float-based version sorting places 4.10 before 4.9. Handle missing roots and sort version components numerically. Commit `6a8740f`.

## Completed small increments

- [x] **Validate occupations and spin labels.** Check orbital counts and array lengths, finite occupations in [0,2], and normalized complete Alpha/Beta labels. Per the approved visualization policy, infer scale from maximum occupancy and use deterministic alpha-majority defaults when labels are absent; densities preserve fractional weights. Commit `3031ff6`.


- [x] **Validate orbital and density requests.** Reject noninteger and out-of-range orbital indices, undefined HOMO-relative requests, and unknown density names while retaining valid relative indexing. Commit `28b3e50`.


- [x] **Harden cumulative-isovalue edge cases.** Signed fractions, zero/non-finite/empty fields, voxel volume, iteration budgets, and tolerances are validated; extreme-scale and tiny-grid regressions pass. Thresholds retain precision and tied samples are treated as groups. Commit `500a232` subsequently added cube field-type selection. Commit `85e8e7d`.


- [x] **Preserve fractional electron counts.** The compiled density constructor stores a public double electron count; rebuilt-extension parity and integrated half-electron density regressions pass. Commit `cd3a8cf`.


- [x] **Preserve orbital energies and symmetry labels.** `MOData` retains metadata in input order with explicit missing entries; frontier selection followed in `0a287ce`. Commit `4b3dd10`.


- [x] **Remove the disabled Babel fallback.** Commit `092fb42`; unsupported formats raise an actionable `ValueError`.
- [x] **Repair test-runner paths and argument handling.** Commit `717c68b`; named tests and discovery run from outside the checkout, shell arguments remain intact, and Blender starts with factory settings. Broader scene cleanup and performance-baseline work in section 10 remains open.

- [x] **Fix Python charge-option defaults.** `default_options` now provides `charge_offset`/`charge_factor`, and a Blender regression imports charge spheres using those defaults. Commit `e8a2ad9`.

- [x] **Free temporary surface BMesh objects.** Mesh welding now frees its BMesh in a `finally` block, including on failure. Commit `34def23`.

- [x] **Distinguish missing optional modules from broken imports.** Missing numerical extensions use Python fallbacks; initialization errors, missing internal dependencies, and broken Blender imports propagate their original errors. Standalone subprocess regressions cover each case. Commit `9bed92f`.

- [x] **Show supported import formats and active backends.** Updated the import menu and addon description; Import Details reports the actual AO and meshing backend separately, with README troubleshooting guidance. Commit `0f16de7`.

- [x] **Honor installer targets and preview replacements.** `--link-target` accepts discovered versions in both copy/link modes; `--dry-run` previews destinations, and invalid targets are rejected before mutation. Tests use temporary installations, including dangling symlinks. README distinguishes addon and wheel entry points. Commit `7a69662`.

- [x] **Run the complete Python fallback suite in CI.** With the fallback surface bugs fixed, the isolated CI job now runs all tests, including backend loading, parsers, runner, and installer regressions; Blender-only and compiled-only checks skip explicitly. Commit `451f479`.

## Scope notes

The first post-import editing increment does not require replacing the renderer. External ESP cubes provide field maps without a Gaussian Coulomb integral evaluator. Full aromaticity, metal chemistry, periodic systems, and native ESP evaluation each deserve separate designs and validation datasets before implementation.

Visualization defaults approved 2026-11-16: infer occupation scale 1/2 from the maximum; classify occupied strictly above half (half is empty). Spatial populations fill the majority spin first; missing labels default to alpha. Preserve actual fractional populations for density calculations.

Verification for `0a287ce` (2026-11-16): all four scientific increments are committed and pushed. Full local Blender/compiled-Python/fallback suites each passed 280 tests with expected skips. GitHub Blender run 43931902956 and wheel run 43931902811 both completed successfully.
