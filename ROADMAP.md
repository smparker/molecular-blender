# Molecular Blender roadmap

An editable menu of possible work, based on a code review on 2026-11-15. These are suggestions, not release commitments. Completed increments are checked; broad items remain open until their full scope is complete. Existing features are mentioned as foundations, not proposed as new work. Items described as audits need focused reproduction and tests before changing behavior.

Keep the current curve-extrusion approach for bonds as the baseline. Keep bond-order drawing and perception opt-in until their behavior has broader validation.

## Suggested starting points

- **Make imports editable:** persistent molecule identity and settings → a selected-molecule panel → restyle and regenerate operators.
- **Build ESP maps incrementally:** general scalar-field sampling → map an imported ESP grid onto a density surface → evaluate ESP directly from Molden wavefunctions.
- **Remove correctness traps:** option defaults, element-dependent bond radii, numerical units, cancellation, and wheel validation.
- **Improve chemistry incrementally:** actionable perception diagnostics and manual overrides before expanding the supported chemistry.

## 1. Post-import editing and scene interaction

Current foundations: linked atom geometry, element materials, molecule collections, endpoint hooks, and a direct `set_bond_order` helper. Most settings and the `Molecule` instance currently live only during import. Start in [molecule.py](molecular_blender/molecule.py), [plotter.py](molecular_blender/plotter.py), and [the importer UI](molecular_blender/__init__.py).

- [ ] **Persist molecule identity and settings in the `.blend`.** Store a versioned settings block, source path, molecule ID, atom IDs, bond endpoints, explicit/perceived order provenance, and surface definitions. Resolve objects by IDs rather than requiring their original names; survive renaming, save/reopen, and duplication.
- [ ] **Add a selected-molecule sidebar panel.** Show source, atom/bond counts, units, charge, representation, and available orbital/field data. Make the owning molecule discoverable when an atom, bond, or surface is selected.
- [ ] **Restyle without reimporting.** Change atom size, bond thickness, multiple-bond spacing, palette, opacity, and surface appearance for one molecule or the current selection. Update custom bevel profiles as well as round bonds, and avoid changing unrelated objects that share materials.
- [ ] **Add an orbital/surface browser.** The sidebar now lists saved orbital indices, energies, occupations, spin labels, and symmetries. Add controls for choosing a density and adding, replacing, hiding, or removing surfaces. Appearance-only changes should not recompute the field.
- [x] **Expose bond editing tools.** Set an order, add/delete a bond, lock a manual assignment against reperception, and adjust strand tilt. Split bonds update both halves. Order changes update junction sizes, and profile rebuilds preserve manual tilt.
- [ ] **Add atom-aware selection and inspection.** The selected atom now shows its element, index, partial charge, and formal charge. Add selection by element, fragment, index, or bonded neighbors, optional atom labels, and distance, angle, and dihedral measurements.
- [ ] **Refresh source data while preserving presentation.** Reload coordinates or orbitals and retain user colors, visibility, transforms, and locked edits where atom identities still match. Report incompatible atom lists rather than silently rebinding objects.
- [ ] **Make edits undoable and repeatable.** Provide explicit Apply, Rebuild Geometry, Recompute Surface, and Reset Style operations, with molecule-scoped cleanup. Define whether a manual move updates the chemistry model and whether moving atoms invalidates wavefunction-derived surfaces.
- [ ] **Save reusable visual presets.** Store/import/export molecular styles and surface color ranges; include a “make this molecule independent” action for linked geometry/materials.

## 2. Electrostatic-potential maps

Current foundations: [MOData.get_density](molecular_blender/orbitals.py) constructs AO density matrices; the Python/Cython evaluators sample orbitals and densities; [isosurfaces.py](molecular_blender/isosurfaces.py) produces surface geometry. There is no ESP evaluator or general scalar-to-surface color mapping yet.

- [ ] **Start with ESP supplied as a grid.** Import an ESP cube alongside a density cube or Molden model, interpolate ESP onto the density-surface vertices, and display a diverging color scale. Validate alignment, grid axes, units, and coverage; do not treat an ESP grid as electron density.
- [ ] **Add a general field-mapped surface material.** Preserve numerical values as mesh attributes, with editable color ramp, explicit min/max, optional symmetric limits, clipping, legend, and units. Sample on the final displayed mesh or transfer attributes correctly through remeshing.
- [ ] **Implement a wavefunction ESP evaluator.** In atomic units, evaluate `V(r) = Σ_A Z_A / |r-R_A| - ∫ ρ(r') / |r-r'| dr'`. Reuse the AO density matrix, but add the electron–point Coulomb integrals; existing pointwise density evaluation does not supply this integral. Compare a native Gaussian-integral backend with an optional external calculation path. [PySCF's MEP example](https://github.com/pyscf/pyscf/blob/master/examples/1-advanced/031-MEP.py) is a useful implementation/reference calculation.
- [ ] **Define the scientific input contract before claiming general Molden ESP support.** Check occupations, orbital completeness, normalization, nuclear charges, ghost centers, and effective-core-potential treatment. Preserve needed atom metadata; require missing core-charge information or explicitly restrict the first implementation to validated all-electron inputs. See the [Molden format](https://www.theochem.ru.nl/molden/molden_format.html).
- [ ] **Keep exact and approximate electrostatics distinct.** A quick map from supplied partial charges could be useful, but label it as a point-charge approximation and retain its charge source. Do not substitute it silently for a wavefunction-derived ESP.
- [ ] **Validate ESP numerically.** Compare selected points and surface values for neutral and ionic reference molecules against an independent code; test the nuclear/electronic sign convention, unit conversion, distance from nuclear singularities, and far-field total-charge behavior.
- [ ] **Make ESP computation practical.** Evaluate surface points in cancellable batches, screen negligible shell pairs where justified, cache the density matrix and evaluated values, and allow coarse preview followed by a final calculation. Benchmark both time and peak memory before optimizing kernels.

## 3. Other scientific visualization features

Entry points: [orbitals.py](molecular_blender/orbitals.py), [compute_ao_py.py](molecular_blender/compute_ao_py.py), [compute_ao_cy.pyx](molecular_blender/compute_ao_cy.pyx), [isosurfaces.py](molecular_blender/isosurfaces.py), and [stylers.py](molecular_blender/stylers.py).

- [x] **Expose existing density modes properly.** Provide explicit total, alpha, beta, and spin-density controls; validate restricted/unrestricted inputs and explain when spin information is absent. These evaluators already exist, but their selection is primarily text-driven. Completed with commit `ba0d216` and visualization defaults approved 2026-11-16.
- [x] **Show orbital metadata and support spin-specific frontier orbitals.** Retain energies and symmetry labels in `MOData`; select HOMO/LUMO separately for alpha/beta sets and handle fractional occupations and unsorted orbitals deliberately. Completed with commit `0a287ce` and visualization defaults approved 2026-11-16.
- [ ] **Support field differences.** Compute density/ESP differences between two compatible calculations, with explicit alignment and grid resampling. Show positive and negative contributions with a shared scale.
- [ ] **Import and animate normal modes.** Read Molden frequency/displacement sections and expose mode, amplitude, phase, and playback controls; keep the equilibrium structure intact. The [format already defines these sections](https://www.theochem.ru.nl/molden/molden_format.html), but the reader currently does not consume them.
- [ ] **Improve publication output.** Add molecule framing, principal-axis/selected-plane alignment, orthographic camera presets, scale bars, legends, transparent-background render presets, and batch figure export without modifying an existing camera unexpectedly.

## 4. Bond perception and chemical representation

Current implementation: covalent-radius connectivity, bounded valence/charge search, explicit SDF precedence, and fixed initial strand orientation. Start in [bond_perception.py](molecular_blender/bond_perception.py), [molecule.py](molecular_blender/molecule.py), and [aromatics.py](molecular_blender/aromatics.py).

- [x] **Make failures explainable.** Report the offending atom IDs/elements, degrees, candidate valences, requested charge, and search budget. Distinguish bad connectivity, unsupported chemistry, missing hydrogens, incompatible charge, and exhausted search. Commit `a4c6961`.
- [ ] **Allow manually specifying bonds and bond orders after import.** Should be able to select two atoms and then click a button to add a single, double, or whatever else.
- [ ] **Scale the solver and expose ambiguity.** Partition independent graph work where possible, propagate charge constraints across fragments, add timing/node diagnostics, and benchmark symmetric/fused systems. Preserve search limits and atomic application of successful results.
- [ ] **Keep multiple-bond orientation coherent.** Use a ring/molecular plane where appropriate, allow manual overrides, and maintain frame-to-frame continuity instead of flipping with the strongest neighbor. Avoid unstable orientation for nearly collinear geometry.

## 5. Import, export, and data fidelity

Entry points: [importers.py](molecular_blender/importers.py), [transform.py](molecular_blender/transform.py), and [turbomole_reader.py](molecular_blender/turbomole_reader.py).

- [x] **Harden Molden basis parsing.** Expands `sp` shells, registers standalone `[7F]`, maps GTO blocks by atom reference, reads Fortran D/d numbers, and validates indexed coefficient blocks and section endings. Existing molecular fixtures and focused synthetic cases pass. Commit `f1547e0`.
- [ ] **Make normalization and AO ordering explicit.** Replace a loosely explained producer workaround with documented conventions and validated transformations for Cartesian/spherical shells. Include independent orbital-overlap and density checks across supported quantum-chemistry writers.
- [x] **Preserve all cube datasets.** All cube fields and dataset IDs are retained, with a zero-based importer/API selector and first-field compatibility. Optional header value counts and malformed lengths are tested. Free-form comments are retained without guessing physical field types or units.
- [ ] **Audit cube coordinates and grid conventions.** Test positive/negative voxel-count conventions, origins, atom coordinates, rotated orthogonal axes, optional header values, and malformed data lengths. Support affine/sheared grids where justified, with voxel volume from the grid determinant rather than orthogonal-axis assumptions.
- [ ] **Consider PDB/mmCIF, MOL2, and periodic inputs.** Choose formats based on actual use cases: biomolecular hierarchy, typed bonds, or unit cells. The disabled Babel fallback has been removed; unsupported formats now report the supported formats.
- [ ] **Make JSON trajectory manifests consistent with native imports.** Validate referenced paths, atom identity, frame counts, and compatible basis/grid data; decide how SDF records and explicit topology fit the manifest instead of silently dropping information.

## 6. Trajectories and interactive performance

Entry points: `AnimateMolecule`, `compute_surface_vertex_sets`, and `create_surface_meshes` in [plotter.py](molecular_blender/plotter.py).

- [ ] **Expose trajectory sampling.** Choose start/end frames, stride, time units, frame rate, and interpolation policy; stream or downsample long trajectories before allocating Blender animation data.
- [ ] **Stabilize topology over time.** Add bond-formation/breaking hysteresis near cutoffs. Design and test per-frame bond-order perception separately from the existing fixed-first-frame assignment, including strand appearance/disappearance and charge conservation.
- [ ] **Improve orbital animation continuity.** Track orbital overlap/phase across frames and expose the limitations around crossings and near-degenerate orbitals rather than trusting an unchanged index alone.
- [ ] **Avoid retaining every surface frame in memory.** Offer lazy calculation, a bounded disk cache, and staged scene updates instead of storing all triangle sets and all mesh objects simultaneously.

## 7. Architecture and data ownership

The existing split between field computation and main-thread mesh creation is a useful starting point; extend it instead of discarding it.

- [ ] **Define one validated options model.** Share names, defaults, units, and validation across the UI, Python API, saved settings, and CLI. The charge-option spelling mismatch is fixed below; sharing and validating the complete model remains open.
- [ ] **Separate scientific data from Blender objects.** Keep parsing, topology, perception, scalar fields, and numerical calculations usable without `bpy`/`mathutils`; let a scene adapter own object/material names and handles. Avoid a missing Blender-only import being indistinguishable from a broken package import.
- [ ] **Split `plotter.py` by responsibility.** Extract atom/bond geometry, surface generation, animation, scene ownership, and driver orchestration behind small interfaces. Factor the duplicated split/unsplit curve setup while preserving the tested extrusion and hook behavior.
- [ ] **Introduce a scalar-field interface.** Standardize point/batch/grid evaluation, bounds, units, sign semantics, provenance, and optional derivatives for cube fields, orbitals, densities, and future ESP. Keep numerical values independent of color/material choices.
- [ ] **Use structured errors and diagnostics.** Carry filename/record/section/atom context from readers and chemistry into Blender reports; retain detailed tracebacks in logs. Replace broad exception swallowing and misleading fallback messages with precise causes.
- [ ] **Make object/material reuse explicit.** Cache by molecule/style identity and validated parameters, not just global names. Provide safe cleanup for molecule-owned meshes, materials, actions, bevel profiles, and collections, including fake-user profiles.
- [ ] **Clarify numerical APIs.** Add type hints and shape/unit contracts at module boundaries; replace ambiguous integer indexing with consistent atom IDs and positional indices. Fix misleading names/docstrings and retain compatibility aliases where public names such as `OrbitalCalculater` are already used.

## 8. Correctness and lifecycle fixes to investigate first

These are concrete targets spotted in the current code, not a claim that every path has been reproduced in a running Blender session.

- [x] **Fix element-dependent bond-radius lookup.** Uses the selected bond-style element consistently; mixed-element Blender regressions cover radius selection and unaffected styles. Commit `dbbafb3`.
- [x] **Audit isosurface method forwarding and resolution units.** Forwards the requested method, converts both adaptive spacings from angstrom to Bohr, and rejects invalid methods/spacings. Analytic Gaussian-sphere regressions check displayed geometry and final sampling spacing with both backends. Commit `537bf83`.
- [x] **Harden cumulative-isovalue calculations.** Validate fractions, signs, all-zero fields, non-finite inputs, tiny grids, and termination tolerances. Define whether a field integrates `ψ²`, density, or absolute signed density; an untyped cube should not silently inherit orbital-square semantics. Completed with commit `500a232` and visualization defaults approved 2026-11-16.
- [x] **Audit occupations and spin densities.** Validate absent/inconsistent spin labels, fractional occupations, electron counts, and alpha/beta frontier indexing. Check the Cython density constructor's integer electron-count argument against fractional-occupation inputs.
- [ ] **Implement real cancellation and worker cleanup.** Explicitly handle Escape, propagate a cancellation signal into long loops, account for orbital trajectory tasks in progress, and prevent stale results from reaching a cancelled/deleted molecule. Restore timers/status and free temporary data on every exit path.
- [ ] **Review the background execution model.** The operator currently runs a daemon thread while Blender remains active; “no `bpy` calls in the worker” alone does not establish support under Blender's [threading guidance](https://docs.blender.org/api/5.0/info_gotchas_threading.html). Evaluate a cancellable subprocess worker or chunked main-thread jobs and make pure-Python/Cython responsiveness measurable.
- [ ] **Make scene generation transactional.** Restore selection, active object, and mode after operations; roll back partial imports on errors. Test translated/rotated/scaled molecule parents and moved atoms for bonds, gradients, and surfaces. Temporary BMesh objects are now freed; transactional scene rollback and context restoration remain open.

## 9. Performance and numerical implementation

Entry points: [molecule.py](molecular_blender/molecule.py), [compute_ao_cy.pyx](molecular_blender/compute_ao_cy.pyx), [marching_cube_cy.pyx](molecular_blender/marching_cube_cy.pyx), and [marching_profile.py](marching_profile.py).

- [ ] **Replace all-pairs connectivity with spatial neighbor search.** Use a cell list or optional spatial index, including trajectory handling, while retaining the existing cutoff semantics and deterministic edge ordering.
- [ ] **Benchmark alternative scene representations.** Offer instanced/Geometry Nodes atom rendering for large systems and compare it with linked objects. Preserve editable curve bonds as a supported mode; measure whether hook creation or object count is the actual bottleneck first.
- [ ] **Reduce operator/context overhead.** Use the Blender data API for predictable construction and parenting where possible; batch linking, material assignment, and dependency-graph updates.
- [ ] **Build indexed surface meshes directly.** Reduce triangle-soup duplication and repeated welding; make weld tolerance and remeshing resolution-aware. Validate seams, normals, disconnected components, surface area, and field-attribute transfer.
- [ ] **Share numerical preparation and cache work.** Reuse grids/AO blocks across compatible orbital/density requests, use memory-bounded batches, and invalidate caches by source/basis/geometry/settings rather than names alone.
- [ ] **Audit Cython memory and floating-point policy.** Check allocation failures, pointer lifetime, array layout/dtype assumptions, and Python/Cython screening differences. Measure a float64 validation mode and the accuracy impact of `-ffast-math` before changing production flags.

## 10. Tests and scientific validation

There is already useful coverage for readers, numerical routines, geometry, bond orders, and an optional RDKit comparison. Extend those foundations rather than adding tests that only mirror implementation details.

- [ ] **Add scene-edit lifecycle tests.** Cover save/reopen, undo/redo, renaming, duplicate imports, shared-material edits, molecule-scoped deletion, source reload, and cancellation/error rollback.
- [ ] **Add controlled render regressions.** Use small fixed scenes for split-color multiple bonds, transparency, gradients, mapped scalar fields, and parent transforms. Pair image checks with geometric/numerical assertions and keep renderer/version differences explicit.
- [ ] **Strengthen Python/Cython parity.** Force each backend in CI and compare shell values, orbital/density grids, screening, bounds, cumulative levels, and surface geometry across s–g functions and diffuse/fractional/open-shell cases.
- [x] **Fix the empty orbital surface in the pure-Python fallback.** Discovered on 2026-11-15: with compiled `.so`/`.pyd` extensions excluded from an isolated checkout, `TestMoldenIsosurface.test_small_orbital` in `test/test_isosurfaces.py` fails because orbital 1 of `examples/fake-atom.molden`, at isovalue `0.05` and resolution `0.3`, produces no triangles. Reproduced with Python 3.14 and NumPy 2.5.3; the full Blender 4.3.2 suite with compiled extensions passes. Fixed the shell-bound formula (`l` → `l + 1`), which collapsed s-shell bounds to zero. Added analytic and s–g backend-parity bounds checks; the complete isolated fallback suite passes, and fallback CI now runs full test discovery. Commit `614be0b`.
- [ ] **Grow independent scientific references.** Add electron-count/overlap tests, exact or analytic fields with known surfaces, ESP reference values, and molecule fixtures covering charges, resonance, distorted geometry, and deliberately invalid chemistry. Run the optional RDKit comparison in a dedicated environment so it is not always skipped in Blender.
- [x] **Test the installed artifacts, not just the checkout.** Standalone wheel tests require both compiled extensions from the installed package and exercise known Gaussian/surface results independently of test ordering. Missing installations/binaries fail; a separate CI job verifies Python fallback. Commit `9f085b1`.
- [ ] **Make test execution reproducible.** Resolve fixtures/imports independently of the current directory, quote shell paths, preserve test arguments as arrays, support Blender paths with spaces, and ensure tests leave scene/datablock state clean. Add targeted malformed/truncated-input tests and performance baselines with peak-memory reporting.

## 11. Packaging, compatibility, and user workflows

Entry points: [pyproject.toml](pyproject.toml), [setup.py](setup.py), [blender_setup.py](blender_setup.py), [build_all.sh](build_all.sh), and [CI workflows](.github/workflows/).

- [x] **Repair and verify wheel-test wiring.** Corrected the command to `test/test_wheel.py`, verified that it ships in a built sdist, and removed CI failure swallowing. The standalone test excludes checkout import precedence. Commit `300b77b`.
- [ ] **Unify build/install/export entry points.** Clarify the roles of `setup.py` and `blender_setup.py`, share extension/compiler configuration, and make Blender target discovery explicit. `--link-target` selection and destination preview are now implemented; sharing compiler configuration remains open.
- [ ] **Maintain an explicit Blender/Python/NumPy compatibility matrix.** CI currently pins Blender 3.6.7 while material code uses newer Principled socket names and transparency settings. Test the versions actually claimed, gate API differences by capability, and verify compiled modules inside each Blender runtime rather than trusting documentation's Python-version mapping.
- [ ] **Validate every supported release artifact.** Install the addon ZIP in a clean Blender configuration, exercise a tiny import/render, and verify binary extension selection on each supported OS/architecture. Review pinned build tools and runner availability against the chosen matrix.
- [ ] **Add a stable headless CLI.** Build on `BlendMoleculeSetup` with commands for inspect, validate, import, render, and export; support a settings file and repeatable batch runs with useful exit codes. Carry forward the existing command-line simplification idea from [improvement_ideas.md](improvement_ideas.md).
- [ ] **Improve discovery and troubleshooting.** The format-aware import label and active-backend diagnostics are complete; add small examples for post-import edits, ESP, trajectories, and common quantum-chemistry writers. Show active numerical backends and actionable import diagnostics without assuming an installed compiler.

## Bugs found during implementation

- [x] **Use the grid determinant for cube voxel volume.** Cube containment work found that column-norm products overestimate rotated anisotropic voxel volumes. VolumeData now uses the absolute determinant, with an independent rotated-grid regression.


- [x] **Count orbital trajectory tasks in surface progress.** Density-control work found missing trajectory ticks and incomplete totals; field counts now include orbital snapshots.


- [x] **Prevent symbolic orbital expressions from wrapping into relative indexing.** Validation review found that `homo-N` resolving to zero was interpreted again as HOMO, and partial expressions were accepted by substring matching. Symbolic requests now require full matches and valid absolute bounds, with a Blender regression. Commit `7c3fbf1`.


- [x] **Implement restricted spatial-orbital spin densities.** Population validation reproduced a scientific bug: a closed-shell Molden input with Alpha labels and occupation 2 was treated as entirely alpha density, producing nonzero spin density. Restricted and fractional spatial inputs now use the approved visualization defaults: majority-spin filling up to 1, with the remainder in the other spin; doubly occupied orbitals split equally. Completed with commit `f5b6429` and visualization defaults approved 2026-11-16.


- [x] **Choose an Ubuntu image that supports the wheel-test Python matrix.** Unblocking macOS builds exposed a later failure in the `test_wheels` matrix: `actions/setup-python` reports `The version '3.9' with architecture 'x64' was not found for Ubuntu 26.04.` The floating `ubuntu-latest` runner now selects Ubuntu 26.04. Choose a supported image for older Python tests or explicitly revise the support policy; keep matrix failures independent so one setup failure does not cancel other platforms. Reproduced on final commit `a9b590e`: https://github.com/smparker/molecular-blender/actions/runs/43757795670/job/132948634610 Commit `466cd05`. Final Blender and wheel CI passed on `7c3fbf1`, including Python 3.9 on Ubuntu 22.04.

- [x] **Update unsupported macOS CI runner labels.** The pushed `451f479` wheel run on 2026-11-15 cancelled its `macos-14` job before executing any steps: GitHub reports `The configuration 'macos-14-arm64-us-default' is not supported`. Its `macos-13` job remains queued, including on the prior run. Updated the matrix and architecture expression to `macos-15-intel` and `macos-15`; both native wheel builds passed again on final commit `a9b590e` in run 43757795670. Run: https://github.com/smparker/molecular-blender/actions/runs/43753231680

- [x] **Fix the Python marching-cubes Z offset.** Discovered during analytic surface checks: `marching_cube_box` reads voxel planes `k-1`/`k` but places them at `k`/`k+1`, translating fallback surfaces by one Z voxel. Reproduce with a 2×2×2 field equal to its Z index and isovalue 0.25; the surface must lie at Z=0.25. Commit `ff53db4`.

- [x] **Correct the duplicated AO index in `fake-atom.molden`.** Stricter Molden coefficient validation exposed a fixture typo: line 384 labelled coefficient 20 as 10. Corrected the label without changing any numerical coefficient. Commit `476f7d1`.

- [x] **Honor the density volume selector.** Inspection while preparing importer diagnostics found that `compute_surface_vertex_sets` chooses the field only from `options["orbital"]`; selecting `volume="density"` still samples the requested orbital. Density selection now samples total density and ignores orbital text, with matching progress counts and an end-to-end Blender regression. Restricted/unrestricted validation and explicit spin controls remain in section 3. Commit `6431742`.

- [x] **Harden Blender configuration discovery.** Installer review found that Linux requires `XDG_CONFIG_HOME` even when unset, absent configuration directories raise instead of being skipped, and float-based version sorting places 4.10 before 4.9. Handle missing roots and sort version components numerically. Commit `6a8740f`.

## Completed small increments

- [x] **Validate occupations and spin labels.** Check orbital counts and array lengths, finite occupations in [0,2], and normalized complete Alpha/Beta labels. Per the approved visualization policy, infer scale from maximum occupancy and use deterministic alpha-majority defaults when labels are absent; densities preserve fractional weights. Commit `3031ff6`.


- [x] **Validate orbital and density requests.** Reject noninteger and out-of-range orbital indices, undefined HOMO-relative requests, and unknown density names while retaining valid relative indexing. Commit `28b3e50`.


- [x] **Harden cumulative-isovalue edge cases.** Signed fractions, zero/non-finite/empty fields, voxel volume, iteration budgets, and tolerances are validated; extreme-scale and tiny-grid regressions pass. Thresholds retain precision and tied samples are treated as groups. Field-type selection for untyped cube inputs remains open. Commit `85e8e7d`.


- [x] **Preserve fractional electron counts.** The compiled density constructor stores a public double electron count; rebuilt-extension parity and integrated half-electron density regressions pass. Commit `cd3a8cf`.


- [x] **Preserve orbital energies and symmetry labels.** `MOData` retains metadata in input order with explicit missing entries; frontier-orbital selection remains open. Commit `4b3dd10`.


- [x] **Remove the disabled Babel fallback.** Commit `092fb42`; unsupported formats raise an actionable `ValueError`.
- [x] **Repair test-runner paths and argument handling.** Commit `717c68b`; named tests and discovery run from outside the checkout, shell arguments remain intact, and Blender starts with factory settings. Broader scene cleanup and performance-baseline work in section 10 remains open.

- [x] **Fix Python charge-option defaults.** `default_options` now provides `charge_offset`/`charge_factor`, and a Blender regression imports charge spheres using those defaults. Commit `e8a2ad9`.

- [x] **Free temporary surface BMesh objects.** Mesh welding now frees its BMesh in a `finally` block, including on failure. Commit `34def23`.

- [x] **Distinguish missing optional modules from broken imports.** Missing numerical extensions use Python fallbacks; initialization errors, missing internal dependencies, and broken Blender imports propagate their original errors. Standalone subprocess regressions cover each case. Commit `9bed92f`.

- [x] **Show supported import formats and active backends.** Updated the import menu and addon description; Import Details reports the actual AO and meshing backend separately, with README troubleshooting guidance. Commit `0f16de7`.

- [x] **Honor installer targets and preview replacements.** `--link-target` accepts discovered versions in both copy/link modes; `--dry-run` previews destinations, and invalid targets are rejected before mutation. Tests use temporary installations, including dangling symlinks. README distinguishes addon and wheel entry points. Commit `7a69662`.

- [x] **Run the complete Python fallback suite in CI.** With the fallback surface bugs fixed, the isolated CI job now runs all tests, including backend loading, parsers, runner, and installer regressions; Blender-only and compiled-only checks skip explicitly. Commit `451f479`.

## Scope notes

The first post-import editing increment does not require replacing the renderer. The first ESP increment does not require implementing Gaussian Coulomb integrals: importing an external ESP grid would validate the field-mapping UI and data model first. Full aromaticity, metal chemistry, periodic systems, native ESP evaluation, and orbital tracking each deserve separate designs and validation datasets before implementation.

Visualization defaults approved 2026-11-16: infer occupation scale 1/2 from the maximum; classify occupied strictly above half (half is empty). Spatial populations fill the majority spin first; missing labels default to alpha. Preserve actual fractional populations for density calculations.

Verification for `0a287ce` (2026-11-16): all four scientific increments are committed and pushed. Full local Blender/compiled-Python/fallback suites each passed 280 tests with expected skips. GitHub Blender run 43931902956 and wheel run 43931902811 both completed successfully.
