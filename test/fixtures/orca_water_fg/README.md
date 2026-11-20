# Spherical f and g orbital reference

ORCA 6.1.1 calculated RHF water with STO-3G and one uncontracted f and one
uncontracted g shell on oxygen. The tilted molecular geometry mixes angular
components in the molecular orbitals. The input requests a cube for each of
the 23 orbitals on the same 13 × 13 × 13 grid. This small basis was chosen to
exercise both shells and their phases.

To reproduce the files, run these commands in a scratch directory containing
`water.inp`.

```sh
orca water.inp > water.out
orca_2mkl water -molden
```

`water.molden` contains the resulting `water.molden.input` with trailing
whitespace removed. Each row of `mo-values.dat` gives x, y and z in Bohr,
followed by the signed amplitudes in Bohr⁻³ᐟ² read from `mo0.cube` through
`mo22.cube` in order. These correspond to Molden orbitals 1 through 23.
The first voxel is at (−3, −3, −3), the spacing is 0.5 Bohr and z varies
fastest. The final two rows probe weak tails near the corners of the grid.
The test suite reads these stored samples directly.
