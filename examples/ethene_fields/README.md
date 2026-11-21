# Ethene from TURBOMOLE and xTB

Both calculations use the neutral singlet at the geometry in [ethene.xyz](ethene.xyz).
The two figures show different physical quantities. The TURBOMOLE Molden file
contains a signed wavefunction and the xTB cube contains the valence electron density.

## TURBOMOLE

With TURBOMOLE 8.1, convert the coordinates with `x2t ethene.xyz > coord`.
In `define`, read `coord`, use C1 symmetry, assign def2-SVP to every atom, and
make a closed-shell neutral EHT guess. Run `dscf`, then `tm2molden batch`.
[turbomole.molden](turbomole.molden) is that default export, renamed from
`molden.input`. Eight orbitals are occupied. The figure renders orbital 8 at
±0.11 bohr⁻³ᐟ² with the legacy Cartesian normalization.

## xTB

Run `xtb ethene.xyz --gfn 2 --input density.inp` in a scratch directory
containing the geometry and [density.inp](density.inp). [xtb-density.cub](xtb-density.cub)
is the resulting `density.cub` from xTB 6.7.1. The grid spacing is about
0.35 bohr. The figure draws the positive 0.035 e/bohr³ contour.
