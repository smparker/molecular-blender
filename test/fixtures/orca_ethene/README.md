# Ethene π orbital reference

ORCA 6.1.1 generated both the RHF/def2-SVP wavefunction and the orbital cube.
The geometry lies in the plane x + y = 2z. This orientation exercises diagonal
and mixed Cartesian terms in the expansion of the carbon spherical d functions.

Run these commands in a scratch directory containing `ethene.inp`.

```sh
orca ethene.inp > ethene.out
orca_2mkl ethene -molden
```

`ethene.molden` retains the header, basis and first eight occupied orbitals
from `ethene.molden.input`. The orbital under test is ORCA orbital 7, operator
0, which is orbital 8 in the add-on. The default **Auto** convention recognizes
the `orca_2mkl` title. The [figure example](../../../examples/render_orca_ethene.py)
renders this orbital from either the fixture or the full `ethene.molden.input`.

`mo8.dat` contains literal voxel samples from `pi.cube`. Positions are in
Bohr and amplitudes in Bohr⁻³ᐟ². The grid starts at (−4, −4, −4), with 33
points per axis and 0.25 Bohr spacing. For each listed position, the integer
index along each axis is `(coordinate + 4) * 4`. The z index varies fastest.
The comparison allows a single overall phase change and the rounding of
ORCA's cube output. Running the tests only requires the checked-in files.
