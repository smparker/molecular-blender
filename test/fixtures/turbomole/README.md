# Paired TURBOMOLE Molden exports

`molden.input` and `molden_std.input` are the output of `tm2molden batch`
and `tm2molden batch norm` on the same one-center test input, using TURBOMOLE
8.1 (revision 44af7e6504). Each file retains the occupied orbital and first
virtual orbital. The test compares every exported Cartesian coefficient.

The input uses He at the origin in C1 symmetry, with one spherical shell
of each type from s through g. Every shell has one primitive of exponent
and coefficient 1. The 25 input orbital rows form an orthonormal cosine
matrix. The occupied row has coefficient 0.2 throughout. The first virtual
row has coefficient `sqrt(2/25) * cos(pi*(j+0.5)/25)` for columns 0–24.
These prescribed orbitals exercise all angular components in both exports.
