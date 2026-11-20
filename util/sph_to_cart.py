#!/usr/bin/env python3

from sympy import *
from sympy.printing.pycode import pycode

import argparse

ANG_LABELS = "spdfg"
MAX_L = 4

def cart_to_sph(l, m, x, y, z):
    j = (x + y - abs(m))
    if j % 2 != 0:
        return 0
    j = j/S(2)

    if x + y + z != l:
        return 0

    NN = factorial(2*x)*factorial(2*y)*factorial(2*z)*factorial(l)*factorial(l-abs(m)) / (
            factorial(2*l)*factorial(x)*factorial(y)*factorial(z)*factorial(l+abs(m)) )
    N = sqrt(NN) / ( S(2)**l * factorial(l) )

    A = S(0)
    for i in range((l-abs(m))//2 + 1):
        a = binomial(l,i) * binomial(i,j) * (-1)**i * factorial(2*l - 2*i) / (
                factorial(l - abs(m) - 2 * i) )
        B = S(0)
        for k in range(j + 1):
            B += binomial(j, k) * binomial(abs(m), x - 2*k) * S(-1)**(sign(m)*(abs(m)-x+2*S(k))/S(2))
        A += a * B

    return N * A

def cart_to_real_sph(l, m, x, y, z):
    mm = abs(m)
    if m > 0:
        cc = (cart_to_sph(l,mm,x,y,z) + cart_to_sph(l,-mm,x,y,z))/sqrt(S(2))
    elif m < 0:
        cc = (cart_to_sph(l,mm,x,y,z) - cart_to_sph(l,-mm,x,y,z))/sqrt(S(-2))
    else:
        cc = cart_to_sph(l,mm,x,y,z)
    return cc

def cart_overlap(x1,y1,z1,x2,y2,z2):
    xx = x1+x2
    yy = y1+y2
    zz = z1+z2
    if (xx % 2) != 0 or (yy % 2) != 0 or (zz % 2) != 0:
        return S(0)

    ax = factorial(xx) / factorial(xx//2)
    ay = factorial(yy) / factorial(yy//2)
    az = factorial(zz) / factorial(zz//2)
    A = ax * ay * az

    bx = factorial(2*x1) * factorial(2*x2) / (factorial(x1) * factorial(x2))
    by = factorial(2*y1) * factorial(2*y2) / (factorial(y1) * factorial(y2))
    bz = factorial(2*z1) * factorial(2*z2) / (factorial(z1) * factorial(z2))
    B = sqrt(bx * by * bz)

    return A/B

def sph_to_cart(l,m,x1,y1,z1):
    out = S(0)
    for x2 in range(l+1):
        for y2 in range(l+1):
            for z2 in range(l+1):
                if x2+y2+z2 == l:
                    out += cart_overlap(x1,y1,z1,x2,y2,z2) * conjugate(
                            cart_to_sph(l,m,x2,y2,z2))
    return out

def real_sph_to_cart(l,m,x1,y1,z1):
    out = S(0)
    for x2 in range(l+1):
        for y2 in range(l+1):
            for z2 in range(l+1):
                if x2+y2+z2 == l:
                    out += cart_overlap(x1,y1,z1,x2,y2,z2) * cart_to_real_sph(l,m,x2,y2,z2)
    return out

def gen_cart_to_real_sph_transform(l, sphorder, cartorder):
    """Return A with Y = A X, for spherical and Cartesian basis functions."""
    xyzs = [ (xyz.count("x"), xyz.count("y"), xyz.count("z")) for xyz in cartorder ]

    transform = []
    for m in sphorder:
        transform.append([cart_to_real_sph(l, m, x, y, z) for x, y, z in xyzs])

    return transform

def gen_real_sph_to_cart_transform(l, sphorder, cartorder):
    """Return overlaps <X_i|Y_m>, the projection of each Cartesian AO onto Y."""
    transform = []
    for i, xyz in enumerate(cartorder):
        x = xyz.count("x")
        y = xyz.count("y")
        z = xyz.count("z")

        transform.append( [ real_sph_to_cart(l, m, x, y, z) for m in sphorder ] )

    return transform

molden_cart_s_order = [ "" ]
molden_cart_p_order = "x y z".split()
molden_cart_d_order = "xx yy zz xy xz yz".split()
molden_cart_f_order = "xxx yyy zzz xyy xxy xxz xzz yzz yyz xyz".split()
molden_cart_g_order = "xxxx yyyy zzzz xxxy xxxz yyyx yyyz zzzx zzzy xxyy xxzz yyzz xxyz yyxz zzxy".split()

molden_cart_order = [ molden_cart_s_order, molden_cart_p_order,
        molden_cart_d_order, molden_cart_f_order, molden_cart_g_order ]

def gen_molden_sph_order(l):
    out = [0]
    for m in range(1,l+1):
        out.append(m)
        out.append(-m)
    return out

molden_sph_order = [ gen_molden_sph_order(l) for l in range(5) ]

def main():
    parser = argparse.ArgumentParser("sph_to_cart",
            description="Computes transformation coefficients between normalized " +
            "real spherical AOs and normalized cartesian AOs (CAOs). Formulas for the " +
            "coefficients are taken from Schlegel, Frisch IJQC 1995 v. 54 p. 83. " +
            "Calculations are done with sympy so the results are exact expressions, " +
            "i.e., have no precision loss"
            )

    parser.add_argument("--cao2sao", action="store_true", help="express spherical basis functions in Cartesian functions")
    parser.add_argument("--sao2cao", action="store_true", help="project Cartesian functions onto spherical functions")
    parser.add_argument("--mo-sph-to-cart", action="store_true",
                        help="generate the Molden importer's d/f/g molecular-orbital coefficient tables")
    parser.add_argument("--latex", action="store_true", help="print pretty latex transforms")

    args = parser.parse_args()

    if args.mo_sph_to_cart:
        for l in range(2, MAX_L+1):
            # If Y = A X for basis functions, a_sph.T Y = (A.T a_sph).T X.
            # The overlap projection S A.T is a different matrix.
            rows = gen_cart_to_real_sph_transform(l, molden_sph_order[l], molden_cart_order[l])
            print(f"shell: {ANG_LABELS[l]}")
            print(pycode(Matrix(rows).T.tolist()).replace("math.sqrt", "sqrt"))

    if args.latex:
        print("Cartesian to spherical transformations")
        g = Function("v")
        for l in range(MAX_L+1):
            for m in range(-l,l+1):
                term = None
                for xyz in molden_cart_order[l]:
                    x = xyz.count("x")
                    y = xyz.count("y")
                    z = xyz.count("z")
                    cc = cart_to_real_sph(l, m, x, y, z)

                    if cc != 0:
                        if term is None:
                            term = cc * g(S(x), S(y), S(z))
                        else:
                            term += cc * g(S(x), S(y), S(z))
                print("v({0},{1}) = {2}".format(l, m, latex(term)))
            print()

        print("Cartesian functions projected onto the spherical subspace")
        for l in range(MAX_L+1):
            for cart in molden_cart_order[l]:
                # looking for coefficients for v(x,y,z)
                x = cart.count("x")
                y = cart.count("y")
                z = cart.count("z")

                term = None
                for m in range(-l,l+1):
                    cc = real_sph_to_cart(l,m,x,y,z)

                    if cc != 0:
                        if term is None:
                            term = cc * g(S(l), S(m))
                        else:
                            term += cc * g(S(l), S(m))
                print("P_l v({0},{1},{2}) = {3}".format(x, y, z, latex(term)))
            print()

    if args.cao2sao:
        print("Spherical basis functions expanded in Cartesian basis functions (rows are spherical):")

        for l in range(MAX_L+1):
            print("shell: {0}".format(ANG_LABELS[l]))
            print(pycode(gen_cart_to_real_sph_transform(l, molden_sph_order[l], molden_cart_order[l])))
            print()


    if args.sao2cao:
        print("Cartesian/spherical overlap projections (rows are Cartesian):")

        for l in range(MAX_L+1):
            print("shell: {0}".format(ANG_LABELS[l]))
            print(pycode(gen_real_sph_to_cart_transform(l, molden_sph_order[l], molden_cart_order[l])))
            print()

if __name__ == "__main__":
    main()
