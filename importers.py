#
#  Molecular Blender
#  Filename: importers.py
#  Copyright (C) 2017 Shane Parker, Joshua Szekely
#
#  This file is part of Molecular Blender.
#
#  Molecular Blender is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 3, or (at your option)
#  any later version.
#
#  Molecular Blender is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU Library General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with Molecular Blender; see COPYING.
#  If not, see <http://www.gnu.org/licenses/>.
#

from .util import stopwatch, Timer
from .periodictable import symbols
from .constants import ang2bohr, bohr2ang

import numpy as np
import re


class Reader(object):
    """Convenience class to read a line, store it, and throw if file is over"""

    def __init__(self, filename):
        self.filename = filename
        self.f = None
        self.mark = 0
        self.marks = {}

    def readline(self):
        """Ignores blanks"""
        self.mark = self.f.tell()
        out = self.f.readline()
        if out == "":
            raise EOFError
        return out

    def set_mark(self, label):
        self.marks[label] = self.mark

    def restore_mark(self, label):
        self.f.seek(self.marks[label])

    def is_marked(self, label):
        return (label in self.marks)

    def __enter__(self):
        self.f = open(self.filename, "r")
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.f.__exit__(exc_type, exc_val, exc_tb)
        self.f = None


def molecule_from_file(filename, options):
    """Tries to guess filetype based on extension, then calls the appropriate reader"""
    if ".xyz" in filename.lower():
        return molecule_from_xyz(filename, options)
    elif ".molden" in filename.lower():
        return molecule_from_molden(filename, options)
    elif ".cube" in filename.lower():
        return molecule_from_cube(filename, options)
    else:  # fall back to some sort of babel import
        return molecule_from_babel(filename, options)


def make_atom_dict(symbol, position, index, charge, gradient, hidden):
    return {"symbol": str(symbol).lower(), "position": [float(x) for x in position],
            "index": int(index), "charge": float(charge),
            "gradient": [float(x) for x in gradient], "trajectory": [], "hidden": hidden}


def make_snap_dict(position, charge, gradient):
    return {"position": [float(x) for x in position], "charge": float(charge),
            "gradient": [float(x) for x in gradient]}


@stopwatch("read xyz")
def molecule_from_xyz(filename, options):
    """Read in xyz file and return a dictionary of results"""
    out = {"atoms": []}
    ignore_H = options.get("ignore_hydrogen", False)

    with Reader(filename) as fh:
        if (options.get("plot_type", "frame") == "frame"):
            # first line contains number of atoms
            natoms = int(fh.readline().split()[0])
            # second line is a comment
            fh.readline()

            for iatom in range(natoms):
                # Expecting:
                #   <symbol> <x> <y> <z> [<charge> [<vx> <vy> <vz>]]
                # e.g.
                #   h 0.0 0.0 1.0 0.0 1.0 2.0 3.0
                tmp = fh.readline().split()
                symb = str(tmp[0]).lower()
                position = [float(x) for x in tmp[1:4]]
                charge = float(tmp[4]) if len(tmp) >= 5 else 0.0
                gradient = [float(x) for x in tmp[5:8]] if len(
                    tmp) >= 8 else [0.0 for x in range(3)]
                hidden = ignore_H and symb == "h"
                out["atoms"].append(make_atom_dict(
                    symb, position, iatom, charge, gradient, hidden))

        elif (options.get("plot_type", "frame") == "animate"):
            # first line contains number of atoms
            natoms = int(fh.readline().split()[0])
            # second line is a comment
            fh.readline()

            for iatom in range(natoms):
                tmp = fh.readline().split()
                symb = str(tmp[0]).lower()
                position = [float(x) for x in tmp[1:4]]
                charge = float(tmp[4]) if len(tmp) >= 5 else 0.0
                gradient = [float(x) for x in tmp[5:8]] if len(
                    tmp) >= 8 else [0.0] * 3
                hidden = ignore_H and symb == "h"
                new_atom = make_atom_dict(
                    symb, position, iatom, charge, gradient, hidden)
                new_atom["trajectory"] = [
                    make_snap_dict(position, charge, gradient)]

                out["atoms"].append(new_atom)

            while (True):
                try:  # gracefully exit if at the end of a file
                    line = fh.readline()
                except EOFError:
                    break

                # painfully fail if EOF is reached when it's not expected
                frame_atoms = int(line.split()[0])  # natoms
                if (frame_atoms != natoms):
                    raise Exception(
                        "All frames in trajectory must have the same number of atoms.")

                fh.readline()  # comment line
                for iatom in range(natoms):
                    tmp = fh.readline().split()
                    symb = str(tmp[0]).lower()
                    if (symb != out["atoms"][iatom]["symbol"]):
                        raise Exception(
                            "The order of the atoms must be the same for each frame in the animation.")
                    position = [float(x) for x in tmp[1:4]]
                    charge = float(tmp[4]) if len(tmp) >= 5 else 0.0
                    gradient = [float(x) for x in tmp[5:8]] if len(
                        tmp) >= 8 else [0.0 for x in range(3)]
                    out["atoms"][iatom]["trajectory"].append(
                        make_snap_dict(position, charge, gradient))

    return out


molden_re = re.compile(r"\[\s*molden\s+format\s*\]", flags=re.IGNORECASE)
atoms_re = re.compile(r"\[\s*atoms\s*\]\s*(angs|au)?", flags=re.IGNORECASE)
gto_re = re.compile(r"\[\s*gto\s*\]", flags=re.IGNORECASE)
mo_re = re.compile(r"\[\s*mo\s*\]", flags=re.IGNORECASE)
newsection_re = re.compile(r"\s*\[")


def molden_read_atoms(f, ignore_H=False):
    """reads [Atoms] section of molden files"""
    out = []
    f.restore_mark("atoms")

    line = f.readline()
    m = atoms_re.search(line)

    # conversion factor incase atomic units are input
    factor = 1.0 if m.group(1).lower() in "angs" else bohr2ang

    # advance until a line starts with [
    try:
        line = f.readline()
        m = newsection_re.match(line)
        while (not m):
            el, index, at, x, y, z = line.split()
            el = el.lower()
            index = int(index)
            pos = [factor * float(a) for a in [x, y, z]]
            charge = 0.0
            grad = [0.0, 0.0, 0.0]
            hidden = ignore_H and el == "h"
            out.append(make_atom_dict(el, pos, index, charge, grad, hidden))
            line = f.readline()
            m = newsection_re.match(line)
    except EOFError:
        pass

    return sorted(out, key=lambda atom: atom["index"])


def molden_read_gto(f):
    """reads through GTO section to collect basis information"""
    f.restore_mark("gto")

    line = f.readline()  # should just say GTO

    out = []

    try:  # allow it to quietly exit if there is no basis
        line = f.readline()
    except EOFError:
        return out

    while True:
        if newsection_re.search(line):
            break
        # specifying basis for atom iatom
        iatom = line.split()[0]

        atombasis = []

        line = f.readline()
        while True:
            shell, nprim = line.split()[0:2]

            if shell not in ["s", "sp", "p", "d", "f", "g"]:
                raise Exception("unsupported shell specified: %s" % shell)

            nprim = int(nprim)

            new_shell = {}

            new_shell["shell"] = shell
            new_shell["exponents"] = []
            new_shell["contractions"] = []

            if shell == "sp":
                raise Exception("sp shells not handled yet")

            for i in range(nprim):
                if shell != "sp":
                    exp, con = f.readline().split()
                    new_shell["exponents"].append(float(exp))
                    new_shell["contractions"].append(float(con))
                else:
                    exp, coeff1, coeff2 = f.readline().split()
                    exp, coeff1, coeff2 = float(
                        exp), float(coeff1), float(coeff2)

            atombasis.append(new_shell)

            # a blank line signals the end of this atom's basis
            try:  # allow safe deaths
                line = f.readline()
                if line.strip() == "":
                    while line.strip() == "":
                        line = f.readline()
                    break
            except EOFError:
                pass

        out.append(atombasis)

    return out


def molden_read_mo(f, nmo):
    """reads through [MO] section to collect MO coefficient information"""
    f.restore_mark("mo")

    moinf = {"sym": re.compile(r"\s*sym\s*=\s*(\S+)", flags=re.IGNORECASE),
             "ene": re.compile(r"\s*ene\s*=\s*(\S+)", flags=re.IGNORECASE),
             "spin": re.compile(r"\s*spin\s*=\s*(\S+)", flags=re.IGNORECASE),
             "occup": re.compile(r"\s*occup\s*=\s*(\S+)", flags=re.IGNORECASE)}

    re_coef = re.compile(r"\s*(\d+)\s+(\S+)", flags=re.IGNORECASE)

    line = f.readline()  # should just say MO

    out = []

    try:  # allow it to quietly exit if there is no basis
        line = f.readline()
    except EOFError:
        return out

    while True:
        if newsection_re.search(line):
            break

        # new MO
        modata = {}
        coef = []

        while True:
            found = False
            for k, srch in moinf.items():
                m = srch.match(line)
                if m:
                    modata[k] = m.group(1)
                    line = f.readline()
                    found = True
                    break
            if not found:
                break

        for i in range(nmo):
            if i != 0:
                line = f.readline()
            m = re_coef.match(line)
            if not m:
                raise Exception(
                    "Expected coefficient line, found something else")
            coef.append(float(m.group(2)))

        modata["coeff"] = coef
        out.append(modata)

        try:  # allow safe deaths
            line = f.readline()
        except EOFError:
            break

    return out


@stopwatch("read molden")
def molecule_from_molden(filename, options):
    """Read molden file to look for coordinates and optionally orbital data"""
    # molden format is not terribly efficient, but hopefully this doesn't matter
    out = {"atoms": []}

    marks = {"molden": molden_re,
             "atoms": atoms_re,
             "gto": gto_re,
             "mo": mo_re,
             "5d": re.compile(r"\[5D\]"),
             "5d10f": re.compile(r"\[5D10F\]"),
             "5d7f": re.compile(r"\[5D7F\]"),
             "9g": re.compile(r"\[9G\]")}

    # Molden defaults to Cartesian basis functions
    shelldegen = {"s": 1, "sp": 4, "p": 3, "d": 6, "f": 10, "g": 15}
    sphdegen = {"s": 1, "sp": 4, "p": 3, "d": 5, "f": 7,  "g": 9}

    ignore_H = options.get("ignore_hydrogen", False)

    with Reader(filename) as f:
        try:
            while(True):
                line = f.readline()
                for x in marks:
                    if marks[x].search(line):
                        f.set_mark(x)
        except EOFError:
            pass

        if not f.is_marked("molden"):
            raise Exception(
                "Molden file is missing [Molden Format] specification")
        if not f.is_marked("atoms"):
            raise Exception("Molden file is missing [Atoms] specification")

        # read through atoms section
        out["atoms"] = molden_read_atoms(f)

        if f.is_marked("gto"):  # read through GTO section to build basis
            make_spherical = []
            if f.is_marked("5d") or f.is_marked("5d7f"):
                make_spherical.extend(["d", "f"])
            if f.is_marked("7f"):
                make_spherical.append("f")
            if f.is_marked("5d10f"):
                make_spherical.append("d")
            if f.is_marked("9g"):
                make_spherical.append("g")

            cartesian = len(make_spherical) == 0
            if not cartesian:
                raise Exception(
                    "MolecularBlender currently can only handle cartesian d-, f-, and g-functions")
            for l in make_spherical:
                shelldegen[l] = sphdegen[l]
            out["basis"] = molden_read_gto(f)

        # read through MO section to build coefs
        if f.is_marked("mo") and "basis" in out:
            nmo = sum([sum([shelldegen[x["shell"]] for x in atom])
                       for atom in out["basis"]])
            out["mo"] = molden_read_mo(f, nmo)

    return out


@stopwatch("read cube")
def molecule_from_cube(filename, options):
    """Read a cube file including its volumetric data"""
    out = {"atoms": [], "volume": {}}

    ignore_H = options.get("ignore_hydrogen", False)

    with Reader(filename) as f:
        f.readline()  # first two lines are comments
        f.readline()  # first two lines are comments

        # 3rd line is <natoms> <origin_x> <origin_y> <origin_z>
        natoms, ox, oy, oz = f.readline().split()
        natoms = int(natoms)
        origin = np.array([ox, oy, oz], dtype=np.float32) * bohr2ang
        out["volume"]["origin"] = origin

        nres = [0, 0, 0]  # number of points in each direction
        # axes[i,:] defines the i-th axis
        axes = np.zeros([3, 3], dtype=np.float32)

        # 4th, 5th, 6th lines describe the axes that define the grid with: <npoints> <axis_x> <axis_y> <axis_z>
        for i in range(3):
            res, vx, vy, vz = f.readline().split()
            nres[i] = int(res)
            axes[i, :] = np.array([vx, vy, vz], dtype=np.float32)
            if nres[i] > 0:  # axes defined in Bohr
                axes[i, :] *= bohr2ang
            else:  # axes defined in Angstrom
                nres[i] *= -1

        out["volume"]["nres"] = nres
        out["volume"]["axes"] = axes

        # make sure axes aren't oblique
        S = np.dot(axes.T, axes)
        for x in (S[0, 1], S[0, 2], S[1, 2]):
            if abs(x) > 1e-6:
                raise Exception(
                    "Sheared axes are not supported with cube files!")

        # next natoms lines define the molecule
        # <atomic number> <charge> <x> <y> <z>
        for i in range(natoms):
            ato, chg, x, y, z = f.readline().split()
            iatom = int(ato)
            position = [bohr2ang * float(xx) for xx in [x, y, z]]
            hidden = ignore_H and symbols[iatom] == "h"
            out["atoms"].append(make_atom_dict(
                symbols[iatom], position, i, chg, [0.0, 0.0, 0.0], hidden))

        # and finally the volumetric data comes, 6 elements per line in z, y, x order
        ndata = np.prod(nres)
        data = np.zeros(ndata)
        i = 0
        while True:
            line = f.readline().split()
            ldata = np.array([float(x) for x in line])
            data[i:i + len(line)] = line
            i += len(line)
            if i == ndata:
                break

    # now should have shape [nx, ny, nz]
    out["volume"]["data"] = data.reshape(nres)

    return out


def molecule_from_babel(filename, options):
    raise Exception("Importing from babel currently disabled")
