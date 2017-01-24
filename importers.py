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

import re

au2ang = 1.8897259885789

class Reader(object):
    """Convenience class to read a line, store it, and throw if file is over"""
    def __init__(self, filename):
        self.filename = filename
        self.f = None

    def readline(self):
        """Ignores blanks"""
        out = self.f.readline()
        if out == "":
            raise EOFError
        return out

    def __iter__(self):
        return self.f.__iter__()

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
    else: # fall back to some sort of babel import
        return molecule_from_babel(filename, options)

def make_atom_dict(symbol, position, index, charge, gradient):
    return { "symbol" : symbol.lower(), "position" : position,
        "index" : index, "charge" : charge, "gradient" : gradient, "trajectory" : [] }

def make_snap_dict(position, charge, gradient):
    return { "position" : position, "charge" : charge, "gradient" : gradient }

@stopwatch("read xyz")
def molecule_from_xyz(filename, options):
    """Read in xyz file and return a dictionary of results"""
    out = { "atoms" : [] }
    with Reader(filename) as fh:
        if (options["plot_type"] == "frame"):
            # first line contains number of atoms
            natoms = int(fh.readline().split()[0])
            # second line is a comment
            fh.readline()

            index = 0
            for line in fh:
                # Expecting:
                #   <symbol> <x> <y> <z> [<charge> [<vx> <vy> <vz>]]
                # e.g.
                #   h 0.0 0.0 1.0 0.0 1.0 2.0 3.0
                tmp = line.split()
                symb = str(tmp[0]).lower()
                position = [ float(x) for x in tmp[1:4] ]
                charge = float(tmp[4]) if len(tmp) >= 5 else 0.0
                gradient = [ float(x) for x in tmp[5:8] ] if len(tmp) >= 8 else [ 0.0 for x in range(3) ]
                out["atoms"].append(make_atom_dict(symb, position, index, charge, gradient))
                index += 1

            assert(index == natoms)
        elif (options["plot_type"] == "animate"):
            # first line contains number of atoms
            natoms = int(fh.readline().split()[0])
            # second line is a comment
            fh.readline()

            for iatom in range(natoms):
                tmp = fh.readline().split()
                symb = str(tmp[0]).lower()
                position = [ float(x) for x in tmp[1:4] ]
                charge = float(tmp[4]) if len(tmp) >= 5 else 0.0
                gradient = [ float(x) for x in tmp[5:8] ] if len(tmp) >= 8 else [ 0.0 ] * 3
                new_atom = make_atom_dict(symb, position, iatom, charge, gradient)
                new_atom["trajectory"] = [ make_snap_dict(position, charge, gradient) ]

                out["atoms"].append(new_atom)

            while (True):
                try: # gracefully exit if at the end of a file
                    line = fh.readline()
                except EOFError:
                    break

                # painfully fail if EOF is reached when it's not expected
                frame_atoms = int(line.split()[0]) # natoms
                if (frame_atoms != natoms):
                    raise Exception("All frames in trajectory must have the same number of atoms.")

                fh.readline() # comment line
                for i in range(natoms):
                    tmp = fh.readline().split()
                    symb = str(tmp[0]).lower()
                    if (symb != out["atoms"][i]["symbol"]):
                        raise Exception("The order of the atoms must be the same for each frame in the animation.")
                    position = [ float(x) for x in tmp[1:4] ]
                    charge = float(tmp[4]) if len(tmp) >= 5 else 0.0
                    gradient = [ float(x) for x in tmp[5:8] ] if len(tmp) >= 8 else [ 0.0 for x in range(3) ]
                    out["atoms"][i]["trajectory"].append(make_snap_dict(position, charge, gradient))

    return out

@stopwatch("read molden")
def molecule_from_molden(filename, options):
    """Read molden file to look for coordinates and optionally orbital data"""
    # molden format is not terribly efficient, but hopefully this doesn't matter
    out = { "atoms" : [] }

    molden_re = re.compile(r"\[\s*molden\s+format\s*\]", flags=re.IGNORECASE)
    atoms_re = re.compile(r"\[\s*atoms\s*\]\s*(angs|au)?", flags=re.IGNORECASE)
    newsection_re = re.compile(r"\s*\[")

    with Reader(filename) as f:
        try:
            line = f.readline()
            m = molden_re.search(line)
            # advance to instance of [Molden]
            while (not m):
                print(line)
                line = f.readline()
                m = molden_re.search(line)

            # advance to instances of [atom]
            line = f.readline()
            m = atoms_re.search(line)
            while (not m):
                line = f.readline()
                m = atoms_re.search(line)
        except:
            print("Exception raised before [atoms] section found")
            raise

        # conversion factor incase atomic units are input
        factor = 1.0 if m.group(1).lower() in "angs" else au2ang

        # advance until a line starts with [
        try:
            line = f.readline()
            m = newsection_re.match(line)
            while (not m):
                el, no, at, x, y, z = line.split()
                out["atoms"].append(make_atom_dict(el,
                    [ factor*float(x), factor*float(y), factor*float(z) ],
                    int(no), 0.0, [ 0.0, 0.0, 0.0] ))
                line = f.readline()
                m = newsection_re.match(line)
        except EOFError:
            pass

    # sort by index for good measure
    out["atoms"] = sorted(out["atoms"], key=lambda atom: atom["index"])

    return out

def molecule_from_babel(filename, options):
    raise Exception("Importing from babel currently disabled")
