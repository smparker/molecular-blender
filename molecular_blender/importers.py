# -*- coding: utf-8 -*-
#
#  Molecular Blender
#  Filename: importers.py
#  Copyright (C) 2017-2025 Shane Parker, Joshua Szekely
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

"""Methods to import XYZ, SDF/MOL, Molden, and cube files."""

import re
import numpy as np
import json
import os
import sys

from .options import trajectory_range
from .util import stopwatch
from .periodictable import symbols
from .constants import bohr2ang
from .transform import transform_sph_to_cart
from .turbomole_reader import read_turbomole_polarizability

class Reader(object):
    """Convenience class to read a line, store it, and throw if file is over"""

    def __init__(self, filename):
        """Build a Reader out of a filename"""
        self.filename = filename
        self.f = None
        self.marks = {}

    def readline(self):
        """Ignores blanks"""
        out = self.f.readline()
        if out == "":
            raise EOFError
        return out

    def set_mark(self, label, line):
        """Sets a mark of where label was found"""
        self.marks[label] = (self.f.tell(), line)

    def restore_mark(self, label):
        """Restores file to where mark was found"""
        mark, line = self.marks[label]
        self.f.seek(mark)
        return line

    def is_marked(self, label):
        """Returns whether mark was found"""
        return label in self.marks

    def __enter__(self):
        """Call upon entrance to a with-as clause"""
        self.f = open(self.filename, "r")
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """Call upon exit from a with-as clause"""
        self.f.__exit__(exc_type, exc_val, exc_tb)
        self.f = None


FORMAT_SUFFIXES = {
    'xyz': ('.xyz',), 'sdf': ('.sdf', '.mol', '.sd'),
    'molden': ('.molden', 'molden.input'), 'cube': ('.cub', '.cube'), 'json': ('.json',),
}


def file_format(filename):
    name = os.fspath(filename).lower()
    return next((kind for kind, suffixes in FORMAT_SUFFIXES.items() if name.endswith(suffixes)), None)


def molecule_from_file(filename, options):
    """Dispatch by filename, using the same rule as the import dialog."""
    readers = dict(xyz=molecule_from_xyz, sdf=molecule_from_sdf, molden=molecule_from_molden,
                   cube=molecule_from_cube, json=molecule_from_json)
    reader = readers.get(file_format(filename))
    if reader is not None:
        data = reader(filename, options)
        if data is None:
            raise ValueError("Missing molecules list in " + os.fspath(filename))
        return data

    raise ValueError(
        f"Unsupported molecular file format: {filename!r}. "
        "Supported formats: XYZ, SDF/MOL/SD, Molden (.molden or molden.input), "
        "cube (.cube or .cub), and JSON."
    )


def make_atom_dict(symbol, position, index, charge, gradient, hidden):
    """Build dictionary for atom"""
    return {"symbol": str(symbol).lower(), "position": [float(x) for x in position],
            "index": int(index), "charge": float(charge),
            "gradient": [float(x) for x in gradient], "trajectory": [], "hidden": hidden}


def make_snap_dict(position, charge, gradient):
    """Build dictionary for snapshot"""
    return {"position": [float(x) for x in position], "charge": float(charge),
            "gradient": [float(x) for x in gradient]}


@stopwatch("read sdf")
def molecule_from_sdf(filename, options):
    """Read a V2000 MOL or one SDF record, preserving explicit bond orders.

    ``sdf_record`` selects a zero-based record (default: first). Coordinates are
    in angstroms. Formal charges are kept separately from partial charges.
    Aromatic/query bond codes and V3000 are rejected rather than approximated.
    Stereo, isotope, radical, S-group and SD data annotations are not interpreted.
    """
    record = options.get("sdf_record", 0)
    if not isinstance(record, int) or isinstance(record, bool) or record < 0:
        raise ValueError("sdf_record must be a nonnegative integer")
    lines = []
    with open(filename) as source:
        current = 0
        for line in source:
            if line.strip() == "$$$$":
                if current == record:
                    break
                current += 1
                lines = []
            elif current == record:
                lines.append(line.rstrip('\r\n'))
    if not lines or not any(line.strip() for line in lines):
        raise ValueError("SDF record {} not found".format(record))
    if len(lines) < 4 or not lines[3].rstrip().endswith("V2000"):
        raise ValueError("Only SDF/MOL V2000 is supported (expected V2000 counts line)")
    try:
        natoms, nbonds = int(lines[3][:3]), int(lines[3][3:6])
    except ValueError as exc:
        raise ValueError("Invalid V2000 atom/bond counts") from exc
    if natoms <= 0 or nbonds < 0:
        raise ValueError("V2000 requires positive atom count and nonnegative bond count")
    end = 4 + natoms + nbonds
    if len(lines) <= end:
        raise ValueError("Truncated V2000 atom/bond blocks or missing M  END")
    out = {"atoms": [], "bonds": [], "title": lines[0]}
    ignore_h = options.get("ignore_hydrogen", False)
    charge_codes = {0: 0, 1: 3, 2: 2, 3: 1, 4: 0, 5: -1, 6: -2, 7: -3}
    for index, line in enumerate(lines[4:4 + natoms]):
        try:
            position = [float(line[i:i + 10]) for i in (0, 10, 20)]
            symbol = line[31:34].strip().lower()
            code = int(line[36:39].strip() or "0")
            formal_charge = charge_codes[code]
            if symbol not in symbols or not np.all(np.isfinite(position)):
                raise ValueError("Invalid element or non-finite coordinates")
        except (ValueError, KeyError) as exc:
            raise ValueError("Invalid V2000 atom {}".format(index + 1)) from exc
        atom = make_atom_dict(symbol, position, index, 0, [0, 0, 0], ignore_h and symbol == "h")
        atom["formal_charge"] = formal_charge
        out["atoms"].append(atom)
    seen = set()
    for number, line in enumerate(lines[4 + natoms:end], 1):
        try:
            i, j, order = (int(line[k:k + 3]) for k in (0, 3, 6))
        except ValueError as exc:
            raise ValueError("Invalid V2000 bond {}".format(number)) from exc
        if not (1 <= i <= natoms and 1 <= j <= natoms) or i == j:
            raise ValueError("Invalid atom reference in V2000 bond {}".format(number))
        pair = tuple(sorted((i, j)))
        if pair in seen:
            raise ValueError("Duplicate V2000 bond {}".format(number))
        if order not in (1, 2, 3):
            raise ValueError("Unsupported V2000 bond type {}: use explicit single/double/triple bonds".format(order))
        seen.add(pair)
        out["bonds"].append({"atoms": [i - 1, j - 1], "order": order})
    found_end = False
    has_charges = False
    for line in lines[end:]:
        if line.strip() == "M  END":
            found_end = True
            break
        if line.startswith("M  CHG"):
            try:
                fields = line.split()
                count = int(fields[2])
                if count < 1 or len(fields) != 3 + 2 * count:
                    raise ValueError("Invalid charge count")
                # M CHG supersedes *all* atom-block charge codes.
                if not has_charges:
                    for atom in out["atoms"]:
                        atom["formal_charge"] = 0
                    has_charges = True
                for k in range(count):
                    index, charge = int(fields[3 + 2*k]), int(fields[4 + 2*k])
                    if not 1 <= index <= natoms:
                        raise ValueError("Invalid charged atom reference")
                    out["atoms"][index - 1]["formal_charge"] = charge
            except (ValueError, IndexError) as exc:
                raise ValueError("Invalid V2000 M  CHG line") from exc
    if not found_end:
        raise ValueError("Missing V2000 M  END")
    options["plot_type"] = "frame"
    return out


@stopwatch("read xyz")
def molecule_from_xyz(filename, options):
    """Stream XYZ snapshots, retaining atom data only for selected source frames."""
    start, end, step = trajectory_range(options)
    animate = options.get("plot_type") in ("animate", "auto")
    atoms, frames = [], []
    count = None
    with Reader(filename) as reader:
        index = 0
        while end < 0 or index <= end:
            line = reader.f.readline()
            if not line:
                break
            natoms = int(line.split()[0])
            if count is None:
                count = natoms
            if natoms != count or natoms < 1:
                raise ValueError("All frames in trajectory must have the same number of atoms")
            reader.readline()
            keep = index >= start and (index - start) % step == 0
            for i in range(natoms):
                line = reader.readline()
                if not keep:
                    continue
                fields = line.split()
                if len(fields) < 4:
                    raise ValueError(f"XYZ frame {index}, atom {i} needs an element and three coordinates")
                symbol = fields[0].lower()
                position = list(map(float, fields[1:4]))
                charge = float(fields[4]) if len(fields) >= 5 else 0.
                gradient = list(map(float, fields[5:8])) if len(fields) >= 8 else [0.] * 3
                if not frames:
                    atoms.append(make_atom_dict(symbol, position, i, charge, gradient,
                                                options.get("ignore_hydrogen", False) and symbol == "h"))
                elif symbol != atoms[i]["symbol"]:
                    raise ValueError("The order of the atoms must be the same for each frame")
                if animate:
                    atoms[i]["trajectory"].append(make_snap_dict(position, charge, gradient))
            if keep:
                frames.append(index)
                if not animate:
                    break
            index += 1
    if not frames:
        raise ValueError("Frame range selects no XYZ snapshots")
    if animate and len(frames) > 1:
        options["plot_type"] = "animate"
    return {"atoms": atoms, "frame_numbers": frames}


ANY_RE = re.compile(r"\[[a-zA-Z0-9 ]*\]")
MOLDEN_RE = re.compile(r"\[\s*molden\s+format\s*\]", flags=re.IGNORECASE)
ATOMS_RE = re.compile(r"\[\s*atoms\s*\]\s*(?:\(?(angs|au)\)?)?", flags=re.IGNORECASE)
GTO_RE = re.compile(r"\[\s*gto\s*\]", flags=re.IGNORECASE)
MO_RE = re.compile(r"\[\s*mo\s*\]", flags=re.IGNORECASE)
NEWSECTION_RE = re.compile(r"\s*\[")


def molden_read_atoms(f, ignore_h=False, animate=False):
    """reads [Atoms] section of molden files"""
    out = []
    line = f.restore_mark("atoms")
    m = ATOMS_RE.search(line)

    # conversion factor incase atomic units are input
    unit = "au" if m.group(1) and m.group(1).lower() in "au" else "angs"
    factor = 1.0 if unit == "angs" else bohr2ang

    # advance until a line starts with [
    try:
        line = f.readline()
        m = NEWSECTION_RE.match(line)
        while not m:
            if not line.strip():
                line = f.readline()
                m = NEWSECTION_RE.match(line)
                continue
            element, index, _atom, x, y, z = line.split()
            element = element.lower()
            index = int(index)
            pos = [factor * molden_float(a) for a in [x, y, z]]
            charge = 0.0
            grad = [0.0, 0.0, 0.0]
            hidden = ignore_h and element == "h"
            new_atom = make_atom_dict(element, pos, index, charge, grad, hidden)
            if animate:
                new_atom["trajectory"] = [make_snap_dict(pos, charge, grad)]
            out.append(new_atom)
            line = f.readline()
            m = NEWSECTION_RE.match(line)
    except EOFError:
        pass

    return sorted(out, key=lambda atom: atom["index"])


def molden_float(value):
    """Read finite Molden numbers, including Fortran D/d exponents."""
    result = float(value.replace("D", "E").replace("d", "e"))
    if not np.isfinite(result):
        raise ValueError("Non-finite Molden number: {}".format(value))
    return result


def molden_read_gto(f, atom_ids=None):
    """Read GTO blocks by atom reference, expanding sp into s and p shells."""
    f.restore_mark("gto")
    blocks = {}
    current = None
    while True:
        try:
            line = f.readline()
        except EOFError:
            break
        if NEWSECTION_RE.match(line):
            break
        fields = line.split()
        if not fields:
            continue
        if fields[0].isdigit():
            reference = int(fields[0])
            if reference < 1 or reference in blocks:
                raise ValueError("Invalid or duplicate GTO atom reference: {}".format(reference))
            if atom_ids is not None and reference not in atom_ids:
                raise ValueError("Unknown GTO atom reference: {}".format(reference))
            current = blocks[reference] = []
            continue
        shell = fields[0].lower()
        if current is None or shell not in ("s", "sp", "p", "d", "f", "g"):
            raise ValueError("Invalid GTO shell header: {}".format(line.strip()))
        try:
            nprim = int(fields[1])
            if nprim < 1:
                raise ValueError("Primitive count must be positive")
        except (IndexError, ValueError) as exc:
            raise ValueError("Invalid GTO primitive count: {}".format(line.strip())) from exc
        exponents, coefficients, pcoefficients = [], [], []
        for _ in range(nprim):
            try:
                row = f.readline().split()
                if len(row) != (3 if shell == "sp" else 2):
                    raise ValueError("Incorrect primitive column count")
                numbers = [molden_float(x) for x in row]
                if numbers[0] <= 0:
                    raise ValueError("Primitive exponents must be positive")
            except (EOFError, ValueError) as exc:
                raise ValueError("Invalid or truncated GTO {} primitive block".format(shell)) from exc
            exponents.append(numbers[0])
            coefficients.append(numbers[1])
            if shell == "sp":
                pcoefficients.append(numbers[2])
        current.append({"shell": "s" if shell == "sp" else shell,
                        "exponents": exponents, "contractions": coefficients})
        if shell == "sp":
            current.append({"shell": "p", "exponents": exponents.copy(),
                            "contractions": pcoefficients})
    references = sorted(blocks) if atom_ids is None else atom_ids
    return [blocks.get(reference, []) for reference in references]


def molden_read_mo(f, nao):
    """Read complete indexed MO coefficient blocks without crossing sections."""
    if nao < 1:
        raise ValueError('Molden MO section requires at least one basis function')
    f.restore_mark("mo")
    metadata = re.compile(r"\s*(Sym|Ene|Spin|Occup)\s*=\s*(\S+)", re.IGNORECASE)
    out, modata, coefficients = [], {}, []

    def finish():
        if modata or coefficients:
            if len(coefficients) != nao:
                raise ValueError("Molden MO expected {} coefficients, found {}".format(
                    nao, len(coefficients)))
            out.append(dict(modata, coeff=coefficients.copy()))

    while True:
        try:
            line = f.readline()
        except EOFError:
            break
        if NEWSECTION_RE.match(line):
            break
        if not line.strip():
            continue
        match = metadata.fullmatch(line.rstrip())
        if match:
            key, value = match.group(1).lower(), match.group(2)
            if coefficients or key in modata:
                finish()
                modata, coefficients = {}, []
            if key in ('ene', 'occup'):
                value = str(molden_float(value))
            modata[key] = value
        else:
            fields = line.split()
            try:
                if len(fields) != 2 or int(fields[0]) != len(coefficients) + 1:
                    raise ValueError("Unexpected coefficient index")
                if len(coefficients) >= nao:
                    raise ValueError("Too many coefficients")
                coefficients.append(molden_float(fields[1]))
            except ValueError as exc:
                raise ValueError("Invalid Molden MO coefficient line: {}".format(line.strip())) from exc
    finish()
    return out


@stopwatch("read molden")
def molecule_from_molden(filename, _options):
    """Read molden file to look for coordinates and optionally orbital data"""
    # molden format is not terribly efficient, but hopefully this doesn't matter
    ignore_h = _options.get("ignore_hydrogen", False)
    animate = _options.get("plot_type", "") == "animate"
    out = {"atoms": []}
    out["denorm"] = not _options.get("orca_fix", False)

    marks = {"molden": MOLDEN_RE,
             "atoms": ATOMS_RE,
             "gto": GTO_RE,
             "mo": MO_RE,
             "5d": re.compile(r"\[5D\]", flags=re.IGNORECASE),
             "5d10f": re.compile(r"\[5D10F\]", flags=re.IGNORECASE),
             "5d7f": re.compile(r"\[5D7F\]", flags=re.IGNORECASE),
             "7f": re.compile(r"\[7F\]", flags=re.IGNORECASE),
             "9g": re.compile(r"\[9G\]", flags=re.IGNORECASE)}

    # Molden defaults to Cartesian basis functions
    shelldegen = {"s": 1, "sp": 4, "p": 3, "d": 6, "f": 10, "g": 15}
    cartdegen = {"s": 1, "sp": 4, "p": 3, "d": 6, "f": 10, "g": 15}
    sphdegen = {"s": 1, "sp": 4, "p": 3, "d": 5, "f": 7, "g": 9}

    with Reader(filename) as f:
        try:
            while True:
                line = f.readline()
                if ANY_RE.search(line):
                    for x in marks:
                        if marks[x].search(line):
                            f.set_mark(x, line)
        except EOFError:
            pass

        if not f.is_marked("molden"):
            raise Exception(
                "Molden file is missing [Molden Format] specification")
        if not f.is_marked("atoms"):
            raise Exception("Molden file is missing [Atoms] specification")

        # read through atoms section
        out["atoms"] = molden_read_atoms(f, ignore_h, animate=animate)

        if f.is_marked("gto"):  # read through GTO section to build basis
            out["basis"] = molden_read_gto(f, [atom["index"] for atom in out["atoms"]])

        # read through MO section to build coefs
        if f.is_marked("mo") and "basis" in out:
            spherical = []
            if f.is_marked("5d") or f.is_marked("5d7f"):
                spherical.extend(["d", "f"])
            if f.is_marked("7f"):
                spherical.append("f")
            if f.is_marked("5d10f"):
                spherical.append("d")
            if f.is_marked("9g"):
                spherical.append("g")

            for l in spherical:
                shelldegen[l] = sphdegen[l]
            nao = sum([sum([shelldegen[x["shell"]] for x in atom])
                       for atom in out["basis"]])
            molden_mo = molden_read_mo(f, nao)

            if spherical and molden_mo:
                molden_mo = transform_sph_to_cart(molden_mo, spherical, out["basis"])

            out["mo"] = molden_mo

    return out


def _cube_origin_header(f):
    header = f.readline().split()
    if len(header) not in (4, 5):
        raise ValueError('Cube atom/origin header requires four or five values')
    nvalues = int(header[4]) if len(header) == 5 else 1
    if nvalues < 1:
        raise ValueError('Cube dataset count must be positive')
    if nvalues > sys.maxsize:
        raise ValueError('Cube dataset count exceeds the supported range')
    return int(header[0]), header[1:4], nvalues


def _cube_labels(f, natoms, nvalues):
    if natoms >= 0:
        return range(nvalues)
    tokens = f.readline().split()
    if not tokens:
        raise ValueError('Missing cube dataset IDs')
    count = int(tokens[0])
    if count < 1 or (nvalues != 1 and nvalues != count):
        raise ValueError('Invalid cube dataset count')
    labels = [int(value) for value in tokens[1:]]
    while len(labels) < count:
        labels.extend(int(value) for value in f.readline().split())
    if len(labels) != count:
        raise ValueError('Incorrect number of cube dataset IDs')
    return labels


def cube_dataset_labels(filename):
    """Read dataset identifiers without loading or allocating the voxel grid."""
    with Reader(filename) as f:
        f.readline()
        f.readline()
        natoms, _, nvalues = _cube_origin_header(f)
        if natoms < 0:
            for _ in range(3 + abs(natoms)):
                f.readline()
        return _cube_labels(f, natoms, nvalues)


@stopwatch("read cube")
def molecule_from_cube(filename, options):
    """Read a cube file including its volumetric data"""
    out = {"atoms": [], "volume": {}}

    ignore_h = options.get("ignore_hydrogen", False)
    animate = options.get("plot_type", "") in [ "animate", "auto" ]

    with Reader(filename) as f:
        out["volume"]["comments"] = [f.readline().rstrip("\r\n"), f.readline().rstrip("\r\n")]

        # 3rd line is <natoms> <origin_x> <origin_y> <origin_z>
        natoms, origin, nvalues = _cube_origin_header(f)

        nres = [0, 0, 0]  # number of points in each direction
        # axes[i,:] defines the i-th axis
        axes = np.zeros([3, 3], dtype=np.float32)

        # 4th, 5th, 6th lines describe the axes: <npoints> <axis_x> <axis_y> <axis_z>
        for i in range(3):
            res, vx, vy, vz = f.readline().split()
            nres[i] = int(res)
            axes[i, :] = np.array([vx, vy, vz], dtype=np.float32)
        if all(n > 0 for n in nres):
            length_unit = bohr2ang
        elif all(n < 0 for n in nres):
            length_unit = 1.0
        else:
            raise ValueError("Cube voxel counts must be nonzero and have the same sign")
        nres = [abs(n) for n in nres]
        axes *= length_unit
        out["volume"]["origin"] = np.array(origin, dtype=np.float32) * length_unit
        if not np.isfinite(out["volume"]["origin"]).all():
            raise ValueError("Cube origin must be finite")
        out["volume"]["nres"] = nres
        out["volume"]["axes"] = axes

        # make sure axes aren't oblique
        lengths = np.linalg.norm(axes, axis=1)
        if not np.all(np.isfinite(lengths) & (lengths > 0)):
            raise ValueError("Cube grid axes must be finite and nonzero")
        directions = axes / lengths[:, None]
        if not np.allclose(directions @ directions.T, np.eye(3), atol=1e-5):
            raise ValueError("Sheared axes are not supported with cube files")

        # next natoms lines define the molecule
        # <atomic number> <charge> <x> <y> <z>
        for i in range(abs(natoms)):
            ato, chg, x, y, z = f.readline().split()
            iatom = int(ato)
            position = [length_unit * float(xx) for xx in [x, y, z]]
            hidden = ignore_h and symbols[iatom] == "h"
            new_atom = make_atom_dict(symbols[iatom], position, i, chg,
                    [0.0, 0.0, 0.0], hidden)

            if animate:
                new_atom["trajectory"] = [make_snap_dict(position, 0.0, [0.0, 0.0, 0.0])]

            out["atoms"].append(new_atom)

        # optionally check for DSET_IDs, basically labels of different data sets
        volume_labels = list(_cube_labels(f, natoms, nvalues))

        # and finally the volumetric data comes, 6 elements per line in z, y, x order
        nvol = len(volume_labels)
        data_layout = nres + [ nvol ]
        ndata = np.prod(data_layout)
        data = np.empty(ndata, dtype=np.float32)
        count = 0
        for line in f.f:
            tokens = line.split()
            end = count + len(tokens)
            if end > ndata:
                raise ValueError('Cube expected {} field values, found at least {}'.format(ndata, end))
            data[count:end] = tokens
            count = end
        if count != ndata:
            raise ValueError('Cube expected {} field values, found {}'.format(ndata, count))
        if not np.all(np.isfinite(data)):
            raise ValueError('Cube field values must be finite')

    datasets = data.reshape(data_layout)
    selected = options.get('cube_dataset', 0)
    if isinstance(selected, bool) or not isinstance(selected, int) or not 0 <= selected < nvol:
        raise ValueError('cube_dataset must be a zero-based index in [0, {})'.format(nvol))
    out['volume']['datasets'] = datasets
    out['volume']['labels'] = volume_labels
    out['volume']['dataset_index'] = selected
    # Keep the existing three-dimensional data entry for callers.
    out['volume']['data'] = np.ascontiguousarray(datasets[..., selected])

    return out


def _check_manifest_frame(first, current):
    atoms = lambda data: [(a["index"], a["symbol"]) for a in data["atoms"]]
    if atoms(first) != atoms(current):
        raise ValueError("Atom count, indices and element order must match the first frame")
    if any((key in first) != (key in current) for key in ("bonds", "basis", "mo", "volume", "polar")):
        raise ValueError("Topology and field types must match the first frame")
    if "bonds" in first:
        def topology(data):
            return sorted((tuple(sorted(b["atoms"])), b["order"]) for b in data["bonds"])
        if (topology(first) != topology(current) or
                [a["formal_charge"] for a in first["atoms"]] != [a["formal_charge"] for a in current["atoms"]]):
            raise ValueError("SDF bond orders, endpoints and formal charges must match the first frame")
    if first.get("basis") != current.get("basis"):
        raise ValueError("Molden basis must match the first frame")
    if "mo" in first:
        def channels(data):
            return [(len(m["coeff"]), (m.get("spin") or "").lower()) for m in data["mo"]]
        if channels(first) != channels(current) or len(first.get("polar", [])) != len(current.get("polar", [])):
            raise ValueError("Molden orbital counts, spin channels and field counts must match the first frame")
    if "volume" in first and first["volume"].get("units") != current["volume"].get("units"):
        raise ValueError("Cube scalar units must match the first frame")
    if "volume" in first and list(first["volume"]["labels"]) != list(current["volume"]["labels"]):
        raise ValueError("Cube dataset IDs and order must match the first frame")


def _linear_combination(entry, directory, options):
    terms = entry["linear_combination"]
    units = entry.get("units")
    if not isinstance(terms, list) or not terms:
        raise ValueError("linear_combination requires a nonempty list of cube terms")
    if not isinstance(units, str) or not units.strip():
        raise ValueError("Declare the common scalar units of the cube terms with units")
    if any(key in entry for key in ("filename", "polar", "sdf_record")):
        raise ValueError("A linear_combination replaces the filename in a manifest entry")
    result = None
    for index, term in enumerate(terms, 1):
        name = term.get("filename") if isinstance(term, dict) else None
        if not isinstance(name, str) or file_format(name) != "cube":
            raise ValueError(f"Term {index} needs a cube filename")
        weight = term.get("coefficient")
        if type(weight) not in (int, float) or not np.isfinite(float(weight)):
            raise ValueError(f"Term {index} ({name}) needs a finite coefficient")
        if term.get("units", units) != units:
            raise ValueError(f"Term {index} ({name}) has different scalar units")
        try:
            cube = molecule_from_cube(os.path.join(directory, name),
                                      dict(options, cube_dataset=term.get("cube_dataset", 0)))
        except (ValueError, OSError, EOFError) as exc:
            raise ValueError(f"Term {index} ({name}): {exc}") from exc
        grid = cube["volume"]
        if result is None:
            result = cube
            template = grid
            total = np.zeros(grid["nres"], dtype=np.float64)
        elif (grid["nres"] != template["nres"]
              or not np.allclose(grid["origin"], template["origin"], atol=1e-5, rtol=0)
              or not np.allclose((grid["axes"] - template["axes"]) * (np.array(grid["nres"]) - 1)[:, None],
                                 0, atol=1e-5, rtol=0)):
            raise ValueError(f"Term {index} ({name}) has a different grid shape, origin or axes")
        total += np.multiply(grid["data"], weight, dtype=np.float64)
        # Only geometry and grid metadata are needed from the first term.
        del grid["data"], grid["datasets"]
    if not np.all(np.isfinite(total)) or np.max(np.abs(total)) > np.finfo(np.float32).max:
        raise ValueError("Linear combination exceeds the finite cube value range")
    values = total.astype(np.float32)
    template.update(data=values, datasets=values[..., None], labels=[0], dataset_index=0,
                    field_type="signed", units=units, linear_combination=terms)
    return result


def molecule_from_json(filename, options):
    """Read one structure per manifest entry, with zero-based timeline indices."""
    with open(filename) as source:
        data = json.load(source)
    if not isinstance(data, dict):
        raise ValueError("Expected a JSON object in " + os.fspath(filename))
    if "molecules" not in data:
        return None  # Historical direct reader API.
    entries = data["molecules"]
    if not isinstance(entries, list) or not entries:
        raise ValueError("Expected a nonempty molecules list in " + os.fspath(filename))

    start, end, step = trajectory_range(options)
    previous = -1
    frames = []
    merged = None
    for i, entry in enumerate(entries):
        label = "{}: entry {}".format(filename, i + 1)
        try:
            if not isinstance(entry, dict):
                raise ValueError("Expected a manifest entry")
            frame = entry.get("frame", previous + 1)
            if type(frame) is not int or frame <= previous:
                raise ValueError("Frame indices must be strictly increasing nonnegative integers")
            previous = frame
            # A single XYZ reference is also used to wrap an entire XYZ trajectory.
            wrapped_xyz = (len(entries) == 1 and isinstance(entry.get("filename"), str)
                           and file_format(entry["filename"]) == "xyz" and "frame" not in entry)
            if not wrapped_xyz and (i < start or (end >= 0 and i > end) or (i - start) % step):
                continue
            frames.append(frame)
            child_options = dict(options)
            if len(entries) > 1 or "frame" in entry:
                child_options.update(plot_type="frame", frame_start=0, frame_end=-1, frame_step=1)
            directory = os.path.dirname(filename)
            if "linear_combination" in entry:
                current = _linear_combination(entry, directory, child_options)
            else:
                if not isinstance(entry.get("filename"), str) or not entry["filename"].strip():
                    raise ValueError("Expected a filename")
                path = os.path.join(directory, entry["filename"])
                label += " (" + entry["filename"] + ")"
                if file_format(path) == "json":
                    raise ValueError("A trajectory cannot reference another json manifest")
                if "sdf_record" in entry:
                    if file_format(path) != "sdf":
                        raise ValueError("sdf_record requires an SDF/MOL file")
                    child_options["sdf_record"] = entry["sdf_record"]
                current = molecule_from_file(path, child_options)
                if "polar" in entry:
                    current["polar"] = read_turbomole_polarizability(os.path.join(directory, entry["polar"]))
            if len(entries) == 1 and "frame" not in entry:
                if "plot_type" in child_options:
                    options["plot_type"] = child_options["plot_type"]
                return current

            snapshots = [make_snap_dict(a["position"], a["charge"], a["gradient"]) for a in current["atoms"]]
            if merged is None:
                merged = current
                for atom, snap in zip(merged["atoms"], snapshots):
                    atom["trajectory"] = [snap]
                continue
            _check_manifest_frame(merged, current)
            for atom, snap in zip(merged["atoms"], snapshots):
                atom["trajectory"].append(snap)
            if "volume" in current:
                merged.setdefault("volume_trajectory", []).append(current["volume"])
            if "mo" in current:
                merged.setdefault("orbital_trajectory", []).append(current)
        except Exception as exc:
            raise ValueError("{}: {}".format(label, exc or "invalid frame")) from exc

    if merged is None:
        raise ValueError("Frame range selects no JSON entries in " + os.fspath(filename))
    merged["frame_numbers"] = frames
    options["plot_type"] = "animate"
    return merged
