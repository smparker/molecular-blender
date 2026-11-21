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
import math
from itertools import islice
import numpy as np
import json
import os
import sys

from .diagnostics import InputError
from .options import molden_convention, trajectory_range
from .esp import add_resampled_cube
from .util import stopwatch
from .periodictable import elements, symbols
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
        self.line_number = 0
        self.section = None
        self.context = {}

    def readline(self):
        """Read the next physical line, counting the expected line at EOF."""
        out = self.f.readline()
        self.line_number += 1
        if out == "":
            raise EOFError
        return out

    def error(self, section, message, **context):
        return InputError(self.filename, section, message, line=self.line_number,
                          **dict(self.context, **context))

    def set_mark(self, label, line):
        """Sets a mark of where label was found"""
        self.marks[label] = (self.f.tell(), line, self.line_number)

    def restore_mark(self, label):
        """Restores file to where mark was found"""
        mark, line, self.line_number = self.marks[label]
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
        if self.section and isinstance(exc_val, (ValueError, EOFError)) and not isinstance(exc_val, InputError):
            raise self.error(self.section, str(exc_val) or "Unexpected end of file") from exc_val


FORMAT_SUFFIXES = {
    'xyz': ('.xyz',), 'sdf': ('.sdf', '.mol', '.sd'),
    'molden': ('.molden', 'molden.input', 'molden_std.input'), 'cube': ('.cub', '.cube'), 'json': ('.json',),
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
        try:
            data = reader(filename, options)
        except InputError:
            raise
        except (OSError, ValueError, EOFError) as exc:
            raise InputError(filename, file_format(filename), str(exc) or "Unexpected end of file") from exc
        if data is None:
            raise InputError(filename, "JSON", "Missing molecules list")
        data.setdefault("source", os.fspath(filename))
        return data

    raise ValueError(
        f"Unsupported molecular file format: {filename!r}. "
        "Supported formats: XYZ, SDF/MOL/SD, Molden (.molden, molden.input or molden_std.input), "
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
    with Reader(filename) as reader:
        label = f"SDF record {record}"
        reader.context = dict(record=record)
        reader.section = label
        for _ in range(record):
            while reader.readline().strip() != "$$$$":
                pass

        def line(section, **context):
            reader.context = dict(record=record, **context)
            reader.section = f"{label}, {section}"
            text = reader.readline().rstrip('\r\n')
            if text.strip() == "$$$$":
                raise ValueError("Unexpected end of record")
            return text

        title = line("title")
        line("header")
        line("header")
        counts = line("V2000 counts")
        if not counts.rstrip().endswith("V2000"):
            raise ValueError("Only SDF/MOL V2000 is supported (expected V2000 counts line)")
        natoms, nbonds = int(counts[:3]), int(counts[3:6])
        if natoms <= 0 or nbonds < 0:
            raise ValueError("V2000 requires positive atom count and nonnegative bond count")

        out = {"atoms": [], "bonds": [], "title": title}
        ignore_h = options.get("ignore_hydrogen", False)
        charge_codes = {0: 0, 1: 3, 2: 2, 3: 1, 4: 0, 5: -1, 6: -2, 7: -3}
        for index in range(natoms):
            row = line(f"atom {index + 1}", atom=index+1)
            position = [_finite_float(row[i:i + 10], axis) for i, axis in zip((0, 10, 20), 'xyz')]
            symbol = row[31:34].strip().lower()
            if symbol not in symbols:
                raise ValueError(f"Unknown element {symbol!r}")
            code = int(row[36:39].strip() or "0")
            if code not in charge_codes:
                raise ValueError(f"Invalid atom charge code {code}")
            atom = make_atom_dict(symbol, position, index, 0, [0, 0, 0], ignore_h and symbol == "h")
            atom["formal_charge"] = charge_codes[code]
            out["atoms"].append(atom)

        seen = set()
        for number in range(1, nbonds + 1):
            row = line(f"bond {number}", bond=number)
            i, j, order = (int(row[k:k + 3]) for k in (0, 3, 6))
            if not (1 <= i <= natoms and 1 <= j <= natoms) or i == j:
                raise ValueError(f"Invalid atom references {i}, {j}")
            pair = tuple(sorted((i, j)))
            if pair in seen:
                raise ValueError(f"Duplicate V2000 bond between {i} and {j}")
            if order not in (1, 2, 3):
                raise ValueError(f"Unsupported V2000 bond type {order}: use explicit single/double/triple bonds")
            seen.add(pair)
            out["bonds"].append({"atoms": [i - 1, j - 1], "order": order})

        has_charges = False
        while True:
            row = line("properties, expected M  END")
            if row.strip() == "M  END":
                break
            if row.startswith("M  CHG"):
                reader.section = f"{label}, M  CHG"
                fields = row.split()
                if len(fields) < 3:
                    raise ValueError("Missing charge count")
                count = int(fields[2])
                if count < 1 or len(fields) != 3 + 2 * count:
                    raise ValueError("Invalid charge count")
                # M CHG supersedes all atom-block charge codes.
                if not has_charges:
                    for atom in out["atoms"]:
                        atom["formal_charge"] = 0
                    has_charges = True
                for k in range(count):
                    index, charge = int(fields[3 + 2*k]), int(fields[4 + 2*k])
                    if not 1 <= index <= natoms:
                        raise ValueError(f"Invalid charged atom reference {index}")
                    out["atoms"][index - 1]["formal_charge"] = charge
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
            reader.context = dict(frame=index)
            reader.section = f"XYZ frame {index}, atom count"
            try:
                line = reader.readline()
            except EOFError:
                break
            fields = line.split()
            if not fields:
                raise ValueError("Missing atom count")
            natoms = int(fields[0])
            if count is None:
                count = natoms
            if natoms != count or natoms < 1:
                raise ValueError("All frames in trajectory must have the same number of atoms")
            reader.section = f"XYZ frame {index}, comment"
            reader.readline()
            keep = index >= start and (index - start) % step == 0
            for i in range(natoms):
                reader.context["atom"] = i
                reader.section = f"XYZ frame {index}, atom {i}"
                line = reader.readline()
                if not keep:
                    continue
                fields = line.split()
                if len(fields) < 4:
                    raise ValueError("Expected an element and three coordinates")
                symbol = fields[0].lower()
                if symbol not in elements:
                    raise ValueError(f"Unknown element {fields[0]!r}")
                position = [_finite_float(token, axis) for token, axis in zip(fields[1:4], 'xyz')]
                charge = _finite_float(fields[4], 'charge') if len(fields) >= 5 else 0.
                gradient = ([_finite_float(token, 'gradient ' + axis) for token, axis in zip(fields[5:8], 'xyz')]
                            if len(fields) >= 8 else [0.] * 3)
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
        raise InputError(filename, "XYZ", "Frame range selects no XYZ snapshots")
    if animate and len(frames) > 1:
        options["plot_type"] = "animate"
    return {"atoms": atoms, "frame_numbers": frames}


def _finite_float(token, field):
    try:
        value = float(token)
    except ValueError:
        value = math.nan
    if not math.isfinite(value):
        raise ValueError(f"Expected a finite {field} value, found {token!r}")
    return value


MOLDEN_RE = re.compile(r"\[\s*molden\s+format\s*\]", flags=re.IGNORECASE)
ATOMS_RE = re.compile(r"\[\s*atoms\s*\]\s*(?:\(?(angs|au)\)?)?", flags=re.IGNORECASE)
GTO_RE = re.compile(r"\[\s*gto\s*\]", flags=re.IGNORECASE)
MO_RE = re.compile(r"\[\s*mo\s*\]", flags=re.IGNORECASE)
NEWSECTION_RE = re.compile(r"\s*\[")


TITLE_RE = re.compile(r"^\s*\[\s*title\s*\]\s*$", re.I)
ORCA_TITLE_RE = re.compile(r"^\s*Molden\s+file\s+created\s+by\s+orca_2mkl\b", re.I)


def convention_from_title(lines):
    for line in lines:
        if NEWSECTION_RE.match(line):
            break
        if ORCA_TITLE_RE.match(line):
            return 'orca'
    return 'standard'


def preview_molden_convention(path):
    """Inspect only the preamble in the file browser. Later titles resolve on import."""
    with open(path) as stream:
        for line in islice(stream, 80):
            if TITLE_RE.match(line):
                return convention_from_title(stream)
            if ATOMS_RE.match(line) or GTO_RE.match(line) or MO_RE.match(line):
                break
    return None


def molden_read_atoms(f, ignore_h=False, animate=False):
    """reads [Atoms] section of molden files"""
    out = []
    line = f.restore_mark("atoms")
    m = ATOMS_RE.search(line)

    # conversion factor incase atomic units are input
    unit = "au" if m.group(1) and m.group(1).lower() in "au" else "angs"
    factor = 1.0 if unit == "angs" else bohr2ang

    seen = set()
    for row, fields in enumerate(_molden_section(f, "atoms"), 1):
        section = f"[Atoms] atom row {row}"
        context = dict(row=row)
        try:
            element, identifier, number, x, y, z = fields
            atomic_number = molden_float(number)
            index = int(identifier)
            section = f"[Atoms] atom {index}"
            context["atom"] = index
            if index in seen:
                raise ValueError("Duplicate atom ID")
            element = element.lower()
            if element not in elements:
                raise ValueError(f"Unknown element {element!r}")
            pos = [factor * molden_float(v) for v in (x, y, z)]
        except ValueError as exc:
            raise f.error(section, str(exc), **context) from exc
        seen.add(index)
        gradient = [0., 0., 0.]
        atom = make_atom_dict(element, pos, index, 0., gradient, ignore_h and element == "h")
        atom["nuclear_charge"] = atomic_number
        if animate:
            atom["trajectory"].append(make_snap_dict(pos, 0., gradient))
        out.append(atom)

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
    reference = None
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
            f.context = dict(atom=reference)
            if reference < 1 or reference in blocks:
                raise f.error(f"[GTO] atom {reference}", "Invalid or duplicate GTO atom reference")
            if atom_ids is not None and reference not in atom_ids:
                raise f.error(f"[GTO] atom {reference}", "Unknown GTO atom reference")
            current = blocks[reference] = []
            continue
        section = "[GTO]" if reference is None else f"[GTO] atom {reference}"
        shell = fields[0].lower()
        if current is None or shell not in ("s", "sp", "p", "d", "f", "g"):
            raise f.error(section, "Invalid GTO shell header: {}".format(line.strip()))
        try:
            nprim = int(fields[1])
            if nprim < 1:
                raise ValueError("Primitive count must be positive")
        except (IndexError, ValueError) as exc:
            raise f.error(section, "Invalid GTO primitive count: {}".format(line.strip())) from exc
        exponents, coefficients, pcoefficients = [], [], []
        for primitive in range(1, nprim + 1):
            try:
                row = f.readline().split()
                if len(row) != (3 if shell == "sp" else 2):
                    raise ValueError("Incorrect primitive column count")
                numbers = [molden_float(x) for x in row]
                if numbers[0] <= 0:
                    raise ValueError("Primitive exponents must be positive")
            except (EOFError, ValueError) as exc:
                raise f.error(section, f"Invalid or truncated GTO {shell} primitive {primitive} of {nprim}") from exc
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
    f.context = {}
    return [blocks.get(reference, []) for reference in references]


def molden_read_mo(f, nao):
    """Read complete indexed MO coefficient blocks without crossing sections."""
    f.restore_mark("mo")
    if nao < 1:
        raise f.error("[MO]", "Molden MO section requires at least one basis function")
    metadata = re.compile(r"\s*(Sym|Ene|Spin|Occup)\s*=\s*(\S+)", re.IGNORECASE)
    out, modata, coefficients = [], {}, []

    def error(message):
        return f.error(f"[MO] orbital {len(out) + 1}", message, orbital=len(out)+1)

    def finish():
        if modata or coefficients:
            if len(coefficients) != nao:
                raise error("Molden MO expected {} coefficients, found {}".format(
                    nao, len(coefficients)))
            if 'occup' not in modata:
                raise error("Missing orbital occupation (Occup)")
            if out and ('spin' in modata) != ('spin' in out[0]):
                raise error("Incomplete spin labels across orbitals")
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
                try:
                    value = str(molden_float(value))
                except ValueError as exc:
                    raise error(f"Invalid {key} value: {value}") from exc
            if key == 'occup' and not 0 <= float(value) <= 2:
                raise error(f"Orbital occupation must be between 0 and 2, found {value}")
            if key == 'spin' and value.lower() not in ('alpha', 'beta'):
                raise error(f"Unknown spin label {value!r}")
            modata[key] = value
        else:
            fields = line.split()
            try:
                if len(fields) != 2 or int(fields[0]) != len(coefficients) + 1:
                    raise ValueError(f"Expected index {len(coefficients) + 1} and one coefficient")
                if len(coefficients) >= nao:
                    raise ValueError("Too many coefficients")
                coefficients.append(molden_float(fields[1]))
            except ValueError as exc:
                raise error(f"Invalid Molden MO coefficient line {line.strip()!r}: {exc}") from exc
    finish()
    if not out:
        raise f.error("[MO]", "Empty orbital section")
    return out


def _molden_section(f, name):
    f.restore_mark(name)
    while True:
        try:
            line = f.readline()
        except EOFError:
            return
        if NEWSECTION_RE.match(line):
            return
        if line.strip():
            yield line.split()


def _bohr_vector(f, section, fields):
    try:
        if len(fields) != 3:
            raise ValueError()
        return [bohr2ang * molden_float(v) for v in fields]
    except ValueError as exc:
        raise f.error(section, 'Expected three finite Cartesian components') from exc


def molden_read_vibrations(f, out, ignore_h, animate):
    atoms = out['atoms']
    equilibrium = []
    if f.is_marked('fr-coord'):
        has_atoms = bool(atoms)
        for i, fields in enumerate(_molden_section(f, 'fr-coord')):
            symbol = fields[0].lower()
            position = _bohr_vector(f, '[FR-COORD]', fields[1:])
            if symbol not in symbols or (has_atoms and (i >= len(atoms) or atoms[i]['symbol'] != symbol)):
                raise f.error('[FR-COORD]', f'Unexpected element {symbol} at atom {i + 1}')
            equilibrium.append(position)
            if not has_atoms:
                atom = make_atom_dict(symbol, position, i+1, 0., [0.] * 3, ignore_h and symbol == 'h')
                if animate:
                    atom['trajectory'].append(make_snap_dict(position, 0., [0.] * 3))
                atoms.append(atom)
        if not equilibrium or len(equilibrium) != len(atoms):
            raise f.error('[FR-COORD]', 'Equilibrium atom count must match the molecule')
    else:
        equilibrium = [atom['position'] for atom in atoms]
    if not f.is_marked('fr-norm-coord'):
        return
    if not f.is_marked('freq'):
        f.restore_mark('fr-norm-coord')
        raise f.error('[FR-NORM-COORD]', 'Missing [FREQ] for vibrational modes')
    frequencies = []
    for fields in _molden_section(f, 'freq'):
        try:
            if len(fields) != 1:
                raise ValueError()
            frequencies.append(molden_float(fields[0]))
        except ValueError as exc:
            raise f.error('[FREQ]', 'Expected one finite frequency per line') from exc

    modes = {}
    number = None
    for fields in _molden_section(f, 'fr-norm-coord'):
        if fields[0].lower() == 'vibration':
            if number is not None and len(modes[number]) != len(atoms):
                raise f.error('[FR-NORM-COORD]', f'Mode {number} needs one displacement per atom')
            if len(fields) != 2 or not fields[1].isdigit():
                raise f.error('[FR-NORM-COORD]', 'Expected vibration followed by the mode number')
            number = int(fields[1])
            if not 1 <= number <= len(frequencies) or number in modes:
                raise f.error('[FR-NORM-COORD]', f'Invalid or repeated mode {number}')
            modes[number] = []
        elif number is None or len(modes[number]) >= len(atoms):
            raise f.error('[FR-NORM-COORD]', 'Expected a vibration header before displacements')
        else:
            modes[number].append(_bohr_vector(f, f'[FR-NORM-COORD] mode {number}', fields))
    if number is None or len(modes[number]) != len(atoms):
        raise f.error('[FR-NORM-COORD]', f'Incomplete displacement block for mode {number}')
    if len(modes) != len(frequencies):
        raise f.error('[FR-NORM-COORD]', f'Found {len(modes)} modes for {len(frequencies)} frequencies')
    out['vibrations'] = dict(equilibrium=equilibrium, frequencies=frequencies,
                             displacements=[modes[i] for i in range(1, len(frequencies)+1)])


@stopwatch("read molden")
def molecule_from_molden(filename, _options):
    """Read molden file to look for coordinates and optionally orbital data"""
    # molden format is not terribly efficient, but hopefully this doesn't matter
    ignore_h = _options.get("ignore_hydrogen", False)
    animate = _options.get("plot_type", "") == "animate"
    out = {"atoms": [], "source": os.fspath(filename)}
    mode = molden_convention(_options)

    marks = {"molden": MOLDEN_RE, "title": TITLE_RE,
             "atoms": ATOMS_RE,
             "core": re.compile(r"\[\s*(?:core|ecp)\s*\]", re.I),
             "freq": re.compile(r"\[\s*FREQ\s*\]", re.I),
             "fr-coord": re.compile(r"\[\s*FR-COORD\s*\]", re.I),
             "fr-norm-coord": re.compile(r"\[\s*FR-NORM-COORD\s*\]", re.I),
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
                if NEWSECTION_RE.match(line):
                    for x in marks:
                        if marks[x].search(line):
                            f.set_mark(x, line)
        except EOFError:
            pass

        if not f.is_marked("molden"):
            raise InputError(filename, "Molden header", "Molden file is missing [Molden Format] specification")
        if mode == 'auto':
            mode = 'standard'
            if f.is_marked('title'):
                f.restore_mark('title')
                mode = convention_from_title(f.f)
        out['molden_convention'] = mode
        out['denorm'] = mode != 'orca'
        out['has_ecp'] = f.is_marked('core')
        if f.is_marked("atoms"):
            out["atoms"] = molden_read_atoms(f, ignore_h, animate=animate)
        elif not f.is_marked("fr-coord"):
            raise InputError(filename, "Molden geometry", "Molden file is missing [Atoms] or [FR-COORD]")

        if f.is_marked("fr-coord") or f.is_marked("fr-norm-coord"):
            molden_read_vibrations(f, out, ignore_h, animate)

        if not out["atoms"]:
            raise f.error("[Atoms]", "Missing atom records")

        if f.is_marked("gto"):  # read through GTO section to build basis
            out["basis"] = molden_read_gto(f, [atom["index"] for atom in out["atoms"]])

        # read through MO section to build coefs
        if f.is_marked("mo") and "basis" not in out:
            f.restore_mark("mo")
            raise f.error("[MO]", "Missing [GTO] basis for orbital coefficients")
        if f.is_marked("mo"):
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

            if mode == 'turbomole_legacy':
                # tm2molden norm scales the Cartesian MO components of each
                # expanded shell together. Explicit spherical shells use Molden normalization.
                cart_factors = {'d': np.sqrt(3.), 'f': np.sqrt(15.), 'g': np.sqrt(105.)}
                factors = []
                for atom_basis in out['basis']:
                    for shell in atom_basis:
                        kind = shell['shell']
                        factor = 1. if kind in spherical else cart_factors.get(kind, 1.)
                        factors.extend([factor] * shelldegen[kind])
                factors = np.asarray(factors)
                for orbital in molden_mo:
                    orbital['coeff'] = (factors * orbital['coeff']).tolist()

            if spherical and molden_mo:
                orca = mode == 'orca'
                if orca:
                    # orca_2mkl includes the xy, xyz and xxyz primitive
                    # factors in d, f and g contractions. Our shared radial
                    # factor is the one for xx, xxx or xxxx respectively.
                    factors = {"d": np.sqrt(3.), "f": np.sqrt(15.), "g": np.sqrt(35.)}
                    for atom_basis in out["basis"]:
                        for shell in atom_basis:
                            if shell["shell"] in spherical:
                                factor = factors[shell["shell"]]
                                shell["contractions"] = [c / factor for c in shell["contractions"]]
                molden_mo = transform_sph_to_cart(molden_mo, spherical, out["basis"], orca=orca)

            out["mo"] = molden_mo

    return out


def _cube_origin_header(f):
    f.section = "Cube atom count and origin"
    header = f.readline().split()
    if len(header) not in (4, 5):
        raise ValueError('Cube atom/origin header requires four or five values')
    nvalues = int(header[4]) if len(header) == 5 else 1
    if nvalues < 1:
        raise ValueError('Cube dataset count must be positive')
    if nvalues > sys.maxsize:
        raise ValueError('Cube dataset count exceeds the supported range')
    return int(header[0]), [_finite_float(t, a) for t, a in zip(header[1:4], 'xyz')], nvalues


def _cube_labels(f, natoms, nvalues):
    f.section = "Cube dataset IDs"
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
        f.section = "Cube comments"
        f.readline()
        f.readline()
        natoms, _, nvalues = _cube_origin_header(f)
        if natoms < 0:
            f.section = "Cube axes and atoms"
            for _ in range(3 + abs(natoms)):
                f.readline()
        return _cube_labels(f, natoms, nvalues)


def _cube_values(f, size):
    """Convert in batches, locating individual bad tokens only on the error path."""
    data = np.empty(size, dtype=np.float32)
    used = 0
    f.section = "Cube field values"
    while lines := list(islice(f.f, 256)):
        tokens = ' '.join(lines).split()
        end = used + len(tokens)
        try:
            if end > size:
                raise ValueError
            with np.errstate(over='ignore'):
                data[used:end] = tokens
            if not np.isfinite(data[used:end]).all():
                raise ValueError
        except ValueError:
            for line in lines:
                f.line_number += 1
                for token in line.split():
                    used += 1
                    if used > size:
                        raise ValueError(f"Cube expected {size} field values, found extra value {token!r}")
                    value = _finite_float(token, f"field #{used}")
                    if abs(value) > float(np.finfo(np.float32).max):
                        raise ValueError(f"Field value #{used} exceeds float32 range: {token!r}")
            raise
        used = end
        f.line_number += len(lines)
    if used != size:
        f.line_number += 1
        raise ValueError(f"Cube expected {size} field values, found {used} ({size-used} missing)")
    return data


@stopwatch("read cube")
def molecule_from_cube(filename, options):
    """Read a cube file including its volumetric data"""
    out = {"atoms": [], "volume": {"source": os.fspath(filename)}}

    ignore_h = options.get("ignore_hydrogen", False)
    animate = options.get("plot_type", "") in [ "animate", "auto" ]

    with Reader(filename) as f:
        f.section = "Cube comments"
        out["volume"]["comments"] = [f.readline().rstrip("\r\n"), f.readline().rstrip("\r\n")]

        # 3rd line is <natoms> <origin_x> <origin_y> <origin_z>
        natoms, origin, nvalues = _cube_origin_header(f)

        nres = [0, 0, 0]  # number of points in each direction
        # axes[i,:] defines the i-th axis
        axes = np.zeros([3, 3], dtype=np.float32)

        # 4th, 5th, 6th lines describe the axes: <npoints> <axis_x> <axis_y> <axis_z>
        for i in range(3):
            f.section = f"Cube grid axis {i+1}"
            res, vx, vy, vz = f.readline().split()
            nres[i] = int(res)
            axes[i, :] = np.array([vx, vy, vz], dtype=np.float32)
            if not (np.isfinite(axes[i]).all() and axes[i].any()):
                raise ValueError("Cube grid axes must be finite and nonzero")
            if nres[i] == 0 or (nres[i] > 0) != (nres[0] > 0):
                raise ValueError("Cube voxel counts must be nonzero and have the same sign")
        if all(n > 0 for n in nres):
            length_unit = bohr2ang
        else:
            length_unit = 1.0
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
            f.context = dict(atom=i)
            f.section = f"Cube atom {i}"
            ato, chg, x, y, z = f.readline().split()
            iatom = int(ato)
            if not 0 <= iatom < len(symbols):
                raise ValueError(f"Unknown atomic number {ato!r}")
            position = [length_unit * _finite_float(t, a) for t, a in zip((x, y, z), 'xyz')]
            chg = _finite_float(chg, "charge")
            hidden = ignore_h and symbols[iatom] == "h"
            new_atom = make_atom_dict(symbols[iatom], position, i, chg,
                    [0.0, 0.0, 0.0], hidden)

            if animate:
                new_atom["trajectory"] = [make_snap_dict(position, 0.0, [0.0, 0.0, 0.0])]

            out["atoms"].append(new_atom)

        # optionally check for DSET_IDs, basically labels of different data sets
        f.context = {}
        volume_labels = list(_cube_labels(f, natoms, nvalues))

        # and finally the volumetric data comes, 6 elements per line in z, y, x order
        nvol = len(volume_labels)
        data_layout = nres + [ nvol ]
        data = _cube_values(f, math.prod(data_layout))

    datasets = data.reshape(data_layout)
    selected = options.get('cube_dataset', 0)
    if isinstance(selected, bool) or not isinstance(selected, int) or not 0 <= selected < nvol:
        raise InputError(filename, 'Cube datasets', f'cube_dataset must be a zero-based index in [0, {nvol})')
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
    resample = entry.get("resample", False)
    if not isinstance(resample, bool):
        raise ValueError("resample must be true or false")
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
            raise InputError(os.path.join(directory, name), f"Term {index} ({name})", str(exc), term=index) from exc
        grid = cube["volume"]
        if result is None:
            result = cube
            template = grid
            total = np.zeros(grid["nres"], dtype=np.float64)
        different_grid = (grid["nres"] != template["nres"]
              or not np.allclose(grid["origin"], template["origin"], atol=1e-5, rtol=0)
              or not np.allclose((grid["axes"] - template["axes"]) * (np.array(grid["nres"]) - 1)[:, None],
                                 0, atol=1e-5, rtol=0))
        if different_grid:
            if not resample:
                raise ValueError(f"Term {index} ({name}) has a different grid shape, origin or axes")
            try:
                add_resampled_cube(total, template, grid, weight)
            except ValueError as exc:
                raise ValueError(f"Term {index} ({name}) must cover the first cube's grid: {exc}") from exc
        else:
            total += np.multiply(grid["data"], weight, dtype=np.float64)
        # Only geometry and grid metadata are needed from the first term.
        del grid["data"], grid["datasets"]
    if not np.all(np.isfinite(total)) or np.max(np.abs(total)) > np.finfo(np.float32).max:
        raise ValueError("Linear combination exceeds the finite cube value range")
    values = total.astype(np.float32)
    template.update(data=values, datasets=values[..., None], labels=[0], dataset_index=0,
                    field_type="signed", units=units, linear_combination=terms)
    if resample:
        template["resample"] = True
    return result


def molecule_from_json(filename, options, *, single_combination=False):
    """Read one structure per manifest entry, with zero-based timeline indices."""
    with open(filename) as source:
        try:
            data = json.load(source)
        except json.JSONDecodeError as exc:
            raise InputError(filename, "JSON", exc.msg, line=exc.lineno, column=exc.colno) from exc
    if not isinstance(data, dict):
        raise InputError(filename, "JSON", "Expected a JSON object")
    entries = data.get("molecules")
    if single_combination and not (isinstance(entries, list) and len(entries) == 1
            and isinstance(entries[0], dict) and "linear_combination" in entries[0]):
        raise InputError(filename, "JSON", "field mapping needs one linear_combination entry in molecules")
    if "molecules" not in data:
        return None  # Historical direct reader API.
    if not isinstance(entries, list) or not entries:
        raise InputError(filename, "JSON", "Expected a nonempty molecules list")

    start, end, step = trajectory_range(options)
    previous = -1
    frames = []
    merged = None
    for i, entry in enumerate(entries):
        label = f"entry {i + 1}"
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
                current["volume"]["source"] = os.fspath(filename)
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
                    if not isinstance(entry["polar"], str) or not entry["polar"].strip():
                        raise ValueError("Expected a DIPL filename in polar")
                    orbital_pairs = None
                    if "mo" in current:
                        occupied = sum(_finite_float(orb.get("occup", ""), "Molden occupation") > 0
                                       for orb in current["mo"])
                        orbital_pairs = occupied * (len(current["mo"]) - occupied)
                    current["polar"] = read_turbomole_polarizability(
                        os.path.join(directory, entry["polar"]), orbital_pairs=orbital_pairs)
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
        except (ValueError, OSError, OverflowError) as exc:
            raise InputError(filename, label, str(exc), entry=i+1) from exc

    if merged is None:
        raise InputError(filename, "JSON", "Frame range selects no JSON entries")
    merged["frame_numbers"] = frames
    options["plot_type"] = "animate"
    return merged
