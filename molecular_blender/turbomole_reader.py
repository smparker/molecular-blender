# -*- coding: utf-8 -*-
#
#  Molecular Blender
#  Filename: turbomole_reader.py
#  Copyright (C) 2025 Shane Parker
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

import math
import re

from typing import Dict, List, Union
from dataclasses import dataclass

import numpy as np

from .diagnostics import InputError

PAIR_HEADER = re.compile(r'(\d+)\s+eigenvalue\s*=\s*(\S+)')


def response_headers(path):
    """Preview the eigenpair IDs and values without allocating response vectors."""
    pairs = []
    inside = False
    with open(path) as source:
        for row in source:
            line = row.strip()
            if line.startswith('$'):
                inside = line == '$eigenpairs'
            elif inside and 'eigenvalue' in line:
                match = PAIR_HEADER.fullmatch(line)
                if not match:
                    raise ValueError('Invalid eigenpair header')
                pairs.append((int(match[1]), DiplParser._parse_scientific_notation(match[2])))
    return tuple(pairs)

@dataclass
class Eigenpair:
    index: int
    eigenvalue: float
    vector: List[float]

class DiplParser:
    def __init__(self, file_content: str, filename="<DIPL>"):
        self.filename = filename
        self.sections: Dict[str, Union[str, int, List[Eigenpair]]] = {}
        self.section_lines = {}
        self.pair_lines = {}
        lines = file_content.splitlines()
        self.eof = len(lines) + 1
        self._parse(lines)

    def error(self, line, section, message):
        return InputError(self.filename, "$" + section, message, line=line)

    @staticmethod
    def _parse_scientific_notation(value: str) -> float:
        number = float(value.replace('D', 'E').replace('d', 'e'))
        if not math.isfinite(number):
            raise ValueError(f"Non-finite number {value!r}")
        return number

    def _parse_fixed_width_numbers(self, line: str, width: int = 20) -> List[float]:
        return [self._parse_scientific_notation(line[i:i + width])
                for i in range(0, len(line), width)]

    def _parse(self, lines):
        section = None
        pair = None
        integers = ('tensor space dimension', 'current subspace dimension')
        inline = integers + ('symmetry', 'scfinstab', 'current iteration')
        for number, raw in enumerate(lines, 1):
            line = raw.strip()
            if not line:
                continue
            if line.startswith('$'):
                header = line[1:].strip()
                section = next((key for key in inline if header == key or header.startswith(key + ' ')
                                or header.startswith(key + '\t')), header)
                self.section_lines[section] = number
                pair = None
                if section == 'eigenpairs':
                    self.sections[section] = []
                line = header[len(section):].strip()
                if not line:
                    continue
            if section == 'eigenpairs':
                match = PAIR_HEADER.fullmatch(line)
                if match:
                    index = int(match[1])
                    if index in self.pair_lines:
                        raise self.error(number, section, f"Duplicate eigenpair {index}")
                    try:
                        value = self._parse_scientific_notation(match[2])
                    except ValueError as exc:
                        raise self.error(number, f"eigenpairs, eigenpair {index}", str(exc)) from exc
                    pair = Eigenpair(index, value, [])
                    self.pair_lines[index] = number
                    self.sections[section].append(pair)
                else:
                    if 'eigenvalue' in line or pair is None:
                        raise self.error(number, section, f"Expected an eigenpair header, found {line!r}")
                    try:
                        pair.vector.extend(self._parse_fixed_width_numbers(line))
                    except ValueError as exc:
                        raise self.error(number, f"eigenpairs, eigenpair {pair.index}", str(exc)) from exc
            elif section in integers:
                try:
                    self.sections[section] = int(line)
                except ValueError as exc:
                    raise self.error(number, section, f"Invalid dimension {line!r}") from exc
            elif section and section != 'end':
                self.sections[section] = line.split()[-1]

    def get_section(self, section_name: str):
        return self.sections.get(section_name)


def read_dipl_file(file_path: str):
    with open(file_path) as source:
        parser = DiplParser(source.read(), file_path)
    sym = parser.get_section('symmetry')
    if sym != "c1":
        raise parser.error(parser.section_lines.get('symmetry', parser.eof), 'symmetry',
                           f"Only c1 symmetry is supported, received {sym}")
    return parser


def read_turbomole_polarizability(dipl_file, *, orbital_pairs=None):
    """Read one response vector per eigenpair, retaining the first RPA component."""
    dipl = read_dipl_file(dipl_file)
    pairs = dipl.get_section('eigenpairs')
    if not pairs:
        raise dipl.error(dipl.section_lines.get('eigenpairs', dipl.eof), 'eigenpairs',
                         'Expected at least one eigenpair')
    nrpa = 1 if dipl.get_section('scfinstab') == "polly" else 2
    size = dipl.get_section('tensor space dimension')
    polarizations = []
    for pair in pairs:
        count = len(pair.vector)
        if not count or count % nrpa or (size is not None and count != nrpa * size):
            expected = f"{nrpa * size} values" if size is not None else f"{nrpa} equally sized nonempty components"
            raise dipl.error(dipl.pair_lines[pair.index], f"eigenpairs, eigenpair {pair.index}",
                             f"Expected {expected}, found {count} vector values")
        if orbital_pairs is not None and count // nrpa != orbital_pairs:
            raise dipl.error(dipl.pair_lines[pair.index], f"eigenpairs, eigenpair {pair.index}",
                             f"Molden file has {orbital_pairs} occupied–virtual pairs, "
                             f"response has {count // nrpa}")
        vec = np.asarray(pair.vector, dtype=np.float32).reshape(nrpa, -1)
        polarizations.append({"index": pair.index, "eigenvalue": pair.eigenvalue, "vector": vec[0]})
    return polarizations
