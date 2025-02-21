import numpy as np

import re
from typing import Dict, List, Union
from dataclasses import dataclass

from dataclasses import dataclass
from typing import Dict, List, Union

@dataclass
class Eigenpair:
    index: int
    eigenvalue: float
    vector: List[float]

class DiplParser:
    def __init__(self, file_content: str):
        """
        Initialize the parser with file content and parse the known sections:
        - title: A simple text description
        - symmetry: The symmetry group
        - tensor space dimension: An integer value
        - scfinstab: A parameter value
        - current subspace dimension: An integer value
        - current iteration: Status of the iteration
        - eigenpairs: The eigenvalues and eigenvectors
        """
        self.content = file_content
        self.sections: Dict[str, Union[str, int, List[Eigenpair]]] = {}
        self._parse()

    def _parse_scientific_notation(self, value: str) -> float:
        """
        Convert scientific notation with 'D' exponent marker to float.
        Handles both positive and negative values.
        """
        return float(value.replace('D', 'E'))

    def _parse_fixed_width_numbers(self, line: str, width: int = 20) -> List[float]:
        """
        Parse a line of fixed-width scientific notation numbers.

        Each number takes up exactly 22 characters, with the sign of the next number
        attached at the end. For example, the line:
        "0.14237842212942D-03-.22577892444500D-03"
        contains two numbers:
        "0.14237842212942D-03" and "-0.22577892444500D-03"

        Args:
            line: String containing concatenated fixed-width numbers
            width: Width of each number field (default 20)

        Returns:
            List of parsed float values
        """
        numbers = []
        i = 0
        while i < len(line):
            # Extract the number string, including any leading negative sign
            if i == 0:
                # First number might have a sign
                num_str = line[i:i + width]
            else:
                # Subsequent numbers' signs are at the end of the previous field
                sign = line[i]
                num_str = line[i:i + width]
                #if sign == '-':
                #    num_str = '-' + num_str

            numbers.append(self._parse_scientific_notation(num_str.strip()))
            i += width  # -1 because signs overlap between fields

        return numbers

    def _parse_eigenpairs(self, content: str) -> List[Eigenpair]:
        """
        Parse the eigenpairs section containing eigenvalues and their vectors.

        The format has two parts:
        1. A header line with index and eigenvalue:
           "        1   eigenvalue =  0.6617507306993266D-01"
        2. Multiple lines of fixed-width scientific notation numbers:
           "0.14237842212942D-03-.22577892444500D-03..."

        Each number in the vector takes exactly 22 characters, with signs serving
        as field separators.
        """
        lines = content.strip().split('\n')

        index = None
        eigenvalue = None

        # Parse the vector values from the remaining lines
        vector_values = None

        eigenpairs = []
        for line in lines:
            # check whether the next line starts with "$" or "   2 eigenvalue" or whatever
            if line.startswith("$") or re.match(r"\s*\d+\s+eigenvalue\s*=\s*(\S+)", line):
                if vector_values:
                    ep = Eigenpair(index, eigenvalue, vector_values)
                    eigenpairs.append(ep)

                index = int(line.split()[0])
                eigenvalue = self._parse_scientific_notation(line.split('=')[1].strip())
                vector_values = []

                continue

            if line.strip():  # Skip empty lines
                values = self._parse_fixed_width_numbers(line.strip())
                vector_values.extend(values)

        if vector_values:
            ep = Eigenpair(index, eigenvalue, vector_values)
            eigenpairs.append(ep)

        return eigenpairs

    def _parse(self):
        """
        Parse each known section in order. We know exactly which sections
        to expect and their format:
        1. Simple text sections (title, symmetry, scfinstab)
        2. Dimension sections that contain integers
        3. Status section (current iteration)
        4. Complex eigenpairs section
        """
        lines = self.content.strip().split('\n')
        current_section = None
        section_content = []

        for line in lines:
            if line.startswith('$symmetry'):
                self.sections['symmetry'] = line.split()[-1]
            elif line.startswith('$tensor space dimension'):
                self.sections['tensor space dimension'] = int(line.split()[-1])
            elif line.startswith('$current subspace dimension'):
                self.sections['current subspace dimension'] = int(line.split()[-1])
            elif line.startswith('$current iteration'):
                self.sections['current iteration'] = line.split()[-1]
            elif line.startswith('$scfinstab'):
                self.sections['scfinstab'] = line.split()[-1]
            elif line.startswith('$'):
                # Save the previous section if we were collecting one
                if current_section == 'eigenpairs' and section_content:
                    self.sections[current_section] = self._parse_eigenpairs('\n'.join(section_content))
                elif current_section and section_content:
                    content = ' '.join(section_content).strip()
                    # Handle different section types
                    if 'dimension' in current_section:
                        self.sections[current_section] = int(content.split()[-1])
                    else:
                        self.sections[current_section] = content.split()[-1]

                # Start new section
                current_section = line[1:].strip()  # Remove $ and whitespace
                section_content = []
            elif current_section and current_section != 'end':
                section_content.append(line)

        # Don't forget to process the last section if it wasn't 'end'
        if current_section == 'eigenpairs' and section_content:
            self.sections[current_section] = self._parse_eigenpairs('\n'.join(section_content))

    def get_section(self, section_name: str) -> Union[str, int, List[Eigenpair]]:
        """Retrieve the parsed content of a specific section."""
        return self.sections.get(section_name)

def read_dipl_file(file_path: str):
    """
    Read and parse a Turbomole DIPL file.

    Args:
        file_path: Path to the DIPL file
    """
    with open(file_path, 'r') as f:
        content = f.read()

    parser = DiplParser(content)

    sym = parser.get_section('symmetry')
    if sym != "c1":
        raise ValueError(f"Only c1 symmetry is supported, received {sym}")

    tsize = parser.get_section('tensor space dimension')

    nvec = parser.get_section('current subspace dimension')

    eigenpairs = parser.get_section('eigenpairs')

    return parser

def read_turbomole_polarizability(dipl_file):
    """
    Read the polarizability tensor from a Turbomole DIPL file.

    Args:
        dipl_file: Path to the DIPL file
    """
    dipl = read_dipl_file(dipl_file)

    nrpa = 1 if dipl.get_section('scfinstab') == "polly" else 2

    nvec = dipl.get_section('current subspace dimension')

    # Extract the polarization tensors
    polarizations = []
    for i, pair in enumerate(dipl.get_section('eigenpairs')):
        eigenvalue = pair.eigenvalue
        vec = np.array(pair.vector).astype(np.float32).reshape(nrpa, -1)
        polarizations.append({ "eigenvalue": eigenvalue, "vector": vec[0] })

    return polarizations
