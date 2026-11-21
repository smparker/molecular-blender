"""SDF/MOL V2000 parsing: explicit bond orders, validation and record selection."""
from pathlib import Path
import tempfile
import unittest

from molecular_blender.importers import molecule_from_file

FIXTURE = Path(__file__).resolve().parents[1] / 'examples' / 'bond_orders.sdf'


class TestSDF(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.directory.cleanup)
        self.lines = FIXTURE.read_text().splitlines()

    def read(self, lines=None, suffix='.sdf', **options):
        path = Path(self.directory.name) / ('model' + suffix)
        path.write_text('\n'.join(self.lines if lines is None else lines) + '\n')
        return molecule_from_file(str(path), options)

    def test_orders_and_coordinates(self):
        options = {}
        model = molecule_from_file(str(FIXTURE), options)
        self.assertEqual(len(model['atoms']), 9)
        self.assertEqual([b['order'] for b in model['bonds']], [2, 1, 2, 1, 3, 1, 1, 1])
        self.assertEqual(model['bonds'][0]['atoms'], [0, 1])
        self.assertEqual(model['atoms'][0]['position'], [-3.4, 1.05, 0])
        self.assertEqual(options['plot_type'], 'frame')
        self.assertTrue(all(not a['trajectory'] for a in model['atoms']))

    def test_mol_and_sd_extensions(self):
        end = self.lines.index('M  END')
        for extension in ('.mol', '.sd', '.SDF'):
            with self.subTest(extension=extension):
                self.assertEqual(len(self.read(self.lines[:end + 1], extension)['bonds']), 8)

    def test_hidden_hydrogens_keep_indices_and_bonds(self):
        model = self.read(ignore_hydrogen=True)
        self.assertEqual([a['index'] for a in model['atoms']], list(range(9)))
        self.assertEqual([a['hidden'] for a in model['atoms']], [False]*6 + [True]*3)
        self.assertEqual(len(model['bonds']), 8)

    def test_records_and_blank_title(self):
        second = self.lines.copy()
        second[0] = ''
        self.assertEqual(self.read(self.lines + second, sdf_record=1)['title'], '')
        self.assertEqual(self.read(self.lines + second)['title'], self.lines[0])
        for record in (2, -1, 0.5, True):
            with self.subTest(record=record), self.assertRaises(ValueError):
                self.read(self.lines + second, sdf_record=record)

    def test_formal_charges(self):
        self.lines[4] = self.lines[4][:36] + '  3' + self.lines[4][39:]
        atom = self.read()['atoms'][0]
        self.assertEqual(atom['formal_charge'], 1)
        self.assertEqual(atom['charge'], 0)
        end = self.lines.index('M  END')
        self.lines[end:end] = ['M  CHG  1   2  -1', 'M  CHG  1   6   1']
        self.assertEqual([a['formal_charge'] for a in self.read()['atoms']], [0, -1, 0, 0, 0, 1, 0, 0, 0])

    def test_empty_bond_table(self):
        lines = self.lines[:13] + ['M  END']
        lines[3] = lines[3][:3] + '  0' + lines[3][6:]
        self.assertEqual(self.read(lines)['bonds'], [])

    def test_invalid_inputs(self):
        cases = []
        for replacement in ('  0  2  1  0', ' 10  2  1  0', '  1  1  1  0',
                            '  1  2  4  0', '  1  2  8  0', 'bad'):
            lines = self.lines.copy()
            lines[13] = replacement
            cases.append(lines)
        lines = self.lines.copy()
        lines[14] = lines[13]
        cases.append(lines)  # duplicate bond
        lines = self.lines.copy()
        lines[3] = lines[3].replace('V2000', 'V3000')
        cases.append(lines)
        for atom_line in ('bad', f'{float("nan"):10}{0:10.4f}{0:10.4f} C   0  0'):
            lines = self.lines.copy()
            lines[4] = atom_line
            cases.append(lines)
        cases.extend([self.lines[:8], [l for l in self.lines if l != 'M  END']])
        lines = self.lines.copy()
        lines.insert(21, 'M  CHG  1  99   1')
        cases.append(lines)
        lines = self.lines.copy()
        lines.insert(21, 'M  CHG')
        cases.append(lines)
        for index, lines in enumerate(cases):
            with self.subTest(case=index), self.assertRaises(ValueError):
                self.read(lines)

    def test_diagnostics_use_physical_line_in_second_record(self):
        for offset, bad, context in [(4, 'bad', 'atom 1'), (13, '  1 88  1', 'bond 1'),
                                     (21, 'M  CHG  1  99   1', 'M  CHG')]:
            second = self.lines.copy()
            second[offset] = bad
            physical = len(self.lines) + offset + 1
            with self.subTest(context=context), self.assertRaisesRegex(
                    ValueError, f'model.sdf:{physical}: SDF record 1, {context}'):
                self.read(self.lines + second, sdf_record=1)
