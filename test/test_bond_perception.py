"""Valence search tests independent of Blender, plus an optional RDKit oracle."""
import unittest

from molecular_blender.bond_perception import perceive_bond_orders, BondPerceptionError


class TestBondPerception(unittest.TestCase):
    def test_simple_molecules(self):
        cases = [
            (['h', 'h'], [(0, 1)], [1]),
            (['n', 'n'], [(0, 1)], [3]),
            (['c', 'o', 'o'], [(0, 1), (0, 2)], [2, 2]),
            (['c', 'o', 'h', 'h'], [(0, 1), (0, 2), (0, 3)], [2, 1, 1]),
            (['c', 'c', 'h', 'h'], [(0, 1), (0, 2), (1, 3)], [3, 1, 1]),
            (['c', 'c'] + ['h']*4, [(0, 1), (0, 2), (0, 3), (1, 4), (1, 5)], [2, 1, 1, 1, 1]),
            (['o', 'h', 'h'], [(0, 1), (0, 2)], [1, 1]),
        ]
        for symbols, edges, expected in cases:
            with self.subTest(symbols=symbols):
                orders, charges = perceive_bond_orders(symbols, edges)
                self.assertEqual(orders, expected)
                self.assertEqual(charges, [0]*len(symbols))

    def test_backtracks_from_shortest_edge(self):
        # A-B-C-D: each carbon lacks one bond. Trying the short central edge
        # first strands the terminal carbons; neutral solution doubles AB/CD.
        edges = [(0, 1), (1, 2), (2, 3), (0, 4), (0, 5),
                 (1, 6), (2, 7), (3, 8), (3, 9)]
        orders, charges = perceive_bond_orders(['c']*4 + ['h']*6, edges,
                                               bond_lengths=[1, .8, 1] + [1]*6)
        self.assertEqual(orders, [2, 1, 2] + [1]*6)
        self.assertEqual(charges, [0]*10)

    def test_ions_and_charge_separation(self):
        for symbols, edges, charge, expected_q in [
            (['n'] + ['h']*4, [(0, i) for i in range(1, 5)], 1, [1, 0, 0, 0, 0]),
            (['o', 'h'], [(0, 1)], -1, [-1, 0]),
            (['cl'], [], -1, [-1]),
        ]:
            with self.subTest(symbols=symbols):
                _, charges = perceive_bond_orders(symbols, edges, charge=charge)
                self.assertEqual(charges, expected_q)
        # Nitromethane requires an N+/O- pair even though the total is neutral.
        symbols = ['c', 'n', 'o', 'o', 'h', 'h', 'h']
        edges = [(0, 1), (1, 2), (1, 3), (0, 4), (0, 5), (0, 6)]
        orders, charges = perceive_bond_orders(symbols, edges)
        self.assertEqual(sorted(orders[1:3]), [1, 2])
        self.assertEqual(charges[1], 1)
        self.assertEqual(sorted(charges[2:4]), [-1, 0])
        self.assertEqual(sum(charges), 0)

    def test_rings(self):
        for count, charge in ((6, 0), (5, -1)):
            symbols = ['c']*count + ['h']*count
            edges = [(i, (i+1) % count) for i in range(count)] + [(i, i+count) for i in range(count)]
            orders, charges = perceive_bond_orders(symbols, edges, charge=charge)
            self.assertEqual(orders.count(2), count//2)
            self.assertEqual(sum(charges), charge)
            self.assertTrue(all(orders[i] + orders[(i-1) % count] <= 3 for i in range(count)))

    def test_impossible_and_unsupported(self):
        cases = [(['c'] + ['h']*5, [(0, i) for i in range(1, 6)]),
                 (['c'] + ['h']*3, [(0, i) for i in range(1, 4)]),
                 (['h', 'h', 'h'], [(0, 1), (0, 2)]),
                 (['fe', 'o'], [(0, 1)])]
        for symbols, edges in cases:
            with self.subTest(symbols=symbols), self.assertRaises(BondPerceptionError):
                perceive_bond_orders(symbols, edges)
        with self.assertRaises(BondPerceptionError):
            perceive_bond_orders(['n'] + ['h']*4, [(0, i) for i in range(1, 5)], charge=0)

    def test_budget_and_input_validation(self):
        edges = [(i, (i+1) % 6) for i in range(6)] + [(i, i+6) for i in range(6)]
        snapshot = edges.copy()
        with self.assertRaisesRegex(BondPerceptionError, 'limit'):
            perceive_bond_orders(['c']*6 + ['h']*6, edges, max_iterations=1)
        self.assertEqual(edges, snapshot)
        for kwargs in ({'charge': .5}, {'max_iterations': 0}, {'bond_lengths': [float('nan')]}):
            with self.subTest(kwargs=kwargs), self.assertRaises(ValueError):
                perceive_bond_orders(['n', 'n'], [(0, 1)], **kwargs)
        for edges in ([(0, 0)], [(0, 2)], [(0, 1), (1, 0)]):
            with self.assertRaises(ValueError):
                perceive_bond_orders(['n', 'n'], edges)

    def test_failure_diagnostics(self):
        cases = [
            (['fe', 'o'], [(0, 1)], 0, 'unsupported_element'),
            (['c'] + ['h']*5, [(0, i) for i in range(1, 6)], 0, 'bad_connectivity'),
            (['c'], [], 0, 'insufficient_connectivity'),
            (['n'] + ['h']*4, [(0, i) for i in range(1, 5)], 0, 'incompatible_charge'),
        ]
        for symbols, edges, charge, reason in cases:
            with self.subTest(reason=reason), self.assertRaises(BondPerceptionError) as raised:
                perceive_bond_orders(symbols, edges, charge=charge,
                                     atom_ids=list(range(10, 10 + len(symbols))))
            error = raised.exception
            self.assertEqual(error.reason, reason)
            self.assertEqual(error.diagnostics['requested_charge'], charge)
            self.assertIn('atom 10', str(error))
            self.assertIn('search nodes', str(error))
            first = error.diagnostics['atoms'][0]
            self.assertEqual(first['element'], symbols[0])
            self.assertEqual(first['degree'], sum(0 in edge for edge in edges))
            self.assertIn('candidate_states', first)
        ring = [(i, (i+1) % 6) for i in range(6)] + [(i, i+6) for i in range(6)]
        with self.assertRaises(BondPerceptionError) as raised:
            perceive_bond_orders(['c']*6 + ['h']*6, ring, max_iterations=1)
        self.assertEqual(raised.exception.reason, 'search_limit')
        self.assertEqual(raised.exception.diagnostics['iterations'], 1)
        self.assertEqual(raised.exception.diagnostics['max_iterations'], 1)


try:
    from rdkit import Chem
    from rdkit.Chem import rdDetermineBonds
    HAS_RDKIT = True
except ImportError:
    HAS_RDKIT = False


@unittest.skipUnless(HAS_RDKIT, 'optional RDKit reference comparison')
class TestRDKitReference(unittest.TestCase):
    def test_common_organic_graphs(self):
        # Compare sanitized canonical structures so equivalent Kekule/resonance
        # representations do not require the same arbitrary edge numbering.
        cases = ['C', 'CC', 'C=C', 'C#C', 'C=O', 'O=C=O', 'N#N', 'CC(=O)O',
                 'CC(=O)[O-]', '[NH4+]', 'C[N+](=O)[O-]', 'C#N',
                 'c1ccccc1', 'c1ccncc1', 'c1cc[nH]c1', 'c1ccc2ccccc2c1',
                 'CC(=O)NC', 'CS(C)=O', 'CS(=O)(=O)C', 'P(=O)(O)(O)O',
                 'B(O)(O)O', '[B-](F)(F)(F)F', 'C[N+](C)(C)C', 'ClCCl']
        for smiles in cases:
            with self.subTest(smiles=smiles):
                reference = Chem.AddHs(Chem.MolFromSmiles(smiles))
                symbols = [a.GetSymbol() for a in reference.GetAtoms()]
                edges = [(b.GetBeginAtomIdx(), b.GetEndAtomIdx()) for b in reference.GetBonds()]
                charge = Chem.GetFormalCharge(reference)
                orders, charges = perceive_bond_orders(symbols, edges, charge=charge)
                mol = Chem.RWMol(reference)
                for atom, q in zip(mol.GetAtoms(), charges):
                    atom.SetIsAromatic(False)
                    atom.SetFormalCharge(q)
                for bond, order in zip(mol.GetBonds(), orders):
                    bond.SetIsAromatic(False)
                    bond.SetBondType([None, Chem.BondType.SINGLE, Chem.BondType.DOUBLE, Chem.BondType.TRIPLE][order])
                Chem.SanitizeMol(mol)
                self.assertEqual(Chem.MolToSmiles(mol), Chem.MolToSmiles(reference))
                # Run RDKit's actual order assignment on the same single graph.
                rdkit_mol = Chem.RWMol(reference)
                for atom in rdkit_mol.GetAtoms():
                    atom.SetIsAromatic(False)
                    atom.SetFormalCharge(0)
                for bond in rdkit_mol.GetBonds():
                    bond.SetIsAromatic(False)
                    bond.SetBondType(Chem.BondType.SINGLE)
                rdDetermineBonds.DetermineBondOrders(rdkit_mol, charge=charge, embedChiral=False)
                Chem.SanitizeMol(rdkit_mol)
                if smiles == 'CS(C)=O':
                    # Our table allows neutral tetravalent sulfur. RDKit's
                    # DetermineBonds table instead selects S+/O- here.
                    self.assertEqual([a.GetFormalCharge() for a in rdkit_mol.GetAtoms()][:4], [0, 1, 0, -1])
                    self.assertEqual(charges[:4], [0, 0, 0, 0])
                else:
                    self.assertEqual(Chem.MolToSmiles(mol), Chem.MolToSmiles(rdkit_mol))
