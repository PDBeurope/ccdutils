import unittest
from fragment_library import FragmentLibrary
from rdkit import Chem


class TestFragmentLibrary(unittest.TestCase):
    def test_number_of_fragments_not_zero(self):
        frag = FragmentLibrary()
        self.assertNotEqual(0, frag.number_of_entries)

    def test_phenyl_in_fragment(self):
        frag = FragmentLibrary()
        self.assertIn('phenyl', frag.smiles_to_fragment_name.values())
        self.assertEqual('phenyl', frag.smiles_to_fragment_name['c1ccccc1'])

    def test_rdkit_mol_smiles_for_phenyl(self):
        frag = FragmentLibrary()
        rdkit_mol_phenyl = frag.smiles_to_rdkit_molecule['c1ccccc1']
        self.assertEqual('c1ccccc1', Chem.MolToSmiles(rdkit_mol_phenyl))


if __name__ == '__main__':
    unittest.main()
