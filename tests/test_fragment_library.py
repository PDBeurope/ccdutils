import unittest
from fragment_library import FragmentLibrary
from utilities import cif_filename
from rdkit import Chem
from pdb_chemical_components_rdkit import PdbChemicalComponentsRDKit


class TestFragmentLibrary(unittest.TestCase):
    def setUp(self):
        self.frag = FragmentLibrary()

    def test_number_of_fragments_not_zero(self):
        self.assertNotEqual(0, self.frag.number_of_entries)

    def test_phenyl_in_fragment(self):
        self.assertIn('phenyl', self.frag.smiles_to_fragment_name.values())
        self.assertEqual('phenyl', self.frag.smiles_to_fragment_name['c1ccccc1'])

    def test_rdkit_mol_smiles_for_phenyl(self):
        frag = FragmentLibrary()
        rdkit_mol_phenyl = frag.smiles_to_rdkit_molecule['c1ccccc1']
        self.assertEqual('c1ccccc1', Chem.MolToSmiles(rdkit_mol_phenyl))

    def test_fragments_for_glu(self):
        cif_file = cif_filename('GLU')
        pdb_cc_rdkit = PdbChemicalComponentsRDKit(file_name=cif_file)
        fragments = self.frag.fragments_for_pdb_chemical_components_rdkit(pdb_cc_rdkit)
        self.assertIn('peptide', fragments)
        self.assertEqual(1, len(fragments['peptide']))  # there is one peptide in GLU
        self.assertEquals(sorted(fragments['peptide'][0]), sorted(['O', 'C', 'CA', 'N']))

    def test_fragments_for_007(self):
        cif_file = cif_filename('007')
        pdb_cc_rdkit = PdbChemicalComponentsRDKit(file_name=cif_file)
        fragments = self.frag.fragments_for_pdb_chemical_components_rdkit(pdb_cc_rdkit)
        self.assertIn('phenyl', fragments)
        self.assertEqual(1, len(fragments['phenyl']))
        self.assertEquals(sorted(fragments['phenyl'][0]), sorted(['C6', 'C7', 'C8', 'C9', 'C10', 'C11']))
        self.assertIn('cyclopentane', fragments)
        self.assertEqual(1, len(fragments['cyclopentane']))
        self.assertEquals(sorted(fragments['cyclopentane'][0]), sorted(['C1', 'C2', 'C3', 'C4', 'C5']))

    def test_fragments_for_bcd(self):
        cif_file = cif_filename('BCD')
        pdb_cc_rdkit = PdbChemicalComponentsRDKit(file_name=cif_file)
        fragments = self.frag.fragments_for_pdb_chemical_components_rdkit(pdb_cc_rdkit)
        self.assertIn('pyranose', fragments)
        self.assertEqual(7, len(fragments['pyranose']))  # there are 7 pyranose rings in BCD

    def test_fragments_for_atp(self):
        cif_file = cif_filename('ATP')
        pdb_cc_rdkit = PdbChemicalComponentsRDKit(file_name=cif_file)
        fragments = self.frag.fragments_for_pdb_chemical_components_rdkit(pdb_cc_rdkit)
        self.assertIn('adenine', fragments)
        self.assertEqual(1, len(fragments['adenine']))
        self.assertEquals(sorted(fragments['adenine'][0]),
                          sorted(['N1', 'C2', 'N3', 'C4', 'C5', 'C6', 'N6', 'N7', 'C8', 'N9']))

if __name__ == '__main__':
    unittest.main()
