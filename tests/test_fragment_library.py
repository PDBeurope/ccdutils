import unittest
from fragment_library import FragmentLibrary
from utilities import cif_filename
from rdkit import Chem
from pdb_chemical_components_rdkit import PdbChemicalComponentsRDKit


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

    def test_fragments_for_glu(self):
        frag = FragmentLibrary()
        cif_file = cif_filename('GLU')
        pdb_cc_rdkit = PdbChemicalComponentsRDKit(file_name=cif_file)
        fragments = frag.fragments_for_pdb_chemical_components_rdkit(pdb_cc_rdkit)
        self.assertTrue('peptide' in fragments)
        self.assertEquals(fragments['peptide'], ['O', 'C', 'CA', 'N'])


if __name__ == '__main__':
    unittest.main()
