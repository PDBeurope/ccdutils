import os
import unittest
from pdb_chemical_components import PdbChemicalComponents


class MyTestCase(unittest.TestCase):
    def test_load_GLC_number_atoms_is_12(self):
        ccd = PdbChemicalComponents()
        glc_cif_file = os.path.join('data', 'cif', 'GLC.cif')
        ccd.read_ccd_from_cif_file(glc_cif_file)
        self.assertEqual(12, ccd.number_atoms())


if __name__ == '__main__':
    unittest.main()
