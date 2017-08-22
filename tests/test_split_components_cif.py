import os
import unittest
from split_components_cif import SplitComponentsCif
from utilities import test_components_cif_first_file_comps

class TestSplitComponentsCif(unittest.TestCase):
    def test_loading_test_file(self):
        test_file = test_components_cif_first_file_comps
        if not os.path.isfile(test_file):
            self.fail('test file {} with 5 ccd in it does not exist?')
        split_cc = SplitComponentsCif(test_file)
        # cif dictionary object should have 5 ccds
        self.assertEqual(len(split_cc.cif_dictionary), 5)

    def test_individual_pdb_ccd_rdkit(self):
        split_cc = SplitComponentsCif(test_components_cif_first_file_comps)
        chem_comp_id_s = []
        for pdb_cc_rdkit in split_cc.individual_pdb_ccd_rdkit():
            chem_comp_id = pdb_cc_rdkit.chem_comp_id
            chem_comp_id_s.append(chem_comp_id)
        self.assertEqual(['000', '001', '002', '003', '004'],chem_comp_id_s)


if __name__ == '__main__':
    unittest.main()
