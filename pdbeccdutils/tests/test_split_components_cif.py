import os
import unittest

from nose.tools import assert_true

from pdbeccdutils.pdb_chemical_components_rdkit import PdbChemicalComponentsRDKit
from pdbeccdutils.split_components_cif import SplitComponentsCif
from pdbeccdutils.utilities import file_name_in_tsts_out, test_components_cif_first_file_comps


class TestSplitComponentsCif(unittest.TestCase):
    def test_loading_test_file(self):
        test_file = test_components_cif_first_file_comps
        if not os.path.isfile(test_file):
            self.fail('test file {} with 5 ccd in it does not exist?')
        split_cc = SplitComponentsCif(test_file)
        # cif dictionary object should have 5 ccds
        self.assertEqual(len(split_cc.cif_dictionary), 5)

    def test_individual_cif_dictionary_blockid(self):
        split_cc = SplitComponentsCif(test_components_cif_first_file_comps)
        block_ids = []
        for individual_dict in split_cc.individual_cif_dictionary():
            self.assertEqual(len(individual_dict), 1)
            block_id = list(individual_dict)[0]
            block_ids.append(block_id)
        self.assertEqual(['000', '001', '002', '003', '004'], block_ids)

    def test_individual_cif_dictionary_with_rdkit_load(self):
        split_cc = SplitComponentsCif(test_components_cif_first_file_comps)
        chem_comp_id_s = []
        for individual_dict in split_cc.individual_cif_dictionary():
            pdb_cc_rdkit = PdbChemicalComponentsRDKit(cif_dictionary=individual_dict)
            chem_comp_id = pdb_cc_rdkit.chem_comp_id
            chem_comp_id_s.append(chem_comp_id)
        self.assertEqual(['000', '001', '002', '003', '004'], chem_comp_id_s)

    def test_loading_file_that_does_not_exist_raises_ioerror(self):
        with self.assertRaises(IOError):
            SplitComponentsCif('/////impossible')


def test_write_individual_cif_dictionary():
        split_cc = SplitComponentsCif(test_components_cif_first_file_comps)
        for individual_dict in split_cc.individual_cif_dictionary():
            block_id = list(individual_dict)[0]
            file_name = file_name_in_tsts_out(block_id + '_split.cif')
            SplitComponentsCif.write_individual_cif_dictionary(individual_dict, file_name)
            yield assert_true, os.path.isfile(file_name), 'individual cif dictionary {} ' \
                                                          'must be written'.format(file_name)


if __name__ == '__main__':
    unittest.main()
