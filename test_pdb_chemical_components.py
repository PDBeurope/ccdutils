#!/usr/bin/env python
import os
import unittest
from pdb_chemical_components import PdbChemicalComponents


class MyTestCase(unittest.TestCase):

    def setUp(self):
        eoh_cif_file = os.path.join('data', 'cif', 'EOH.cif')
        self.ccd_eoh = PdbChemicalComponents(file_name=eoh_cif_file, cif_parser='CifFile')
        self.ccd_hard_code_cmo = PdbChemicalComponents()
        self.ccd_hard_code_cmo.test_hard_code_CMO()

    def tearDown(self):
        pass

    def test_hard_code_CMO_chem_comp_id_is_CMO(self):
        self.assertEqual('CMO', self.ccd_hard_code_cmo.chem_comp_id)

    def test_hard_code_CMO_chem_comp_name(self):
        self.assertEqual('CARBON MONOXIDE', self.ccd_hard_code_cmo.chem_comp_name)

    def test_hard_code_CMO_atom_ids(self):
        self.assertEqual(('C', 'O'), self.ccd_hard_code_cmo.atom_ids)

    def test_hard_code_CMO_has_2_atoms(self):
        self.assertEqual(2, self.ccd_hard_code_cmo.number_atoms)

    def test_load_EOH_cif_chem_comp_id_is_EOH(self):
        self.assertEqual('EOH', self.ccd_eoh.chem_comp_id)

    def test_load_EOH_cif_chem_comp_name(self):
        self.assertEqual('ETHANOL', self.ccd_eoh.chem_comp_name)

    def test_load_EOH_cif_atom_ids_are_correct(self):
        self.assertEqual(('C1', 'C2',  'O', 'H11', 'H12', 'H21', 'H22', 'H23', 'HO'), self.ccd_eoh.atom_ids)

    def test_load_EOH_number_atoms_is_9(self):
        self.assertEqual(9, self.ccd_eoh.number_atoms)


if __name__ == '__main__':
    unittest.main()
