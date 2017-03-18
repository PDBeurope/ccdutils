#!/usr/bin/env python
import os
import unittest
from pdb_chemical_components import PdbChemicalComponents


class MyTestCase(unittest.TestCase):

    def test_load_EOH_cif_chem_comp_id_is_EOH(self):
        eoh_cif_file = os.path.join('data', 'cif', 'EOH.cif')
        ccd = PdbChemicalComponents(file_name=eoh_cif_file)
        self.assertEqual('EOH', ccd.chem_comp_id)

    def test_load_EOH_cif_atom_ids_are_correct(self):
        eoh_cif_file = os.path.join('data', 'cif', 'EOH.cif')
        ccd = PdbChemicalComponents(file_name=eoh_cif_file)
        self.assertEqual(['C1', 'C2',  'O', 'H11', 'H12', 'H21', 'H22', 'H23', 'HO'], ccd.chem_comp_id)

    def test_load_EOH_number_atoms_is_9(self):
        eoh_cif_file = os.path.join('data', 'cif', 'EOH.cif')
        ccd = PdbChemicalComponents()
        ccd.read_ccd_from_cif_file(eoh_cif_file)
        self.assertEqual(9, ccd.number_of_atoms)


if __name__ == '__main__':
    unittest.main()
