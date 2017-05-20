#!/usr/bin/env python
import os
import unittest
from pdb_chemical_components import PdbChemicalComponents


class MyTestCase(unittest.TestCase):

    def load_carbon_monoxide_hard_coded(self):
        self.cmo = PdbChemicalComponents()
        self.cmo.load_carbon_monoxide_hard_coded()

    def test_hard_code_CMO_chem_comp_id_is_CMO(self):
        self.load_carbon_monoxide_hard_coded()
        self.assertEqual('CMO', self.cmo.chem_comp_id)

    def test_hard_code_CMO_chem_comp_name(self):
        self.load_carbon_monoxide_hard_coded()
        self.assertEqual('CARBON MONOXIDE', self.cmo.chem_comp_name)

    def test_hard_code_CMO_inchikey(self):
        self.load_carbon_monoxide_hard_coded()
        self.assertEqual('UGFAIRIUMAVXCW-UHFFFAOYSA-N', self.cmo.inchikey)

    def test_hard_code_CMO_atom_ids(self):
        self.load_carbon_monoxide_hard_coded()
        self.assertEqual(('C', 'O'), self.cmo.atom_ids)

    def test_hard_code_CMO_has_2_atoms(self):
        self.load_carbon_monoxide_hard_coded()
        self.assertEqual(2, self.cmo.number_atoms)

    def test_hard_code_CMO_has_1_bond(self):
        self.load_carbon_monoxide_hard_coded()
        self.assertEqual(1, self.cmo.number_bonds)

    def test_hard_code_CMO_bond_atom_ids(self):
        self.load_carbon_monoxide_hard_coded()
        the_bond = self.cmo.bonds[0]
        self.assertEqual('C', the_bond.atom_id_1)
        self.assertEqual('O', the_bond.atom_id_2)


    def load_EOH_cif(self):
        eoh_cif_file = os.path.join('data', 'cif', 'EOH.cif')
        try:
            self.ccd_eoh_from_ciffile = PdbChemicalComponents(file_name=eoh_cif_file, cif_parser='CifFile')
        except ImportError:
            self.ccd_eoh_from_ciffile = None
        try:
            self.ccd_eoh_from_mmccifio = PdbChemicalComponents(file_name=eoh_cif_file, cif_parser='mmcifIO')
        except ImportError:
            self.ccd_eoh_from_mmccifio = None

    def test_CifFile_parser_loads(self):
        eoh_cif_file = os.path.join('data', 'cif', 'EOH.cif')
        try:
            PdbChemicalComponents(file_name=eoh_cif_file, cif_parser='CifFile')
        except ImportError:
            self.fail('failed to load module for CifFile parser other tests will not test this.')

    def test_mmcifIO_parser_loads(self):
        eoh_cif_file = os.path.join('data', 'cif', 'EOH.cif')
        try:
            PdbChemicalComponents(file_name=eoh_cif_file, cif_parser='mmcifIO')
        except ImportError:
            self.fail('failed to load module for mmcifIO parser other tests will not test this.')

    def test_load_EOH_cif_chem_comp_id_is_EOH(self):
        self.load_EOH_cif()
        for eoh_from_cif in self.ccd_eoh_from_ciffile, self.ccd_eoh_from_mmccifio:
            if eoh_from_cif is not None:
                self.assertEqual('EOH', eoh_from_cif.chem_comp_id)

    def test_load_EOH_cif_chem_comp_name(self):
        self.load_EOH_cif()
        for eoh_from_cif in self.ccd_eoh_from_ciffile, self.ccd_eoh_from_mmccifio:
            if eoh_from_cif is not None:
                self.assertEqual('ETHANOL', eoh_from_cif.chem_comp_name)

    def test_load_EOH_cif_atom_ids_are_correct(self):
        self.load_EOH_cif()
        for eoh_from_cif in self.ccd_eoh_from_ciffile, self.ccd_eoh_from_mmccifio:
            if eoh_from_cif is not None:
                self.assertEqual(('C1', 'C2',  'O', 'H11', 'H12', 'H21', 'H22', 'H23', 'HO'), eoh_from_cif.atom_ids)

    def test_load_EOH_number_atoms_is_9(self):
        self.load_EOH_cif()
        for eoh_from_cif in self.ccd_eoh_from_ciffile, self.ccd_eoh_from_mmccifio:
            if eoh_from_cif is not None:
                self.assertEqual(9, eoh_from_cif.number_atoms)

if __name__ == '__main__':
    unittest.main()
