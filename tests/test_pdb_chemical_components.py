#!/usr/bin/env python
# software from PDBe: Protein Data Bank in Europe; http://pdbe.org
#
# Copyright 2017 EMBL - European Bioinformatics Institute
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing,
# software distributed under the License is distributed on
# an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
# KIND, either express or implied. See the License for the
# specific language governing permissions and limitations
# under the License.
#
import os
import unittest
from os.path import expanduser

from nose.tools import assert_equals, assert_true

from pdb_chemical_components import PdbChemicalComponents
from utilities import cif_filename, test_cif_path_name

# by default just use the 'auto' cif_parse (which ever can be imported).
cif_parser_list = ('auto',)
# but if environment variable TEST_CIF_PARSER_ALL test all possible
if 'TEST_CIF_PARSER_ALL' in os.environ:
    cif_parser_list = ('PDBeCIF', 'CifFile')


def test_check_test_file_path_name_is_accessible_after_changing_dir():
    home = expanduser("~")
    os.chdir(home)
    assert_true(os.path.isdir(test_cif_path_name))


def test_hard_code_cmo():
    cmo = PdbChemicalComponents(cif_parser='test_hard_code_cmo')
    yield assert_equals, 'CMO', cmo.chem_comp_id, 'chem_comp_id'
    yield assert_equals, 'CARBON MONOXIDE', cmo.chem_comp_name, 'chem_comp_name'
    yield assert_equals, 'REL', cmo.chem_comp_pdbx_release_status, 'chem_comp_pdbx_release_status'
    yield assert_equals, 'UGFAIRIUMAVXCW-UHFFFAOYSA-N', cmo.inchikey, 'chem_inchikey'
    yield assert_equals, 2, cmo.number_atoms, 'number_atoms'
    yield assert_equals, ('C', 'O'), cmo.atom_ids, 'atom_ids'
    yield assert_equals, ('C', 'O'), cmo.atom_elements, 'atom_elements'
    yield assert_equals, ('N', 'N'), cmo.atom_stereo_configs, 'atom_stereo_configs'
    yield assert_equals, (-1, 1), cmo.atom_charges, '(property) atom_charges'
    yield assert_equals, 1, cmo.number_bonds, 'number_bonds'
    the_bond = cmo.bonds[0]
    yield assert_equals, 'C', the_bond.atom_id_1, 'bond atom_id_1'
    yield assert_equals, 'O', the_bond.atom_id_2, 'bond atom_id_2'
    yield assert_equals, [0], cmo.bond_atom_index_1, '(generated) bond_atom_index_1'
    yield assert_equals, [1], cmo.bond_atom_index_2, '(generated) bond_atom_index_2'
    yield assert_equals, [3], cmo.bond_order, '(generated) bond_order'
    yield assert_equals, [False], cmo.bond_aromatic, '(generated) bond_aromatic'
    yield assert_equals, ((0.607, 0., 0.), (-0.6, 0., 0.)), cmo.ideal_xyz, '(generated) ideal_xyz'
    yield assert_equals, ((-0.296, 8.526, 17.112), (0.023, 7.997, 18.053)), cmo.model_xyz, '(generated) model_xyz'


def test_eoh_loads_with_parser():
    for cif_parser in cif_parser_list:
        try:
            PdbChemicalComponents(file_name=cif_filename('EOH'), cif_parser=cif_parser)
        except ImportError as msg:
            yield assert_equals, 0, 1, \
                "Import problem using cif_parser='{}' message='{}'. " \
                "Other tests will skip using testing this parser.".format(cif_parser, msg)


def test_load_eoh_from_cif():
    for cif_parser in cif_parser_list:
        description = ', with cif_parser={}'.format(cif_parser)
        try:
            eoh = PdbChemicalComponents(file_name=cif_filename('EOH'), cif_parser=cif_parser)
            yield assert_equals, 'EOH', eoh.chem_comp_id, 'chem_comp_id' + description
            yield assert_equals, 'ETHANOL', eoh.chem_comp_name, 'chem_comp_name' + description
            yield assert_equals, 'REL', eoh.chem_comp_pdbx_release_status, 'chem_comp_pdbx_release_status'
            yield assert_equals, 'LFQSCWFLJHTTHZ-UHFFFAOYSA-N', eoh.inchikey, 'inchikey' + description
            yield assert_equals, 9, eoh.number_atoms, 'number_atoms' + description
            yield assert_equals, ('C1', 'C2',  'O', 'H11', 'H12', 'H21', 'H22', 'H23', 'HO'), \
                eoh.atom_ids, 'atom_ids' + description
            yield assert_equals, ('C', 'C',  'O', 'H', 'H', 'H', 'H', 'H', 'H'), \
                eoh.atom_elements, 'atom_elements' + description
            yield assert_equals, tuple(['N'] * 9), eoh.atom_stereo_configs, 'atom_stereo_configs' + description
            yield assert_equals, tuple([0] * 9), eoh.atom_charges, '(property) atom_charges'
            yield assert_equals, 8, eoh.number_bonds, 'number_bonds' + description
            third_bond = eoh.bonds[2]
            yield assert_equals, 'C1', third_bond.atom_id_1, 'third bond atom_id_1' + description
            yield assert_equals, 'H11', third_bond.atom_id_2, 'third bond atom_id_2' + description
            yield assert_equals, 0, eoh.bond_atom_index_1[2], 'third bond (generated) bond_atom_index_1' + description
            yield assert_equals, 3, eoh.bond_atom_index_2[2], 'third bond (generated) bond_atom_index_2' + description
            yield assert_equals, [1] * 8, eoh.bond_order,  '(generated) bond_order' + description
            yield assert_equals, [False] * 8, eoh.bond_aromatic,  '(generated) bond_aromatic' + description
            yield assert_equals, (1.130, 0.315, 0.), eoh.ideal_xyz[2], 'property ideal_xyz - check atom #3 the oxygen'
            yield assert_equals, (15.861, 48.185, 8.256), eoh.model_xyz[2], 'property model_xyz - check atom #3 the oxygen'
        except ImportError:
            pass


def test_load_hem_from_cif():
    for cif_parser in cif_parser_list:
        description = ', with cif_parser={}'.format(cif_parser)
        try:
            hem_file = os.path.join(test_cif_path_name, 'HEM_edited.cif')
            hem = PdbChemicalComponents(hem_file, cif_parser=cif_parser)
            yield assert_equals, 'HEM', hem.chem_comp_id, 'chem_comp_id' + description
            yield assert_equals, True, 'Fe' in hem.atom_elements, 'Fe in hem.atom_elements' + description
        except ImportError:
            pass


class DummyTestCaseSoPycharmRecognizesNoseTestsAsTests(unittest.TestCase):
    pass
