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
from nose.tools import assert_equals
import os
from pdb_chemical_components import PdbChemicalComponents


def load_carbon_monoxide_hard_coded():
    cmo = PdbChemicalComponents()
    cmo.load_carbon_monoxide_hard_coded()
    return cmo


def test_hard_code_cmo_chem_comp_id_is_cmo():
    cmo = load_carbon_monoxide_hard_coded()
    assert_equals('CMO', cmo.chem_comp_id)


def test_hard_code_cmo_chem_comp_name():
    cmo = load_carbon_monoxide_hard_coded()
    assert_equals('CARBON MONOXIDE', cmo.chem_comp_name)


def test_hard_code_cmo_inchikey():
    cmo = load_carbon_monoxide_hard_coded()
    assert_equals('UGFAIRIUMAVXCW-UHFFFAOYSA-N', cmo.inchikey)


def test_hard_code_has_2_atoms():
    cmo = load_carbon_monoxide_hard_coded()
    assert_equals(2, cmo.number_atoms)


def test_hard_code_cmo_atom_ids():
    cmo = load_carbon_monoxide_hard_coded()
    assert_equals(('C', 'O'), cmo.atom_ids)


def test_hard_code_has_1_bond():
    cmo = load_carbon_monoxide_hard_coded()
    assert_equals(1, cmo.number_bonds)


def test_hard_code_bond_atom_ids():
    cmo = load_carbon_monoxide_hard_coded()
    the_bond = cmo.bonds[0]
    assert_equals('C', the_bond.atom_id_1)
    assert_equals('O', the_bond.atom_id_2)


def eoh_cif_file():
    return os.path.join('data', 'cif', 'EOH.cif')


def test_ciffile_parser_loads():
    try:
        PdbChemicalComponents(file_name=eoh_cif_file(), cif_parser='CifFile')
    except ImportError as msg:
        assert False, 'Failed to load module for CifFile parser "{}". ' \
                      'Other tests will skip using testing this parser'.format(msg)


def test_load_eoh_from_cif():
    for cif_parser in 'mmcifIO', 'CifFile':
        try:
            eoh = PdbChemicalComponents(file_name=eoh_cif_file(), cif_parser=cif_parser)
            description = ', with cif_parser={}'.format(cif_parser)
            yield assert_equals, 'EOH', eoh.chem_comp_id, 'chem_comp_id' + description
            yield assert_equals, 'ETHANOL', eoh.chem_comp_name, 'chem_comp_name' + description
            yield assert_equals, 'LFQSCWFLJHTTHZ-UHFFFAOYSA-N', eoh.inchikey, 'inchikey' + description
            yield assert_equals, 9, eoh.number_atoms, 'number_atoms' + description
            yield assert_equals, ('C1', 'C2',  'O', 'H11', 'H12', 'H21', 'H22', 'H23', 'HO'), \
                eoh.atom_ids, 'atom_ids' + description
            yield assert_equals, 8, eoh.number_bonds, 'number_bonds' + description
        except ImportError:
            pass
