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
import unittest

from nose.tools import assert_equals, assert_not_equals

from pdbeccdutils.pdb_chemical_components_rdkit import PdbChemicalComponentsRDKit
from pdbeccdutils.utilities import cif_filename, supply_list_of_sample_cifs


def test_hard_code_cmo():
    cmo = PdbChemicalComponentsRDKit(cif_parser='test_hard_code_cmo')
    yield assert_equals, 'CMO', cmo.chem_comp_id, 'chem_comp_id'
    yield assert_equals, 'CARBON MONOXIDE', cmo.chem_comp_name, 'chem_comp_name'
    yield assert_equals, 'UGFAIRIUMAVXCW-UHFFFAOYSA-N', cmo.inchikey, 'chem_inchikey'
    # rdkit should be able to get an inchikey that is the same that read from the CCD cif file
    yield assert_equals, 2, cmo.rwmol_original.GetNumAtoms(), 'rdkit_mol.GetNumAtoms() should give 2'
    yield assert_equals, cmo.inchikey, cmo.inchikey_from_rdkit, 'inchikey from cif file should match the rdkit inchikey'


def test_load_eoh_from_cif():
    eoh = PdbChemicalComponentsRDKit(file_name=cif_filename('EOH'))
    yield assert_equals, 'EOH', eoh.chem_comp_id, 'chem_comp_id'
    yield assert_equals, eoh.inchikey, eoh.inchikey_from_rdkit, 'inchikey from cif file should match the rdkit inchikey'


def test_inchikey_match_for_all_sample_cifs():
    for ciffile in supply_list_of_sample_cifs():
        pdb_cc = PdbChemicalComponentsRDKit(file_name=ciffile)
        inchikey_from_rdkit = pdb_cc.inchikey_from_rdkit
        if pdb_cc.chem_comp_id in ('CDL', 'ASX', '7OM'):
            yield assert_not_equals, pdb_cc.inchikey, inchikey_from_rdkit, \
                'know inchikeys do not match for ' + pdb_cc.chem_comp_id
        else:
            yield assert_equals, pdb_cc.inchikey, inchikey_from_rdkit, \
                'check inchikeys match for ' + pdb_cc.chem_comp_id


class DummyTestCaseSoPycharmRecognizesNoseTestsAsTests(unittest.TestCase):
    pass
