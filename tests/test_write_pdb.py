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

from nose.tools import assert_true
from pdb_chemical_components_rdkit import PdbChemicalComponentsRDKit
from utilities import supply_list_of_sample_cifs, file_name_in_tsts_out


def test_pdb_write_for_all_sample_cifs():
    for ciffile in supply_list_of_sample_cifs():
        pdb_cc = PdbChemicalComponentsRDKit(file_name=ciffile)
        pdb_ideal_with_h = file_name_in_tsts_out(pdb_cc.chem_comp_id + '.ideal_withH.pdb')
        pdb_cc.pdb_file_or_string(file_name=pdb_ideal_with_h, ideal=True)
        yield assert_true, os.path.isfile(pdb_ideal_with_h) and os.path.getsize(pdb_ideal_with_h) > 0, \
            '{} call to pdb_cc.sdf_file_or_string(file="{}") must create a non-empty file.'.\
            format(pdb_cc.chem_comp_id, pdb_ideal_with_h)
        pdb_model_with_h = file_name_in_tsts_out(pdb_cc.chem_comp_id + '.model_withH.pdb')
        pdb_cc.pdb_file_or_string(file_name=pdb_model_with_h, ideal=False)
        yield assert_true, os.path.isfile(pdb_model_with_h) and os.path.getsize(pdb_model_with_h) > 0, \
            '{} call to pdb_cc.sdf_file_or_string(file="{}") must create a non-empty file.'.\
            format(pdb_cc.chem_comp_id, pdb_model_with_h)


class DummyTestCaseSoPycharmRecognizesNoseTestsAsTests(unittest.TestCase):
    pass
