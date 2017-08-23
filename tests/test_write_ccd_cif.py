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

from nose.tools import assert_equals, assert_true

from pdb_chemical_components import PdbChemicalComponents
from utilities import cif_filename, file_name_in_tsts_out

# currently limited to PDBeCIF!
cif_parser_list = ('PDBeCIF',)

def test_load_eoh_from_cif_write_it_out():
    for cif_parser in cif_parser_list:
        description = ', with cif_parser={}'.format(cif_parser)
        try:
            orig_pdbccd = PdbChemicalComponents(file_name=cif_filename('EOH'), cif_parser=cif_parser)
            out_cif_file_name = file_name_in_tsts_out('EOH_test_write.cif')
            orig_pdbccd.write_ccd_cif(out_cif_file_name)
            yield assert_true, os.path.isfile(out_cif_file_name), 'have written cif file {}'.format(out_cif_file_name)
            read_back_pdb_ccd = PdbChemicalComponents(file_name=out_cif_file_name, cif_parser=cif_parser)
            yield assert_equals, 'EOH', read_back_pdb_ccd.chem_comp_id, 'readback chem_comp_id' + description
            yield assert_equals, 'LFQSCWFLJHTTHZ-UHFFFAOYSA-N', read_back_pdb_ccd.inchikey, 'readback inchikey'
            yield assert_equals, orig_pdbccd, read_back_pdb_ccd, 'check equality of orig_pdbccd and read_back_pdb_ccd'
        except ImportError:
            pass

class DummyTestCaseSoPycharmRecognizesNoseTestsAsTests(unittest.TestCase):
    pass
