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
import traceback
import unittest

from nose.tools import assert_equals, assert_true, assert_false

from pdbeccdutils.utilities import cif_filename, supply_list_of_sample_cifs, file_name_in_tsts_out

try:
    from pdb_ccd_mogul import PdbCCDMogul


    def test_eoh_mogul():
        eoh = PdbCCDMogul(file_name=cif_filename('EOH'))
        eoh.run_mogul()
        # Mogul result on ethanol two analyzed bonds, 1 angle
        assert_equals(len(eoh.store_bonds), 2)
        assert_equals(len(eoh.store_angles), 1)
        assert_equals(len(eoh.store_torsions), 0)
        assert_equals(len(eoh.store_rings), 0)


    def test_html_write_for_all_sample_cifs():
        for cif_file in supply_list_of_sample_cifs():
            pdb_ccd_mog = PdbCCDMogul(file_name=cif_file)
            chem_comp_id = pdb_ccd_mog.pdb_ccd_rdkit.chem_comp_id
            html_out_file = file_name_in_tsts_out(chem_comp_id + '.pdb_ccd_ideal_mogul_report.html')
            try:
                pdb_ccd_mog.run_mogul()
                pdb_ccd_mog.prepare_file_html(html_out_file)
                yield assert_true, os.path.isfile(html_out_file) and os.path.getsize(html_out_file) > 0, \
                    '{} call to pdb_ccd_mogul.prepare_file_html("{}") must create a non-empty file.'.\
                    format(chem_comp_id, html_out_file)
            except Exception as e_mess:
                yield assert_false, 'exception {}\ntraceback{}'.format(e_mess, traceback.format_exc())
except ImportError:

    def test_import_error_skip_test():
        yield assert_true, 'skipping test_pdb_ccd_mogul.py tests because of ImportError (for ccdc)'


class DummyTestCaseSoPycharmRecognizesNoseTestsAsTests(unittest.TestCase):
    pass
