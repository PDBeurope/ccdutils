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
from pdb_ccd_mogul import PdbCCDMogul
from utilities import cif_filename
from nose.tools import assert_equals


def test_eoh_mogul():
    eoh = PdbCCDMogul(file_name=cif_filename('EOH'))
    # Mogul result on ethanol two analyzed bonds, 1 angle
    assert_equals(len(eoh.analysed_bonds), 2)
    assert_equals(len(eoh.analysed_angles), 1)
    assert_equals(len(eoh.analysed_torsions), 0)
    assert_equals(len(eoh.analysed_rings), 0)



class DummyTestCaseSoPycharmRecognizesNoseTestsAsTests(unittest.TestCase):
    pass
