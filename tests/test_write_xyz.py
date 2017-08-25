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

from nose.tools import assert_equals, assert_in
from pdb_chemical_components_rdkit import PdbChemicalComponentsRDKit
from utilities import cif_filename


def test_xyz_for_eoh_ideal():
    eoh = PdbChemicalComponentsRDKit(file_name=cif_filename('EOH'))
    xyz_string_ideal = eoh.xyz_file_or_string(ideal=True)
    lines = xyz_string_ideal.split('\n')
    yield assert_equals, len(lines), 11, 'EOH xyz format should have 11 lines'
    yield assert_equals, lines[0], 9, 'EOH xyz should start with 9 the number of atoms'
    yield assert_equals, lines[1], 'EOH', 'EOH xyz 2nd line should be "EOH"'
    yield assert_in, 'C', lines[2], 'EOH xyz first atom is C'
    yield assert_in, '0.007', lines[2], 'EOH xyz first atom x coordinate'

def test_xyz_for_eoh_model():
    eoh = PdbChemicalComponentsRDKit(file_name=cif_filename('EOH'))
    xyz_string_model = eoh.xyz_file_or_string(ideal=False)
    lines = xyz_string_model.split('\n')
    yield assert_in, '7.491', lines[2], 'EOH xyz first atom record z coordinate'
