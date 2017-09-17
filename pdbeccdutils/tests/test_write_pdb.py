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

from nose.tools import assert_true, assert_equals
from pdb_chemical_components_rdkit import PdbChemicalComponentsRDKit
from utilities import supply_list_of_sample_cifs, file_name_in_tsts_out, cif_filename, test_comparison_files_path


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

def test_XXX_atom_records_against_previous_pdbechem():
    for three_letter_code, ideal in {'ATP': True, '00O': False}.items():
        yield assert_true, XXX_atom_records_against_previous_pdbechem(three_letter_code, ideal), \
              'check atom record against previous for {} ideal {}'.format(three_letter_code, ideal)

def XXX_atom_records_against_previous_pdbechem( three_letter_code, ideal):
    """
    compare ccd PDB HETATM records (ideal coordinates) with ATOM from the previous PDBeChem pdb file for ATP
    """
    ciffile = cif_filename(three_letter_code)
    pdb_cc = PdbChemicalComponentsRDKit(file_name=ciffile)
    pdb_string = pdb_cc.pdb_file_or_string(ideal=ideal)
    lines = pdb_string.split('\n')
    lines_hetatm = list(filter(lambda x: x.startswith('HETATM'), lines))
    comparison_file = os.path.join(test_comparison_files_path, three_letter_code + '.pdb')
    with open(comparison_file, 'r') as comp_file:
        lines = comp_file.read().splitlines()
    lines_atom = list(filter(lambda x: x.startswith('ATOM'), lines))
    for line_no in range(len(lines_atom)):
        line = lines_atom[line_no]
        line = line.replace('ATOM  ', 'HETATM')  # comparison file used ATOM rather than HETATM
        line = line.replace('    0  ', 'A   1  ')  # chain id and residue number have changed to A 1
        line = line.replace('+0', '  ')  # charge at end of line (incorrect in prvious)
        lines_atom[line_no] = line
    assert_equals(len(lines_hetatm), len(lines_atom))  # number HETATM/ATOM line ccd_utils/comparison file ==
    for line_no in range(len(lines_hetatm)):
        assert_equals(lines_hetatm[line_no], lines_atom[line_no])  # HETATM/ATOM line cf ccd_utils/comparison file
    return True


class DummyTestCaseSoPycharmRecognizesNoseTestsAsTests(unittest.TestCase):
    pass
