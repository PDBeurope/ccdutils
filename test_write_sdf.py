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
import glob
import os
import unittest

from nose.tools import assert_equals, assert_in, assert_not_equal, assert_is_instance, assert_true, assert_not_in
from test_pdb_chemical_components import cif_filename
from pdb_chemical_components_rdkit import PdbChemicalComponentsRDKit
from test_pdb_ccd_rdkit_loading import supply_list_of_sample_cifs


def test_hard_code_cmo():
    cmo = PdbChemicalComponentsRDKit(cif_parser='test_hard_code_cmo')
    sdf_string_ideal = cmo.sdf_file_or_string(ideal=True)
    yield assert_is_instance, sdf_string_ideal, str, 'sdf string must be a string'
    if isinstance(sdf_string_ideal, str):
        yield assert_not_equal, 0, len(sdf_string_ideal), \
            'zero characters in sdf string ideal={}'.format(sdf_string_ideal)
        yield assert_in, 'V2000', sdf_string_ideal, 'the sdf string ideal should contain V2000'
        yield assert_in, 'END', sdf_string_ideal, 'the sdf string ideal should contain END'
        yield assert_in, '0.607', sdf_string_ideal, 'the sdf string ideal should contain x coordinate for C'
    sdf_string_model = cmo.sdf_file_or_string(ideal=False)
    yield assert_in, '-0.296', sdf_string_model, 'the sdf string model should contain x coordinate for C'
    sdf_ideal_with_h = file_name_in_subdir_for_output_files('CMO.hard_coded.ideal_withH.sdf')
    sdf_model_with_h = file_name_in_subdir_for_output_files('CMO.hard_coded.model_withH.sdf')
    cmo.sdf_file_or_string(file_name=sdf_ideal_with_h)
    cmo.sdf_file_or_string(file_name=sdf_model_with_h,ideal=False)
    yield assert_true, os.path.isfile(sdf_ideal_with_h) and os.path.getsize(sdf_ideal_with_h) > 0, \
        'call to cmo.sdf_file_or_string(file="{}") must create a non-empty file.'.format(sdf_ideal_with_h)
    yield assert_true, os.path.isfile(sdf_model_with_h) and os.path.getsize(sdf_model_with_h) > 0, \
        'call to cmo.sdf_file_or_string(file="{}") must create a non-empty file.'.format(sdf_model_with_h)

def test_load_eoh_from_cif():
    eoh = PdbChemicalComponentsRDKit(file_name=cif_filename('EOH'))
    sdf_string_ideal_h = eoh.sdf_file_or_string(ideal = True)
    sdf_string_model_h = eoh.sdf_file_or_string(ideal = False)
    sdf_string_ideal_no_h = eoh.sdf_file_or_string(hydrogen = False)
    sdf_string_model_no_h = eoh.sdf_file_or_string(ideal = False, hydrogen = False)
    yield assert_true, len(sdf_string_ideal_no_h) > 0, 'sdf_file_or_string must create a non-empty str'
    yield assert_true, len(sdf_string_model_no_h) > 0, 'sdf_file_or_string must create a non-empty str'
    yield assert_not_in, ' H ', sdf_string_ideal_no_h, 'sdf_file_or_string must create a non-empty str without H atom'
    yield assert_not_in, ' H ', sdf_string_model_no_h, 'sdf_file_or_string must create a non-empty str without H atom'
    yield assert_true, sdf_string_ideal_h.startswith('EOH'), 'ideal_h: sdf_file_or_string must start with EOH'
    yield assert_true, sdf_string_model_h.startswith('EOH'), 'model_h: sdf_file_or_string must start with EOH'
    yield assert_true, sdf_string_model_no_h.startswith('EOH'), 'model_no_h: sdf_file_or_string must start with EOH'
    yield assert_true, sdf_string_ideal_no_h.startswith('EOH'), 'ideal_no_h: sdf_file_or_string must start with EOH'
    # alias stuff
    yield assert_not_in, 'H23', sdf_string_ideal_h, 'sdf_string_ideal_h should not have atom alias for H23'
    sdf_string_ideal_h_with_alias = eoh.sdf_file_or_string(ideal=True, alias=True)
    yield assert_in,  'H23', sdf_string_ideal_h_with_alias, \
        'sdf produce with alias=False should have an atom alias record for H23'


def test_sdf_write_for_all_sample_cifs():
    for ciffile in supply_list_of_sample_cifs():
        #if 'HEM' in ciffile:
        #   continue
        pdb_cc = PdbChemicalComponentsRDKit(file_name=ciffile)
        sdf_ideal_with_h = file_name_in_subdir_for_output_files(pdb_cc.chem_comp_id + '.ideal_withH.sdf')
        pdb_cc.sdf_file_or_string(file_name=sdf_ideal_with_h)
        yield assert_true, os.path.isfile(sdf_ideal_with_h) and os.path.getsize(sdf_ideal_with_h) > 0, \
            '{} call to pdb_cc.sdf_file_or_string(file="{}") must create a non-empty file.'.\
            format(pdb_cc.chem_comp_id, sdf_ideal_with_h)
        sdf_model_with_h = file_name_in_subdir_for_output_files(pdb_cc.chem_comp_id + '.model_withH.sdf')
        pdb_cc.sdf_file_or_string(file_name = sdf_model_with_h,ideal=False)
        yield assert_true, os.path.isfile(sdf_model_with_h) and os.path.getsize(sdf_model_with_h) > 0, \
            '{} call to pdb_cc.sdf_file_or_string(file="{}") must create a non-empty file.'.\
            format(pdb_cc.chem_comp_id, sdf_model_with_h)
        sdf_ideal_no_h = file_name_in_subdir_for_output_files(pdb_cc.chem_comp_id + '.ideal_noH.sdf')
        pdb_cc.sdf_file_or_string(file_name = sdf_ideal_no_h, hydrogen = False)
        yield assert_true, os.path.isfile(sdf_model_with_h) and os.path.getsize(sdf_ideal_no_h) > 0, \
            '{} call to pdb_cc.sdf_file_or_string(file="{}") must create a non-empty file.'.\
            format(pdb_cc.chem_comp_id, sdf_ideal_no_h)
        sdf_model_no_h = file_name_in_subdir_for_output_files(pdb_cc.chem_comp_id + '.model_noH.sdf')
        pdb_cc.sdf_file_or_string(file_name = sdf_model_no_h, ideal = False, hydrogen = False)
        yield assert_true, os.path.isfile(sdf_model_with_h) and os.path.getsize(sdf_model_no_h) > 0, \
            '{} call to pdb_cc.sdf_file_or_string(file="{}") must create a non-empty file.'.\
            format(pdb_cc.chem_comp_id, sdf_model_no_h)


def file_name_in_subdir_for_output_files(file_name):
    """
    creates the subdirectory "tests_out" if necessary and returns the file_name in this directory, If the
    file already exists it will remove it.

    Args:
        file_name (str):  the name for the file

    Returns:
        str: the filename in the subdirectory tests_out
    """
    subdir = 'tests_out'
    if not os.path.isdir(subdir):
        os.mkdir(subdir)
    out_file_name = os.path.join(subdir, file_name)
    if os.path.isfile(out_file_name):
        os.remove(out_file_name)
    return out_file_name

class DummyTestCaseSoPycharmRecognizesNoseTestsAsTests(unittest.TestCase):
    pass
