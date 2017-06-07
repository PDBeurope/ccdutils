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
from test_pdb_chemical_components import cif_filename, test_file_path_name
from pdb_chemical_components_rdkit import PdbChemicalComponentsRDKit


def test_load_eoh_from_cif():
    eoh = PdbChemicalComponentsRDKit(file_name=cif_filename('EOH'))
    yield assert_equals, 'EOH', eoh.chem_comp_id, 'chem_comp_id'
    yield assert_equals, eoh.inchikey, eoh.inchikey_from_rdkit, 'inchikey from cif file should match the rdkit inchikey'
    pdb_string_ideal_h = eoh.pdb_file_or_string(ideal = True)
    pdb_string_model_h = eoh.sdf_file_or_string(ideal = False)
    yield assert_in, '0.007', pdb_string_ideal_h, 'pdb_file_or_string must contain x coordinate for ideal EOH'
    yield assert_in, '15.2120', pdb_string_model_h, 'pdb_file_or_string must contain x coordinate for model EOH'
    yield assert_in, 'HEADER', pdb_string_ideal_h, 'pdb_file_or_string must contain title section'
    img_with_h = file_name_in_subdir_for_output_files('EOH.img_withH.png')
    eoh.image_file(file_name=img_with_h, hydrogen=True)
    yield assert_true, os.path.isfile(img_with_h) and os.path.getsize(img_with_h) > 0, \
        '{} call to eoh.image_file(file="{}") must create a non-empty file.'.\
        format('EOH', img_with_h)

def test_inchikey_match_for_all_sample_cifs():
    for ciffile in supply_list_of_sample_cifs():
        pdb_cc = PdbChemicalComponentsRDKit(file_name=ciffile)
        yield assert_equals, pdb_cc.inchikey, pdb_cc.inchikey_from_rdkit, \
            'check inchikeys match for ' + pdb_cc.chem_comp_id

def test_sdf_write_for_all_sample_cifs():
    for ciffile in supply_list_of_sample_cifs():
        #if 'HEM' in ciffile:
        #   continue
        pdb_cc = PdbChemicalComponentsRDKit(file_name=ciffile)
        pdb_ideal_with_h = file_name_in_subdir_for_output_files(pdb_cc.chem_comp_id + '.ideal_withH.pdb')
        pdb_cc.pdb_file_or_string(file_name = pdb_ideal_with_h, ideal = True)
        yield assert_true, os.path.isfile(pdb_ideal_with_h) and os.path.getsize(pdb_ideal_with_h) > 0, \
            '{} call to pdb_cc.sdf_file_or_string(file="{}") must create a non-empty file.'.\
            format(pdb_cc.chem_comp_id, pdb_ideal_with_h)
        pdb_model_with_h = file_name_in_subdir_for_output_files(pdb_cc.chem_comp_id + '.model_withH.pdb')
        pdb_cc.pdb_file_or_string(file_name = pdb_model_with_h, ideal = False)
        yield assert_true, os.path.isfile(pdb_model_with_h) and os.path.getsize(pdb_model_with_h) > 0, \
            '{} call to pdb_cc.sdf_file_or_string(file="{}") must create a non-empty file.'.\
            format(pdb_cc.chem_comp_id, pdb_model_with_h)

def supply_list_of_sample_cifs():
    """
    returns the list of sample pdb ccd cifs for test.

    Args:
        None

    Returns:
        list of filenames
    """
    return sorted(glob.glob(os.path.join(test_file_path_name, '*.cif')))


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
