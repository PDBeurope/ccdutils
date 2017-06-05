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

from nose.tools import assert_equals, assert_in, assert_not_equal, assert_is_instance, assert_true
from test_pdb_chemical_components import cif_filename, test_file_path_name
from pdb_chemical_components_rdkit import PdbChemicalComponentsRDKit


def test_hard_code_cmo():
    cmo = PdbChemicalComponentsRDKit(cif_parser='test_hard_code_cmo')
    yield assert_equals, 'CMO', cmo.chem_comp_id, 'chem_comp_id'
    yield assert_equals, 'CARBON MONOXIDE', cmo.chem_comp_name, 'chem_comp_name'
    yield assert_equals, 'UGFAIRIUMAVXCW-UHFFFAOYSA-N', cmo.inchikey, 'chem_inchikey'
    # rdkit should be able to get an inchikey that is the same that read from the CCD cif file
    yield assert_equals, 2, cmo.rdkit_mol.GetNumAtoms(), 'rdkit_mol.GetNumAtoms() should give 2'
    yield assert_equals, cmo.inchikey, cmo.inchikey_from_rdkit, 'inchikey from cif file should match the rdkit inchikey'
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
    yield assert_equals, 'EOH', eoh.chem_comp_id, 'chem_comp_id'
    yield assert_equals, eoh.inchikey, eoh.inchikey_from_rdkit, 'inchikey from cif file should match the rdkit inchikey'
    sdf_string_ideal_no_h = eoh.sdf_file_or_string(hydrogen = False)
    sdf_string_model_no_h = eoh.sdf_file_or_string(ideal = False, hydrogen = False)
    yield assert_true, len(sdf_string_ideal_no_h) > 0, 'sdf_file_or_string must create a non-empty string'
    yield assert_true, len(sdf_string_model_no_h) > 0, 'sdf_file_or_string must create a non-empty string'

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
        sdf_ideal_with_h = file_name_in_subdir_for_output_files(pdb_cc.chem_comp_id + '.ideal_withH.sdf')
        pdb_cc.sdf_file_or_string(file_name=sdf_ideal_with_h)
        yield assert_true, os.path.isfile(sdf_ideal_with_h) and os.path.getsize(sdf_ideal_with_h) > 0, \
            '{} call to pdb_cc.sdf_file_or_string(file="{}") must create a non-empty file.'.\
            format(pdb_cc.chem_comp_id, sdf_ideal_with_h)
        sdf_model_with_h = file_name_in_subdir_for_output_files(pdb_cc.chem_comp_id + '.model_withH.sdf')
        pdb_cc.sdf_file_or_string(file_name=sdf_model_with_h,ideal=False)
        yield assert_true, os.path.isfile(sdf_model_with_h) and os.path.getsize(sdf_model_with_h) > 0, \
            '{} call to pdb_cc.sdf_file_or_string(file="{}") must create a non-empty file.'.\
            format(pdb_cc.chem_comp_id, sdf_model_with_h)


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
