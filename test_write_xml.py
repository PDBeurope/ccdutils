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
    xml = file_name_in_subdir_for_output_files('EOH.xml')
    eoh.xml_file_or_string(file_name=xml)
    yield assert_true, os.path.isfile(xml) and os.path.getsize(xml) > 0, \
        '{} call to eoh.image_file(file="{}") must create a non-empty file.'.\
        format('EOH', xml)

#def test_sdf_write_for_all_sample_cifs():
#    for ciffile in supply_list_of_sample_cifs():
        #if 'HEM' in ciffile:
        #   continue
#        pdb_cc = PdbChemicalComponentsRDKit(file_name=ciffile)
#        img_no_h_no_label = file_name_in_subdir_for_output_files(pdb_cc.chem_comp_id + '.img_noH_no_label.svg')
#        pdb_cc.image_file(file_name=img_no_h_no_label, hydrogen=False, atom_labels=False)
#        yield assert_true, os.path.isfile(img_no_h_no_label) and os.path.getsize(img_no_h_no_label) > 0, \
#            '{} call to eoh.image_file(file="{}") must create a non-empty file.'.\
#            format(pdb_cc.chem_comp_id, img_no_h_no_label)
#        img_no_h_label_wedge = file_name_in_subdir_for_output_files(pdb_cc.chem_comp_id + '.img_label_wedge.svg')
#        pdb_cc.image_file(file_name=img_no_h_label_wedge, hydrogen=False, atom_labels=True, wedge=True)
#        yield assert_true, os.path.isfile(img_no_h_label_wedge) and os.path.getsize(img_no_h_label_wedge) > 0, \
#            '{} call to eoh.image_file(file="{}") must create a non-empty file.'.\
#            format(pdb_cc.chem_comp_id, img_no_h_label_wedge)
#
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
