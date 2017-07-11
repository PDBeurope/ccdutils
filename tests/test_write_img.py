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

from nose.tools import assert_in, assert_true

from pdb_chemical_components_rdkit import PdbChemicalComponentsRDKit
from utilities import cif_filename, test_file_path_name, file_name_in_tsts_out


def test_load_eoh_from_cif():
    eoh = PdbChemicalComponentsRDKit(file_name=cif_filename('EOH'))
    img_with_h_no_label = file_name_in_tsts_out('EOH.img_withH_no_label.svg')
    eoh.image_file_or_string(file_name=img_with_h_no_label, hydrogen=True, atom_labels=False)
    yield assert_true, os.path.isfile(img_with_h_no_label) and os.path.getsize(img_with_h_no_label) > 0, \
        '{} call to eoh.image_file_or_string(file="{}") must create a non-empty file.'.\
        format('EOH', img_with_h_no_label)
    img_with_h_label = file_name_in_tsts_out('EOH.img_withH_label.svg')
    eoh.image_file_or_string(file_name=img_with_h_label, hydrogen=True, atom_labels=True)
    yield assert_true, os.path.isfile(img_with_h_label) and os.path.getsize(img_with_h_label) > 0, \
        '{} call to eoh.image_file_or_string(file="{}") must create a non-empty file.'.\
        format('EOH', img_with_h_label)

def test_atp_svg():
    atp = PdbChemicalComponentsRDKit(file_name=cif_filename('ATP'))
    img_no_h_label_wdege = file_name_in_tsts_out('ATP.img_noH_no_label.svg')
    atp_svg = atp.image_file_or_string(hydrogen=False, atom_labels=True, wedge=True)
    yield assert_in, 'PA', atp_svg, 'The svg file must contain the atom name PA'
    

def test_sdf_write_for_all_sample_cifs():
    for ciffile in supply_list_of_sample_cifs():
        #if 'HEM' in ciffile:
        #   continue
        pdb_cc = PdbChemicalComponentsRDKit(file_name=ciffile)
        img_no_h_no_label = file_name_in_tsts_out(pdb_cc.chem_comp_id + '.img_noH_no_label.svg')
        pdb_cc.image_file_or_string(file_name=img_no_h_no_label, hydrogen=False, atom_labels=False)
        yield assert_true, os.path.isfile(img_no_h_no_label) and os.path.getsize(img_no_h_no_label) > 0, \
            '{} call to eoh.image_file_or_string(file="{}") must create a non-empty file.'.\
            format(pdb_cc.chem_comp_id, img_no_h_no_label)
        img_no_h_label_wedge = file_name_in_tsts_out(pdb_cc.chem_comp_id + '.img_label_wedge.svg')
        pdb_cc.image_file_or_string(file_name=img_no_h_label_wedge, hydrogen=False, atom_labels=True, wedge=True)
        yield assert_true, os.path.isfile(img_no_h_label_wedge) and os.path.getsize(img_no_h_label_wedge) > 0, \
            '{} call to eoh.image_file_or_string(file="{}") must create a non-empty file.'.\
            format(pdb_cc.chem_comp_id, img_no_h_label_wedge)

def supply_list_of_sample_cifs():
    """
    returns the list of sample pdb ccd cifs for test.

    Args:
        None

    Returns:
        list of filenames
    """
    return sorted(glob.glob(os.path.join(test_file_path_name, '*.cif')))


class DummyTestCaseSoPycharmRecognizesNoseTestsAsTests(unittest.TestCase):
    pass
