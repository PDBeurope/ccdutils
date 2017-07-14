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
import collections
import os
import unittest

from nose.tools import assert_in, assert_true

from pdb_chemical_components_rdkit import PdbChemicalComponentsRDKit
from utilities import cif_filename, supply_list_of_sample_cifs, file_name_in_tsts_out


def test_eoh_svg():
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
    

def test_svg_write_for_all_sample_cifs():
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

def test_svg_highlight_bonds():
    mol = PdbChemicalComponentsRDKit(file_name=cif_filename('GOL'))
    this_img = file_name_in_tsts_out('test_svg_highlight_bonds.svg')
    # highlight first four bonds
    highlight_bonds = collections.OrderedDict()
    highlight_bonds[(0,1)] = (69./255., 117./255., 180./255.)  # C1-O1 blue
    highlight_bonds[(0,2)] = (145./255., 191./255., 219./255.)  # C1-C2 mid blue
    highlight_bonds[(2,3)] = (254./255., 224./255., 144./255.)  # C2-O2 light orange
    highlight_bonds[(2,4)] = (252./255., 141./255., 89./255.)   # C2-C3 mid orange
    highlight_bonds[(4,5)] = (215./255., 48./255., 39./255.)   # C3-O3 blood orange

    mol.image_file_or_string(file_name=this_img, hydrogen=False, atom_labels=False, wedge=False,
                             highlight_bonds=highlight_bonds, black=True)

class DummyTestCaseSoPycharmRecognizesNoseTestsAsTests(unittest.TestCase):
    pass
