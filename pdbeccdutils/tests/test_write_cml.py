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
from lxml import etree

from nose.tools import assert_equals, assert_true
from utilities import cif_filename, file_name_in_tsts_out, test_cif_path_name
from pdb_chemical_components_rdkit import PdbChemicalComponentsRDKit


def test_load_eoh_from_cif():
    eoh = PdbChemicalComponentsRDKit(file_name=cif_filename('EOH'))
    cml = file_name_in_tsts_out('EOH.cml')
    cml_string = eoh.cml_file_or_string()
    tree = etree.fromstring(cml_string)
    yield assert_equals, tree[0][1].text, eoh.chem_comp_name,\
        'cml_file_or_string must provide correct systematic name'
    cml_file = open (os.path.join(test_cif_path_name, 'EOH.cml'), 'r')
    cml_tree = etree.parse(cml_file)
    cml_root = cml_tree.getroot()
    formal_charge_cml = cml_root[0].attrib['formalCharge']
    yield assert_equals, tree[0].attrib['formalCharge'], formal_charge_cml,\
        'cml_file_or_string must provide correct formal charge'
    eoh.cml_file_or_string(file_name=cml)
    yield assert_true, os.path.isfile(cml) and os.path.getsize(cml) > 0, \
        '{} call to eoh.cml_file_or_string(file="{}") must create a non-empty file.'.\
        format('EOH', cml)

def test_load_atp_from_cif():
    atp = PdbChemicalComponentsRDKit(file_name=cif_filename('ATP'))
    cml = file_name_in_tsts_out('ATP.cml')
    cml_string = atp.cml_file_or_string()
    tree = etree.fromstring(cml_string)
    cml_file = open(os.path.join(test_cif_path_name, 'ATP.cml'), 'r')
    cml_tree = etree.parse(cml_file)
    cml_root = cml_tree.getroot()
    for item in cml_root[0].findall('formula'):
        if item.attrib['dictRef'] == 'ebiMolecule:stereoSmiles':
            cml_stereo = item.text
        if item.attrib['dictRef'] == 'ebiMolecule:nonStereoSmiles':
            cml_nonstereo = item.text
    yield assert_equals, tree[0][2].text, cml_stereo,\
        'cml_file_or_string must provide correct stereo smiles'
    yield assert_equals, tree[0][3].text, cml_nonstereo,\
        'cml_file_or_string must provide correct nonstereo smiles'

class DummyTestCaseSoPycharmRecognizesNoseTestsAsTests(unittest.TestCase):
    pass
