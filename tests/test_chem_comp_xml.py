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
from nose.tools import assert_true, assert_false, assert_in, assert_not_in

from chem_comp_xml import ChemCompXMl
from pdb_chemical_components_rdkit import PdbChemicalComponentsRDKit
from utilities import cif_filename, file_name_in_tsts_out


def test_chem_comp_for_eoh_and_glu():
    cc_xml = ChemCompXMl()
    for chem_comp_id in 'EOH', 'GLU':
        ccd = PdbChemicalComponentsRDKit(file_name=cif_filename(chem_comp_id))
        cc_xml.store_ccd(ccd)
    cc_xml_string = cc_xml.to_string()
    old_records = ['<chemComp>', '<id>EOH</id>', '<name>ETHANOL</name>', '<formula>C2 H6 O</formula>',
                   '<systematicName>ethanol</systematicName>', '<stereoSmiles>CCO</stereoSmiles>',
                   '<nonStereoSmiles>CCO</nonStereoSmiles>', '<InChi>InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3</InChi>',
                   '</chemComp>',
                   # now want inchikey - not in old
                   '<InChIKey>LFQSCWFLJHTTHZ-UHFFFAOYSA-N</InChIKey>',
                   '<chemComp>', '<id>EOH</id>', '<name>ETHANOL</name>', '<formula>C2 H6 O</formula>',
                   '<systematicName>ethanol</systematicName>', '<stereoSmiles>CCO</stereoSmiles>',
                   '<nonStereoSmiles>CCO</nonStereoSmiles>', '<InChi>InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3</InChi>',
                   '</chemComp>',
                   # now want inchikey - not in old
                   '<InChIKey>LFQSCWFLJHTTHZ-UHFFFAOYSA-N</InChIKey>'
                   ]
    lines = cc_xml_string.splitlines()
    any_fail = False
    for old_record in old_records:
        match = False
        for line in lines:
            if old_record in line:
                match = True
                break
        yield assert_true, match, 'cc_xml_string should contain "{}"'.format(old_record)
        if not match:
            any_fail = True
    yield assert_false, any_fail, 'echo cc_xml_string after any failure: "{}"'.format(cc_xml_string)


def test_chem_comp_for_sy9_should_not_have_systematic_name():
    cc_xml = ChemCompXMl()
    ccd = PdbChemicalComponentsRDKit(file_name=cif_filename('SY9'))
    cc_xml.store_ccd(ccd)
    cc_xml_string = cc_xml.to_string()
    yield assert_in, '<id>SY9</id>', cc_xml_string
    yield assert_not_in, 'systematicName', cc_xml_string


def test_to_file():
    cc_xml = ChemCompXMl()
    ccd = PdbChemicalComponentsRDKit(file_name=cif_filename('SY9'))
    cc_xml.store_ccd(ccd)
    file_name = file_name_in_tsts_out('test_chem_comp_xml.xml')
    cc_xml.to_file(file_name)
    yield assert_true, os.path.isfile(file_name) and os.path.getsize(file_name) > 0, \
        'call to cc_xml.to_file({}) must create a non-empty file.'.format(file_name)
