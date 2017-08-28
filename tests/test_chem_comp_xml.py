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

from nose.tools import assert_true, assert_false

from pdb_chemical_components_rdkit import PdbChemicalComponentsRDKit
from utilities import cif_filename


def test_chem_comp_for_eoh():
    # record from old /nfs/ftp/pub/databases/msd/pdbechem/chem.xml
    #   <chemComp>
    #     <id>EOH</id>
    #     <name>ETHANOL</name>
    #     <formula>C2 H6 O</formula>
    #     <systematicName>ethanol</systematicName>
    #     <stereoSmiles>CCO</stereoSmiles>
    #     <nonStereoSmiles>CCO</nonStereoSmiles>
    #     <InChi>InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3</InChi>
    #   </chemComp>
    old_eoh_records = ('<chemComp>', '<id>EOH</id>', '<name>ETHANOL</name>', '<formula>C2 H6 O</formula>',
                       '<systematicName>ethanol</systematicName>', '<stereoSmiles>CCO</stereoSmiles>',
                       '<nonStereoSmiles>CCO</nonStereoSmiles>', '<InChi>InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3</InChi>',
                       '</chemComp>')
    eoh = PdbChemicalComponentsRDKit(file_name=cif_filename('EOH'))
    chem_comp_xml = eoh.chem_comp_xml()
    lines = chem_comp_xml.splitlines()
    any_fail = False
    for old_record in old_eoh_records:
        match = False
        for line in lines:
            if old_record in line:
                match = True
                break
        yield assert_true, match, 'EOH chem_comp_xml should contain "{}"'.format(old_record)
        if not match:
            any_fail = True
    yield assert_false, any_fail, 'echo EOH chem_comp_xml after any failure: "{}"'.format(chem_comp_xml)
