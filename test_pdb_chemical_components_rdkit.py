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
from nose.tools import assert_equals
from test_pdb_chemical_components import cif_filename
from pdb_chemical_components_rdkit import PdbChemicalComponentsRDKit


def test_hard_code_cmo():
    cmo = PdbChemicalComponentsRDKit(cif_parser='test_hard_code_cmo')
    cmo.load_carbon_monoxide_hard_coded()
    yield assert_equals, 'CMO', cmo.chem_comp_id, 'chem_comp_id'
    yield assert_equals, 'CARBON MONOXIDE', cmo.chem_comp_name, 'chem_comp_name'
    yield assert_equals, 'UGFAIRIUMAVXCW-UHFFFAOYSA-N', cmo.inchikey, 'chem_inchikey'
    # rdkit should be able to get an inchikey that is the same that read from the CCD cif file
    yield assert_equals, 2, cmo.rdkit_mol.GetNumAtoms(), 'rdkit_mol.GetNumAtoms() should give 2'
    yield assert_equals,  cmo.inchikey, cmo.rdkit_inchikey, 'inchikey from cif file should match the rdkit inchikey'


def test_load_eoh_from_cif():
    eoh = PdbChemicalComponentsRDKit(file_name=cif_filename('EOH'))
    yield assert_equals, 'EOH', eoh.chem_comp_id, 'chem_comp_id'
    yield assert_equals, eoh.inchikey, eoh.rdkit_inchikey, 'inchikey from cif file should match the rdkit inchikey'
