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

import logging
import pprint
import sys

from rdkit import Chem

from utilities import fragment_library_file_path


class FragmentLibrary(object):
    def __init__(self):
        self.smiles_to_fragment_name = {}
        """{} smiles to fragment name dict {str,str}"""
        self.smiles_to_rdkit_molecule = {}

    def load(self):
        """
        loads fragment library from the file in data.
        """
        self._load_smiles_to_fragment_name_from_file(fragment_library_file_path)
        self._create_smiles_to_rdkit_mol()

    def fragments_for_pdb_chemical_components_rdkit(self, pdb_ccd_rdkit):
        """

        Args:
            pdb_ccd_rdkit: A PdbChemicalComponentsRDKit molecule

        Returns:
            the fragments in the molecule ????
        """
        raise NotImplementedError('to be written')

    def _load_smiles_to_fragment_name_from_file(self, fragment_file_name):
        """
        loads the smiles to fragment name dictionary from the given file.

        Args:
            fragment_file_name (str): the fragment file name (normally smi.text in the data directory)

        Note:
            fragment_file_name has format
            cyclopropane:C1CC1
            phenyl:c1ccccc1
        """
        try:
            fragment_file = open(fragment_file_name, 'r')
        except IOError as err:
            print("Error cannot open fragment file {} error is '{}'".format(fragment_file_name, err.strerror))
            sys.exit(1)
        self.smiles_to_fragment_name = {}
        lines = fragment_file.read().splitlines()
        for line in lines:
            name, smile = line.split(':')
            smile = smile.replace('\t', '')  # take out tabs
            self.smiles_to_fragment_name[smile] = name
        number_of_entries = len(self.smiles_to_fragment_name)
        logging.debug('method load_smiles_from_smi_text_file:')
        logging.debug('Have loaded smiles_to_fragment_name dictionary with {} '
                      'entries from file {}'.format(number_of_entries, fragment_file_name))
        logging.debug('\tdump entries:\n' + pprint.pformat(self.smiles_to_fragment_name))

    def _create_smiles_to_rdkit_mol(self):
        """
        creates a dictionary with an rdkit molecule for each SMILES string in the input list.
        """
        self.smiles_to_rdkit_molecule = {}
        for smile in self.smiles_to_fragment_name.keys():
            rdkit_mol = Chem.MolFromSmiles(smile)
            self.smiles_to_rdkit_molecule[smile] = rdkit_mol
