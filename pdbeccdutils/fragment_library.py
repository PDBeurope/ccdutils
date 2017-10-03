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
import logging
import sys
from collections import OrderedDict

import pandas as pd
from rdkit import Chem

from pdbeccdutils.utilities import default_fragment_library_file_path


class FragmentLibrary(object):
    """
    store common named fragments and provide substructure search for them
    """
    def __init__(self, override_fragment_library_file_path=None):
        self.fragments_pd = None
        """panda data frame containing the information"""
        self.fragment_name_to_rdkit_molecule = OrderedDict()
        if override_fragment_library_file_path is None:
            fragment_file_name = default_fragment_library_file_path
        else:
            fragment_file_name = override_fragment_library_file_path
        self._load(fragment_file_name)

    @property
    def number_of_entries(self):
        return len(self.fragment_name_to_rdkit_molecule)

    def fragments_for_pdb_chemical_components_rdkit(self, pdb_ccd_rdkit):
        """
        Args:
            pdb_ccd_rdkit: A PdbChemicalComponentsRDKit molecule

        Returns:
            the fragments in the molecule - a dictionary fragment name: list of (list of atom ids/names)

        Notes
            There can be more than one match to a fragment in the same molecule - for instance cyclodextrin has 7
            pyranose rings. Each one will have a separate list of atom ids (names) that match the fragment.
            ence there is a list of these separate lists.
        """
        if pdb_ccd_rdkit.rwmol_cleaned_remove_h is not None:
            rdkit_mol = pdb_ccd_rdkit.rdkit_mol_with_any_metals_disconnected_from_n
        elif pdb_ccd_rdkit.rwmol_original_remove_h is not None:
            rdkit_mol = pdb_ccd_rdkit.rwmol_original_remove_h
        else:
            return {}
        fragments = OrderedDict()
        for fragment_name, frag_rdkit_mol in self.fragment_name_to_rdkit_molecule.items():
            if rdkit_mol.HasSubstructMatch(frag_rdkit_mol):
                fragments[fragment_name] = []
                # logging.debug('match to "{}"'.format(fragment_name))
                for all_index in rdkit_mol.GetSubstructMatches(frag_rdkit_mol):
                    # logging.debug('all_index={}'.format(all_index))
                    atom_names = []
                    for atom_index in all_index:
                        atom_names.append(pdb_ccd_rdkit.atom_ids[atom_index])
                    fragments[fragment_name].append(atom_names)
        return fragments

    def html_report_of_fragments(self, pdb_ccd_rdkit, html_file_name):
        """
        idea produce an html report to a file of the fragments identified within a molecule

        Args:
            pdb_ccd_rdkit: A PdbChemicalComponentsRDKit molecule
            html_file_name (str): file bane for the output html file
        """
        raise NotImplementedError('html_report_of_fragments to be written')


    def _load(self, fragment_file_name):
        """
        loads fragment library from the file in data.

        Args:
            fragment_file_name (str): the fragment file name
        """
        self._load_tsv_into_fragments_pd(fragment_file_name)
        self._create_frag_name_to_rdkit_mol()

    def _load_tsv_into_fragments_pd(self, fragment_file_name):
        logging.debug('parse {} into panda dataframe:'.format(fragment_file_name))
        self.fragments_pd = pd.read_csv(fragment_file_name, sep='\t')
        logging.debug('dataframe:\n{}'.format(self.fragments_pd.to_string()))
        columns = self.fragments_pd.columns.values.tolist()
        logging.debug('columns: {}'.format(columns))
        for required in 'name', 'kind', 'query':
            if not required in columns:
                raise ValueError('fragment library file lacks required column "{}"'.format(required))

    def _create_frag_name_to_rdkit_mol(self):
        """
        creates a dictionary with an rdkit molecule for each fragment in the input list.
        """
        self.fragment_name_to_rdkit_molecule = OrderedDict()
        for row in self.fragments_pd.itertuples():
            fragment_name = row.name
            query = row.query
            kind = row.kind
            if kind == 'SMARTS':
                rdkit_mol = Chem.MolFromSmarts(query)
            else:
                rdkit_mol = Chem.MolFromSmiles(query)
                if kind == 'LIKE':
                    for bond in rdkit_mol.GetBonds():
                        bond.SetBondType(Chem.rdchem.BondType.UNSPECIFIED)
            self.fragment_name_to_rdkit_molecule[fragment_name] = rdkit_mol

    def smiles_for_fragment_name(self, query_name):
        """
        provides the query SMILES or SMARTS string for a given fragment name
         
        Returns:
            str: the SMILES or SMARTS string or None if fragment name not found. 
        """
        df = self.fragments_pd
        select = df.loc[df['name'] == query_name]
        if select.shape[0] == 0:
            return None
        else:
            return select['query'].values[0]



def produce_png_img_of_all_fragments():
    from rdkit.Chem import Draw
    frag_library = FragmentLibrary()
    ms = list(frag_library.fragment_name_to_rdkit_molecule.values())
    labels = list(frag_library.fragment_name_to_rdkit_molecule.keys())
    img = Draw.MolsToGridImage(ms, legends=labels, molsPerRow=10)
    img_file_name = 'fragment_library.png'
    try:
        img.save(img_file_name)
        print('Have written image of all the fragments to {}'.format(img_file_name))
    except IOError as e_mess:
        print('Error cannot write image file to current directory - problem:\n{}'.format(e_mess))
        sys.exit(1)


if __name__ == "__main__":
    # logging.basicConfig(level=logging.DEBUG, format='%(levelname)s: %(message)s', )
    produce_png_img_of_all_fragments()
