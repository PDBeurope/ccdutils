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
from yattag import Doc

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

    def html_report_of_fragments(self, pdb_ccd_rdkit, fragments, html_file_name):
        """
        idea produce an html report to a file of the fragments identified within a molecule

        Args:
            pdb_ccd_rdkit: A PdbChemicalComponentsRDKit molecule
            fragments: the fragments in the molecule - a dictionary fragment name: list of (list of atom ids/names)
            html_file_name (str): file name for the output html file
        """

        html_text = self.prepare_html(pdb_ccd_rdkit, fragments)
        try:
            with open(html_file_name, "w") as text_file:
                text_file.write(html_text)
        except IOError as e_mess:
            print('Error cannot open or write html file - problem:\n{}'.format(e_mess))
            sys.exit(1)

    def prepare_html(self, pdb_ccd_rdkit, fragments):
        """
        prepares an html document of the fragments within a chem_comp
        
        Args:
            pdb_ccd_rdkit: A PdbChemicalComponentsRDKit molecule
            fragments: the fragments in the molecule - a dictionary fragment name: list of (list of atom ids/names)
        Returns:
            str: containing html document
        """
        fragments_ordered_by_number_atoms = OrderedDict(sorted(fragments.items(), key=lambda x: -len(x[1][0])))
        doc, tag, text, line = Doc().ttl()
        title = 'fragments for {}'.format(pdb_ccd_rdkit.chem_comp_id)
        with tag('html'):
            with tag('head'):
                with tag('title'):
                    text(title)
            with tag('body'):
                with tag('h1'):
                    text(title)
                with tag('p'):
                    text(' chem_comp_name: {}'.format(pdb_ccd_rdkit.chem_comp_name))
                for (name, list_of_atom_list) in fragments_ordered_by_number_atoms.items():
                    with tag('h3'):
                        text('{}'.format(name))
                    with tag('p'):
                        for atoms in list_of_atom_list:
                            svg_diagram = pdb_ccd_rdkit.image_file_or_string(highlight_atoms=atoms,
                                                                             atom_labels=False, black=True)
                            svg_diagram = svg_diagram.replace('svg:', '')
                            doc.asis(svg_diagram)
                    with tag('p'):
                        description = self.information_for_fragment_name(name, query_type='description')
                        if description is not None:
                            doc.asis(description + ' ')
                        url = self.information_for_fragment_name(name, query_type='url')
                        if url is not None:
                            doc.stag('br')
                            with tag('a', href=url):
                                text(url)

        result = doc.getvalue()
        return result

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

    def information_for_fragment_name(self, fragment_name, query_type='query'):
        """
        provides the information (default the query SMILES or SMARTS string) for a given fragment name

        
        Args:
            fragment_name (str): name of fragment 
            query_type (str): the tsv column header for example url
         
        Returns:
            str: the required information from the tsv file or None if fragment name not found. 
        """
        df = self.fragments_pd
        select = df.loc[df['name'] == fragment_name]
        if select.shape[0] == 0:
            return None
        else:
            value = select[query_type].values[0]
            if str(value) == 'nan':
                return None
            return value

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
