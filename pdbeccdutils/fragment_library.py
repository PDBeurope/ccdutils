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
from collections import OrderedDict
from rdkit import Chem
from pdbeccdutils.utilities import fragment_library_file_path, this_script_dir


class FragmentLibrary(object):
    """
    store common fragments from
    """
    def __init__(self):
        self.fragment_name_to_smiles = OrderedDict()
        self.fragment_name_to_rdkit_molecule = OrderedDict()
        self._load()

    @property
    def number_of_entries(self):
        return len(self.fragment_name_to_smiles)

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
        fragments = {}
        rdkit_mol = pdb_ccd_rdkit.mol_remove_h
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

    def _load(self):
        """
        loads fragment library from the file in data.
        """
        self._load_fragment_name_to_smiles_from_file(fragment_library_file_path)
        self._create_frag_name_to_rdkit_mol()

    def _load_fragment_name_to_smiles_from_file(self, fragment_file_name):
        """
        loads the fragment name to file dictionary from the given file.

        Args:
            fragment_file_name (str): the fragment file name 

        Note:
            fragment_file_name has SMILES format:
            C1CC1 cyclopropane
            c1ccccc1 phenyl
        """
        with open(fragment_file_name, 'r') as fragment_file:
            self.fragment_name_to_smiles = OrderedDict()
            lines = fragment_file.read().splitlines()
            for line in lines:
                smile, name = line.split(' ')
                smile = smile.replace('\t', '')  # take out tabs
                self.fragment_name_to_smiles[name] = smile

    def _create_frag_name_to_rdkit_mol(self):
        """
        creates a dictionary with an rdkit molecule for each fragment in the input list.
        """
        self.fragment_name_to_rdkit_molecule = OrderedDict()
        for fragment_name, smile in self.fragment_name_to_smiles.items():
            rdkit_mol = Chem.MolFromSmiles(smile)
            self.fragment_name_to_rdkit_molecule[fragment_name] = rdkit_mol


def produce_png_img_of_all_fragments():
    from rdkit.Chem import Draw
    frag_library = FragmentLibrary()
    ms = list(frag_library.fragment_name_to_rdkit_molecule.values())
    labels = list(frag_library.fragment_name_to_smiles.keys())
    img = Draw.MolsToGridImage(ms, legends=labels, molsPerRow=10)
    img_file_name = os.path.join(this_script_dir(), 'data', 'fragment_library.png')
    img.save(img_file_name)
    print('Have written image of all the fragments to {}'.format(img_file_name))


if __name__ == "__main__":
    produce_png_img_of_all_fragments()
