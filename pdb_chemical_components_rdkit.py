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
from pdb_chemical_components import PdbChemicalComponents
from rdkit import Chem


class PdbChemicalComponentsRDKit(PdbChemicalComponents):
    """ PdbChemicalComponents class with extension to produce RDKit Mol"""
    def __init__(self, file_name=None, cif_parser='auto'):
        super(PdbChemicalComponentsRDKit, self).__init__(file_name, cif_parser)
        self.rdkit_mol = None
        self.setup_rdkit_mol()
        self._inchikey_from_rdkit = None

    def setup_rdkit_mol(self):
        # use method from http://rdkit-discuss.narkive.com/RVC3HZjy/building-mol-manually
        empty_mol = Chem.Mol()  # creates a blank molecule
        self.rdkit_mol = Chem.RWMol(empty_mol)
        for atom_index in range(self.number_atoms):
            element = self.atom_elements[atom_index]
            name = self.atoms[atom_index]['atom_id']
            rdkit_atom = Chem.Atom(element)
            rdkit_atom.SetProp('name', name)  # from sameer_prototype_chem.py
            # set the name of the atom to be included in sdf file Alias lines
            # https://gist.github.com/ptosco/6e4468350f0fff183e4507ef24f092a1#file-pdb_atom_names-ipynb
            rdkit_atom.SetProp('molFileAlias', name)
            self.rdkit_mol.AddAtom(rdkit_atom)
        for bond_index in range(self.number_bonds):
            index_1 = self.bond_atom_index_1[bond_index]
            index_2 = self.bond_atom_index_2[bond_index]
            order = Chem.rdchem.BondType(self.bond_order[bond_index])
            self.rdkit_mol.AddBond(index_1, index_2, order)

    @property
    def inchikey_from_rdkit(self):
        if self._inchikey_from_rdkit is None:
            inchi = Chem.inchi.MolToInchi(self.rdkit_mol)
            self._inchikey_from_rdkit = Chem.inchi.InchiToInchiKey(inchi)
        return self._inchikey_from_rdkit
