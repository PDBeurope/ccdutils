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
from rdkit.Geometry import rdGeometry
from rdkit.Chem.rdmolops import AssignAtomChiralTagsFromStructure


class PdbChemicalComponentsRDKit(PdbChemicalComponents):
    """ PdbChemicalComponents class with extension to produce RDKit Mol"""

    def __init__(self, file_name=None, cif_parser='auto'):
        super(PdbChemicalComponentsRDKit, self).__init__(file_name, cif_parser)
        self.rdkit_mol = None
        self.rdkit_model_mol = None
        self.setup_rdkit_mol()
        self._inchikey_from_rdkit = None

    def setup_rdkit_mol(self):
        # use method from http://rdkit-discuss.narkive.com/RVC3HZjy/building-mol-manually
        empty_mol = Chem.Mol()  # creates a blank molecule
        self.rdkit_mol = Chem.RWMol(empty_mol)
        self.rdkit_model_mol = Chem.RWMol(empty_mol)
        for atom_index in range(self.number_atoms):
            element = self.atom_elements[atom_index]
            atom_name = self.atom_ids[atom_index]
            rdkit_atom = Chem.Atom(element)
            rdkit_atom.SetProp('name', atom_name)  # from sameer_prototype_chem.py
            # set the name of the atom to be included in sdf file Alias lines
            # https://gist.github.com/ptosco/6e4468350f0fff183e4507ef24f092a1#file-pdb_atom_names-ipynb
            rdkit_atom.SetProp('molFileAlias', atom_name)
            charge = self.atom_charges[atom_index]
            rdkit_atom.SetFormalCharge(charge)
            self.rdkit_mol.AddAtom(rdkit_atom)
            self.rdkit_model_mol.AddAtom(rdkit_atom)
        for bond_index in range(self.number_bonds):
            index_1 = self.bond_atom_index_1[bond_index]
            index_2 = self.bond_atom_index_2[bond_index]
            order = Chem.rdchem.BondType(self.bond_order[bond_index])
            self.rdkit_mol.AddBond(index_1, index_2, order)
            self.rdkit_model_mol.AddBond(index_1, index_2, order)
        Chem.SanitizeMol(self.rdkit_mol, catchErrors=True)
        Chem.SanitizeMol(self.rdkit_model_mol, catchErrors=True)
        Chem.Kekulize(self.rdkit_mol)
        Chem.Kekulize(self.rdkit_model_mol)
        ideal_conformer = Chem.Conformer(self.number_atoms)
        model_conformer = Chem.Conformer(self.number_atoms)
        for atom_index in range(self.number_atoms):
            (ideal_x, ideal_y, ideal_z) = self.ideal_xyz[atom_index]
            (model_x, model_y, model_z) = self.model_xyz[atom_index]
            rdkit_xyz = rdGeometry.Point3D(ideal_x, ideal_y, ideal_z)
            rdkit_model_xyz = rdGeometry.Point3D(model_x, model_y, model_z)
            ideal_conformer.SetAtomPosition(atom_index, rdkit_xyz)
            model_conformer.SetAtomPosition(atom_index, rdkit_model_xyz)
        self.rdkit_mol.AddConformer(ideal_conformer)
        AssignAtomChiralTagsFromStructure(self.rdkit_mol)
        self.rdkit_model_mol.AddConformer(model_conformer)
        AssignAtomChiralTagsFromStructure(self.rdkit_model_mol)

    @property
    def inchikey_from_rdkit(self):
        if self._inchikey_from_rdkit is None:
            inchi = Chem.inchi.MolToInchi(self.rdkit_mol)
            self._inchikey_from_rdkit = Chem.inchi.InchiToInchiKey(inchi)
        return self._inchikey_from_rdkit

    def sdf_file_or_string(self, file_name=None, ideal=True, hydrogen=True, alias=True):
        """
        write a sdf file or return a string containing the molecule as a sdf file

        Args:
            file_name (str): optional filename
            ideal (bool): write the ideal coordinates (if True) or model coordinates (False)
            hydrogen (bool): include hydrogen atoms in the sdf
            alias (bool): use the alias feature to include atom names in the sdf

        Returns:
            None or a string containing the molecule converted to sdf

        Notes:
            TODO currently limited to writing the ideal coordinates with hydrogen atoms
            This method should not alter self.rdkit_mol by removing hydrogen atoms etc.
        """
        if not ideal:
            #raise NotImplementedError('sdf_file_or_string ideal=False to be coded')  # TODO implement ideal
            sdf_string = Chem.MolToMolBlock(self.rdkit_model_mol)
            if file_name is None:
                return sdf_string
            else:
                with open(file_name,'w') as sdf_file:
                    sdf_file.write(sdf_string)
        if not hydrogen:
            raise NotImplementedError('sdf_file_or_string hydrogen=False to be coded')  # TODO implement hydrogen
        if not alias:
            raise NotImplementedError('sdf_file_or_string alias=False to be coded')  # TODO implement alias
        if ideal and hydrogen and alias:   
            sdf_string = Chem.MolToMolBlock(self.rdkit_mol)
            if file_name is None:
                return sdf_string
            else:
                with open(file_name, "w") as sdf_file:
                    sdf_file.write(sdf_string)
            return None

    def pdb_file_or_string(self):
        # TODO implement pdb_file_or_string - most options like sdf_file_or_string
        raise NotImplementedError('to be coded')

    def cml_file_or_string(self):
        # TODO implement cml_file_or_string - not sure about options!
        raise NotImplementedError('to be coded')

    def image_file(self, file_name=None, wedge=False, atom_labels=False, hydrogen=False):
        """
        writes an image of the molecule to a file using rdkit.Chem.Draw

        Args:
            file_name (str): the name of the file. Type normally got from the filename ending for instance .png or .svg
            wedge (bool):  wedge the bonds in the image
            atom_labels (str): include atom labels in the image
            hydrogen (str): include hydrogen atoms in the image.

        Returns:

        """
        raise NotImplementedError('to be coded')
