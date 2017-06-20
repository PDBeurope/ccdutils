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
from rdkit.Chem import Draw, AllChem, rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D, IPythonConsole
from IPython.display import SVG
import image
from lxml import etree


class PdbChemicalComponentsRDKit(PdbChemicalComponents):
    """ PdbChemicalComponents class with extension to produce RDKit Mol"""

    def __init__(self, file_name=None, cif_parser='auto'):
        super(PdbChemicalComponentsRDKit, self).__init__(file_name, cif_parser)
        self.rdkit_mol = None
        self.rdkit_mol_conformer_id_ideal = None
        """The RDKit conformer ID for the ideal cooordinate (int)."""
        self.rdkit_mol_conformer_id_model = None
        """The RDKit conformer ID for the model cooordinate (int)."""
        self.__setup_rdkit_mol()
        self._inchikey_from_rdkit = None

    def __setup_rdkit_mol(self):
        """
        setup the rdkit mol called by the initializer only

        Returns:
            None

        """
        # use method from http://rdkit-discuss.narkive.com/RVC3HZjy/building-mol-manually
        empty_mol = Chem.Mol()
        self.rdkit_mol = Chem.RWMol(empty_mol)
        self.__setup_atoms()
        self.__setup_bonds()
        self.__setup_conformers()
        Chem.SanitizeMol(self.rdkit_mol, catchErrors=True)
        Chem.Kekulize(self.rdkit_mol)
        AssignAtomChiralTagsFromStructure(self.rdkit_mol)
        self.mol_remove_h = Chem.RWMol(self.rdkit_mol)
        try:
            self.mol_remove_h = Chem.RemoveHs(self.mol_remove_h)
        except ValueError:
            pass  # TODO deal with this error properly

    def __setup_atoms(self):
        """
        sets up the atoms in the rdkit mol - elements, atom names, sdf alias name, charges.

        Returns:
            None
        """
        for atom_index in range(self.number_atoms):
            element = self.atom_elements[atom_index]
            atom_name = self.atom_ids[atom_index]
            rdkit_atom = Chem.Atom(element)
            rdkit_atom.SetProp('name', atom_name)  # from sameer_prototype_chem.py
            charge = self.atom_charges[atom_index]
            rdkit_atom.SetFormalCharge(charge)
            self.rdkit_mol.AddAtom(rdkit_atom)

    @staticmethod
    def __sdf_alias_on_off(rwmol=None, alias=False):
        """
        turns on or off the labels used for atom alias records in sdf output

        Args:
            rwmol: The RDkit RWMol molecule to be written
            alias (bool: turn on (True) or off (False)

        Returns:
            None

        Notes:
            procedure used based on
            https://gist.github.com/ptosco/6e4468350f0fff183e4507ef24f092a1#file-pdb_atom_names-ipynb

        """
        if rwmol is None:
            raise ValueError('call to __sdf_alias with rwmol is None')
        for rdkit_atom in rwmol.GetAtoms():
            if alias:
                atom_name = rdkit_atom.GetProp('name')
                rdkit_atom.SetProp('molFileAlias', atom_name)
            else:
                rdkit_atom.ClearProp('molFileAlias')

    def __setup_bonds(self):
        """
        sets up the bonds in the rdkit mol - index numbers and bond orders

        Returns:
            None

        Notes:
            we do not use the CCD aromatic records pdbx_aromatic_flag
        """
        for bond_index in range(self.number_bonds):
            index_1 = self.bond_atom_index_1[bond_index]
            index_2 = self.bond_atom_index_2[bond_index]
            order = Chem.rdchem.BondType(self.bond_order[bond_index])
            self.rdkit_mol.AddBond(index_1, index_2, order)

    def __setup_conformers(self):
        """
        loads the ideal and model xyz coordinates as separate rdkit conformers.

        Returns:
            None
        """
        ideal_conformer = Chem.Conformer(self.number_atoms)
        model_conformer = Chem.Conformer(self.number_atoms)
        for atom_index in range(self.number_atoms):
            (ideal_x, ideal_y, ideal_z) = self.ideal_xyz[atom_index]
            (model_x, model_y, model_z) = self.model_xyz[atom_index]
            rdkit_xyz = rdGeometry.Point3D(ideal_x, ideal_y, ideal_z)
            rdkit_model_xyz = rdGeometry.Point3D(model_x, model_y, model_z)
            ideal_conformer.SetAtomPosition(atom_index, rdkit_xyz)
            model_conformer.SetAtomPosition(atom_index, rdkit_model_xyz)
        self.rdkit_mol_conformer_id_ideal = self.rdkit_mol.AddConformer(ideal_conformer, assignId=True)
        self.rdkit_mol_conformer_id_model = self.rdkit_mol.AddConformer(model_conformer, assignId=True)

    @property
    def inchikey_from_rdkit(self):
        if self._inchikey_from_rdkit is None:
            inchi = Chem.inchi.MolToInchi(self.rdkit_mol)
            self._inchikey_from_rdkit = Chem.inchi.InchiToInchiKey(inchi)
        return self._inchikey_from_rdkit

    def sdf_file_or_string(self, file_name=None, ideal=True, hydrogen=True, alias=False):
        """
        write a sdf file or return a string containing the molecule as a sdf file

        Args:
            file_name (str): optional filename
            ideal (bool): write the ideal coordinates (True) or model coordinates (False)? Default True: ideal.
            hydrogen (bool): include hydrogen atoms in the sdf? Default True: yes)
            alias (bool): use the alias feature to include atom names in the sdf? Default False no.

        Returns:
            None or a string containing the molecule converted to sdf
        """
        fname = self.chem_comp_id + '.sdf'
        if ideal:
            conformer_id = self.rdkit_mol_conformer_id_ideal
        else:
            conformer_id = self.rdkit_mol_conformer_id_model
        if hydrogen:
            mol_h_select = self.rdkit_mol
        else:
            mol_h_select = self.mol_remove_h

        self.__sdf_alias_on_off(mol_h_select, alias=alias)

        sdf_string = Chem.MolToMolBlock(mol_h_select, confId=conformer_id)
        sdf_string = fname + sdf_string
        if file_name is None:
            return sdf_string
        else:
            with open(file_name, "w") as sdf_file:
                sdf_file.write(sdf_string)
        return None

    def pdb_file_or_string(self, file_name=None, ideal=True):
        """
        write a pdb file or return a string containing the molecule as a pdb file

        Args:
            file_name (str): optional filename
            ideal (bool): write the ideal coordinates (True) or model coordinates (False)? Default True: ideal.

        Returns:
            None or a string containing the molecule converted to pdb
        """
        atom_pdb_residue_info = Chem.rdchem.AtomPDBResidueInfo()
        atom_pdb_residue_info.SetResidueName(self.chem_comp_id)
        atom_pdb_residue_info.SetTempFactor(20.0)
        atom_pdb_residue_info.SetOccupancy(1.0)
        atom_pdb_residue_info.SetChainId('A')
        atom_pdb_residue_info.SetResidueNumber(1)
        for atom in self.rdkit_mol.GetAtoms():
            atom_index = atom.GetIdx()
            pdbx_align = self.atom_pdbx_align[atom_index]
            element = self.atom_elements[atom_index]
            atom_name = self.atom_ids[atom_index]
            if len(atom_name) < 4:
                atom_name = atom_name + ' ' * (3 - len(atom_name) + (pdbx_align == '0'))
            if pdbx_align == '1':
                atom_name = ' ' + atom_name
            atom_pdb_residue_info.SetName(atom_name)
            atom.SetMonomerInfo(atom_pdb_residue_info)
        if ideal:
            conformer_id = self.rdkit_mol_conformer_id_ideal
        else:
            conformer_id = self.rdkit_mol_conformer_id_model
        pdb_title1 = 'HEADER    NONAME 07-Jun-17\nTITLE     Produced by PDBeChem\nCOMPND    {}\n'.format(
            self.chem_comp_id)
        pdb_title2 = 'AUTHOR    EBI-PDBe Generated\nREVDAT   1  07-Jun-17     0\n'
        pdb_string = pdb_title1 + pdb_title2 + Chem.MolToPDBBlock(self.rdkit_mol, conformer_id)
        if file_name is None:
            return pdb_string
        else:
            with open(file_name, "w") as pdb_file:
                pdb_file.write(pdb_string)
        return None

    def xml_file_or_string(self, file_name=None):
        # TODO implement xml_file_or_string - not sure about options!
        #raise NotImplementedError('to be coded')
        top = etree.Element('xml')
        top.set('dictRef','ebiMolecule:ebiMoleculeDict.xml')
        top.set('ebiMolecule','http://www.ebi.ac.uk/felics/molecule')
        mol = etree.SubElement(top, 'molecule', id=self.chem_comp_id, formalcharge='0')#Need charge
        id_inchi = etree.SubElement(mol,'identifier', dictRef='ebiMolecule:inchi')
        id_inchi.text = self.inchikey
        id_systematic = etree.SubElement(mol,'identifier', dictRef='ebiMolecule:systematicName')
        id_systematic.text = self.chem_comp_name
        id_formula1 = etree.SubElement(mol, 'identifier',dictRef="ebiMolecule:stereoSmiles")
        #TODO add methods for getting smile file
        atom_array = etree.SubElement(mol, 'atomArray')
        for atom_index in range(self.number_atoms):
            element = self.atom_elements[atom_index]
            atom_name = self.atom_ids[atom_index]
            atom_entry = etree.SubElement(atom_array, 'atom', id=atom_name, elementType=element)
        bond_array = etree.SubElement(mol, 'bondArray')
        for bond_index in range(self.number_bonds):
            atom_1 = self.bond_atom_name_1[bond_index]
            atom_2 = self.bond_atom_name_2[bond_index]
            bond_order = Chem.rdchem.BondType(self.bond_order[bond_index])
            #print bond_order
            bond_entry = etree.SubElement(bond_array, 'bond')
            bond_entry.set('atomsRefs2',atom_1+' '+atom_2)
            bond_entry.set('order', str(bond_order))
        tree = etree.ElementTree(top)
        xml_string = etree.tostring(top, pretty_print=True)
        if file_name is None:
            return xml_string
        else:
            with open (file_name, 'w') as xml_file:
                xml_file.write(xml_string)
                xml_file.close()
        return None

    def image_file_or_string(self, file_name=None, wedge=True, atom_labels=True, hydrogen=False):
        """
        writes an image of the molecule to a file using rdkit.Chem.Draw

        Args:
            file_name (str): the name of the file. Type normally got from the filename ending for instance .png or .svg
            wedge (bool):  wedge the bonds in the image
            atom_labels (str): include atom labels in the image
            hydrogen (str): include hydrogen atoms in the image.

        Returns:
            None or a string containing the svg string of the molecule
        """
        if hydrogen:
            mol_h_select = self.rdkit_mol
        else:
            mol_h_select = self.mol_remove_h
        AllChem.GenerateDepictionMatching3DStructure(mol_h_select, mol_h_select)
        drawer = rdMolDraw2D.MolDraw2DSVG(400,200)
        opts = drawer.drawOptions()
        if not atom_labels:
            molecule_to_draw = rdMolDraw2D.PrepareMolForDrawing(mol_h_select, wedgeBonds=wedge)
        else:
            for atom in self.rdkit_mol.GetAtoms():
                atom_index = atom.GetIdx()
                atom_name = self.atom_ids[atom_index]
                opts.atomLabels[atom_index] = atom_name
            molecule_to_draw = rdMolDraw2D.PrepareMolForDrawing(mol_h_select)
        drawer.DrawMolecule(molecule_to_draw)
        drawer.FinishDrawing()
        svg = drawer.GetDrawingText().replace('svg:','')
        if file_name is None:
            return svg
        else:
            with open(file_name, 'w') as img_file:
                img_file.write(svg)
                img_file.close()
        return None
