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
from pdbeccdutils.pdb_chemical_components import PdbChemicalComponents
from pdbeccdutils.ring2dtemplates import supply_ring_2d_templates
# noinspection PyPackageRequirements
from rdkit import Chem
# noinspection PyPackageRequirements
from rdkit.Geometry import rdGeometry
# noinspection PyPackageRequirements
from rdkit.Chem.rdmolops import AssignAtomChiralTagsFromStructure
# noinspection PyPackageRequirements
from rdkit.Chem import AllChem
# noinspection PyPackageRequirements
from rdkit.Chem.Draw import rdMolDraw2D
from lxml import etree


class PdbChemicalComponentsRDKit(PdbChemicalComponents):
    """ PdbChemicalComponents class with extension to produce RDKit Mol"""

    def __init__(self, file_name=None, cif_dictionary=None, cif_parser='auto'):
        super(PdbChemicalComponentsRDKit, self).__init__(file_name, cif_dictionary, cif_parser)
        logging.debug('RDKit setup for PDB CCD chem_comp_id {}'.format(self.chem_comp_id))
        self.rwmol_original = None
        """The 'original' rdkit RWMol with bonds from PDB CCD cif not sanitized"""
        self.rwmol_original_remove_h = None
        """rwmol_original with hydrogen atoms removed - not sanitized"""
        self.rwmol_cleaned = None
        """Cleaned up rdkit RWMol with problematic bonds removed then sanitized """
        self.rwmol_cleaned_remove_h = None
        """rwmol_cleaned with hydrogen atoms removed, sanitized"""
        self.conformer_id_ideal = None
        """The RDKit conformer ID for the ideal cooordinate (int)."""
        self.conformer_id_model = None
        """The RDKit conformer ID for the model cooordinate (int)."""
        self.__setup_rdkit_mol()
        self.__create_rdkit_mol_cleaned()
        self.__xyz_2d = None
        self.__xyz_2d_no_hydrogen = None
        self._inchi_from_rdkit = None
        self._inchikey_from_rdkit = None

    def __setup_rdkit_mol(self):
        """
        setup the rdkit mol called by the initializer only

        Returns:
            None

        """
        # use method from http://rdkit-discuss.narkive.com/RVC3HZjy/building-mol-manually
        # noinspection PyArgumentList
        empty_mol = Chem.Mol()
        self.rwmol_original = Chem.RWMol(empty_mol)
        logging.debug('_setup_atoms:')
        self.__setup_atoms()
        logging.debug('_setup_bonds:')
        self.__setup_bonds()
        logging.debug('_setup_conformers:')
        self.__setup_conformers()
        logging.debug('removeH from original:')
        self.rwmol_original_remove_h = Chem.RemoveHs(self.rwmol_original, sanitize=False)

    def __setup_atoms(self):
        """
        sets up the atoms in the rdkit mol - elements, atom names, sdf alias name, charges.

        Returns:
            None
        """
        for atom_index in range(self.number_atoms):
            element = self.atom_elements[atom_index]
            isotope = None
            if element == 'D':
                element = 'H'
                isotope = 2
            elif element == 'X':
                element = 'O'
                isotope = 14
            atom_name = self.atom_ids[atom_index]
            rdkit_atom = Chem.Atom(element)
            rdkit_atom.SetProp('name', atom_name)  # from sameer_prototype_chem.py
            charge = self.atom_charges[atom_index]
            if isotope is not None:
                rdkit_atom.SetIsotope(isotope)
            if charge is None:
                charge = 0
            rdkit_atom.SetFormalCharge(charge)
            self.rwmol_original.AddAtom(rdkit_atom)

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
            self.rwmol_original.AddBond(index_1, index_2, order)

    def __create_rdkit_mol_cleaned(self):
        rwmol = Chem.RWMol(self.rwmol_original)
        # take metal to N bonds where N has valence 4
        smarts_metal_n_bond = Chem.MolFromSmarts('[Fe]~[N]')
        metal_n_bonds = rwmol.GetSubstructMatches(smarts_metal_n_bond)
        Chem.SanitizeMol(rwmol, sanitizeOps=Chem.SanitizeFlags.SANITIZE_CLEANUP)
        for (metal_index, n_index) in metal_n_bonds:
            metal_atom = rwmol.GetAtomWithIdx(metal_index)
            n_atom = rwmol.GetAtomWithIdx(n_index)
            n_valence = n_atom.GetExplicitValence()
            if n_valence == 4:
                # logging.debug('Remove bond between metal {} and nitrogen {} that has valence 4'
                #               .format(metal_atom.GetProp('name'), n_atom.GetProp('name')))
                # rwmol.RemoveBond(metal_index, n_index)
                logging.debug('Adjust charges for bond between metal {} and nitrogen {} that has valence 4'
                              .format(metal_atom.GetProp('name'), n_atom.GetProp('name')))
                n_atom.SetFormalCharge(1)
                metal_atom.SetFormalCharge(metal_atom.GetFormalCharge()-1)
        logging.debug('SanitizeMol rdkit_mol_cleaned')
        success = Chem.SanitizeMol(rwmol, catchErrors=True)
        if success != 0:
            logging.error('Sanitization failed with code {}'.format(success))
            return
        logging.debug('Kekulize rdkit_mol_cleaned')
        Chem.Kekulize(rwmol)
        AssignAtomChiralTagsFromStructure(rwmol)
        self.rwmol_cleaned = rwmol
        try:
            self.rwmol_cleaned_remove_h = Chem.RemoveHs(rwmol, sanitize=True)
        except ValueError as e_mess:
            logging.error('Problem in cleaned removeHs for {} with sanitize=True: {}'.
                          format(self.chem_comp_id, e_mess))
            self.rwmol_cleaned_remove_h = None

    def __setup_conformers(self):
        """
        loads the ideal and model xyz coordinates as separate rdkit conformers.

        Returns:
            None
        """
        ideal_conformer = Chem.Conformer(self.number_atoms)
        model_conformer = Chem.Conformer(self.number_atoms)
        for atom_index in range(self.number_atoms):
            (ideal_x, ideal_y, ideal_z) = self.replace_none_with_0_0_0(self.ideal_xyz[atom_index])
            (model_x, model_y, model_z) = self.replace_none_with_0_0_0(self.model_xyz[atom_index])
            rdkit_xyz = rdGeometry.Point3D(ideal_x, ideal_y, ideal_z)
            rdkit_model_xyz = rdGeometry.Point3D(model_x, model_y, model_z)
            ideal_conformer.SetAtomPosition(atom_index, rdkit_xyz)
            model_conformer.SetAtomPosition(atom_index, rdkit_model_xyz)
        self.conformer_id_ideal = self.rwmol_original.AddConformer(ideal_conformer, assignId=True)
        self.conformer_id_model = self.rwmol_original.AddConformer(model_conformer, assignId=True)

    def load_xyz_into_conformer(self, their_xyz):
        """
        loads a set of 3D coordinations as a rdkit conformer
        
        Args:
            their_xyz: the xyz coordinates of (x, y, z) one for each atom

        Returns:
            (int) the conformer id 
        """
        """  """
        new_conformer = Chem.Conformer(self.number_atoms)
        for atom_index in range(self.number_atoms):
            (x, y, z) = their_xyz[atom_index]
            rdkit_xyz = rdGeometry.Point3D(x, y, z)
            new_conformer.SetAtomPosition(atom_index, rdkit_xyz)
        conformer_id = self.rwmol_original.AddConformer(new_conformer, assignId=True)
        return conformer_id

    def __load_2d_coordinates_in_rwmol(self, rwmol):
        """
        loads 2D coordinates into a rwmol conformer

        Args:
            rwmol: RDKit rwmol object

        Returns:
            conformer id if successful or None on an error.
        """
        this_n_atoms = rwmol.GetNumAtoms()
        new_conformer = Chem.Conformer(this_n_atoms)
        if this_n_atoms == self.number_atoms:  # with hydrogens
            xyz = self.xyz_2d
        else:
            xyz = self.xyz_2d_no_hydrogen
        if xyz is None:
            return None
        else:
            for atom_index in range(this_n_atoms):
                (x, y, z) = xyz[atom_index]
                rdkit_xyz = rdGeometry.Point3D(x, y, z)
                new_conformer.SetAtomPosition(atom_index, rdkit_xyz)
            conformer_id = rwmol.AddConformer(new_conformer, assignId=True)
            return conformer_id

    @property
    def xyz_2d(self):
        """
        A set of 2D coordinates for pictures hydrogen atoms are included.

        Returns:
            tuple of tuple( x, y, z) for each atom. x, y, z are floats,
            z will be zero.
        """

        if self.__xyz_2d is None:
            self.__xyz_2d = self.__create_xyz_2d(hydrogen=True)
        return self.__xyz_2d

    @property
    def xyz_2d_no_hydrogen(self):
        """
        A set of 2D coordinates for pictures hydrogen atoms are not included.

        Returns:
            tuple of tuple( x, y, z) for each atom. x, y, z are floats,
            z will be zero.

        Notes:
            list will only be for non-hydrogen atoms.
        """

        if self.__xyz_2d_no_hydrogen is None:
            self.__xyz_2d_no_hydrogen = self.__create_xyz_2d(hydrogen=False)
        return self.__xyz_2d_no_hydrogen

    def __create_xyz_2d(self, hydrogen=False):
        """
        uses rdkit to create a set of 2D coords for the molecule.
        
        Args:
            hydrogen (bool):  include hydrogen atoms

        Returns:
            tuple of tuple( x, y, z) for each atom. x, y, z are floats, z will be zero.
            
        Notes:
            if hydrogen is False then returned list will only have coordinates for the non-hydrogen atoms.
        """
        this_xyz = []
        logging.debug('call to __create_xyz_2d with hydrogen={}'.format(hydrogen))
        if self.rwmol_cleaned is None:
            logging.warning('using original rwmol for {} image generation as clean up/sanitize failed'.
                            format(self.chem_comp_id))
            if hydrogen:
                mol_to_draw = self.rwmol_original
            else:
                mol_to_draw = self.rwmol_original_remove_h
        else:
            if hydrogen:
                mol_to_draw = self.rwmol_cleaned
            else:
                mol_to_draw = self.rwmol_cleaned_remove_h
        # make a copy as GenerateDepictionMatching3DStructure wipes existing conformations!
        mol_to_draw = Chem.RWMol(mol_to_draw)
        for template_name, template_mol in supply_ring_2d_templates().items():
            if mol_to_draw.HasSubstructMatch(template_mol):
                AllChem.GenerateDepictionMatching2DStructure(mol_to_draw, template_mol)
                logging.debug('have generated 2D coords from template "{}"'.format(template_name))
                break
        # noinspection PyArgumentList
        if mol_to_draw.GetNumConformers() != 1:
            try:
                AllChem.GenerateDepictionMatching3DStructure(mol_to_draw, mol_to_draw)
            except Exception as e_mess:
                logging.error('Problem for {} in generating 2D coords: {}'.
                              format(self.chem_comp_id, e_mess))
                return
        # noinspection PyArgumentList
        n_conformers = mol_to_draw.GetNumConformers()
        assert n_conformers == 1, 'should have one conformer'
        conformer = mol_to_draw.GetConformer(0)
        assert not conformer.Is3D(), 'conformer must be 2D'
        for position in conformer.GetPositions():
            float_x = float(position[0])
            float_y = float(position[1])
            float_z = float(position[2])
            this_xyz.append((float_x, float_y, float_z))
        this_xyz = tuple(this_xyz)
        return this_xyz

    @property
    def inchi_from_rdkit(self):
        """
        provides the InChI

        Returns:
            str: the InChI or 'ERROR' if there was an error finding it.
        """
        if self._inchi_from_rdkit is None:
            try:
                self._inchi_from_rdkit = Chem.inchi.MolToInchi(self.rwmol_original)
            except ValueError:
                self._inchi_from_rdkit = 'ERROR'
        return self._inchi_from_rdkit

    @property
    def inchikey_from_rdkit(self):
        """
        provides the InChIKey

        Returns:
            str: the InChIKey or 'ERROR' if there was an error finding it.
        """
        if self._inchikey_from_rdkit is None:
            inchi = self.inchi_from_rdkit
            if inchi != 'ERROR':
                self._inchikey_from_rdkit = Chem.inchi.InchiToInchiKey(inchi)
            else:
                self._inchikey_from_rdkit = 'ERROR'
        return self._inchikey_from_rdkit

    def inchikey_from_rdkit_matches_ccd(self, connectivity_only=False):
        """
        checks whether inchikey matches between ccd and rdkit
        
        Args:
            connectivity_only (bool): restrict to the first 14 character - the connectivity information.

        Returns:
            bool: True for match 

        """
        if self.inchikey is None or self.inchikey_from_rdkit == 'ERROR':
            return False
        if connectivity_only:
            if len(self.inchikey) < 14 or len(self.inchikey_from_rdkit) < 14:
                return False
            elif self.inchikey[:14] != self.inchikey_from_rdkit[:14]:
                return False
        elif self.inchikey != self.inchikey_from_rdkit:
            return False
        return True

    def sdf_file_or_string(self, file_name=None, ideal=True, hydrogen=True, alias=False,
                           xyz=None, raise_exception=False):
        """
        write a sdf file or return a string containing the molecule as a sdf file

        Args:
            file_name (str): optional filename
            ideal (bool): write the ideal coordinates (True) or model coordinates (False)? Default True: ideal.
            hydrogen (bool): include hydrogen atoms in the sdf? Default True: yes)
            alias (bool): use the alias feature to include atom names in the sdf? Default False no.
            xyz: list of (x,y,z) coordinates to be written instead of ideal or model
            raise_exception (bool): raise an exception on RDKit problems. Defaults to silently returning

        Returns:
            None or a string containing the molecule converted to sdf
        """
        fname = self.chem_comp_id + '.sdf'
        if xyz is not None:
            conformer_id = self.load_xyz_into_conformer(their_xyz=xyz)
        elif ideal:
            conformer_id = self.conformer_id_ideal
        else:
            conformer_id = self.conformer_id_model
        if hydrogen:
            mol_for_sdf = self.rwmol_original
        else:
            mol_for_sdf = self.rwmol_original_remove_h
        self.__sdf_alias_on_off(mol_for_sdf, alias=alias)

        try:
            sdf_string = Chem.MolToMolBlock(mol_for_sdf, confId=conformer_id)
        except Exception:
            if raise_exception:
                raise  # re-raise
            else:
                return
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
        # noinspection PyArgumentList
        atom_pdb_residue_info = Chem.rdchem.AtomPDBResidueInfo()
        atom_pdb_residue_info.SetResidueName(self.chem_comp_id)
        atom_pdb_residue_info.SetTempFactor(20.0)
        atom_pdb_residue_info.SetOccupancy(1.0)
        atom_pdb_residue_info.SetChainId('A')
        atom_pdb_residue_info.SetResidueNumber(1)
        atom_pdb_residue_info.SetIsHeteroAtom(True)
        for atom in self.rwmol_original.GetAtoms():
            atom_index = atom.GetIdx()
            pdbx_align = self.atom_pdbx_align[atom_index]
            atom_name = self.atom_ids[atom_index]
            if pdbx_align == '1':
                atom_name = ' ' + atom_name
            atom_name = '{:<4}'.format(atom_name)  # make sure it is 4 characters
            atom_pdb_residue_info.SetName(atom_name)
            atom.SetMonomerInfo(atom_pdb_residue_info)
        if ideal:
            conformer_id = self.conformer_id_ideal
        else:
            conformer_id = self.conformer_id_model
        ideal_or_model = 'ideal' if ideal else 'model'
        pdb_title = 'TITLE     {} coordinates'.format(ideal_or_model)
        pdb_title += ' for PDB Chemical Component Definition {}\n'.format(self.chem_comp_id)
        pdb_title += 'AUTHOR    ccd_utils using RDKit\n'
        pdb_string = pdb_title + Chem.MolToPDBBlock(self.rwmol_original, conformer_id)
        if file_name is None:
            return pdb_string
        else:
            with open(file_name, "w") as pdb_file:
                pdb_file.write(pdb_string)
        return None

    def cml_file_or_string(self, file_name=None):
        top = etree.Element('cml')
        top.set('dictRef', 'ebiMolecule:ebiMoleculeDict.cml')
        top.set('ebiMolecule', 'http://www.ebi.ac.uk/felics/molecule')
        f_charge = 0
        for atom_index in range(self.number_atoms):
            charge = self.atom_charges[atom_index]
            if charge is None:
                charge = 0
            f_charge += charge
        mol = etree.SubElement(top, 'molecule', id=self.chem_comp_id, formalCharge=str(f_charge))
        id_inchi = etree.SubElement(mol, 'identifier', dictRef='ebiMolecule:inchi')
        id_inchi.text = self.inchikey
        id_systematic = etree.SubElement(mol, 'identifier', dictRef='ebiMolecule:systematicName')
        id_systematic.text = self.chem_comp_name
        id_formula1 = etree.SubElement(mol, 'formula', dictRef="ebiMolecule:stereoSmiles")
        id_formula2 = etree.SubElement(mol, 'formula', dictRef="ebiMolecule:nonStereoSmiles")
        id_formula1.text = self.smiles_canonical_cactvs
        id_formula2.text = self.smiles_acdlabs
        atom_array = etree.SubElement(mol, 'atomArray')
        for atom_index in range(self.number_atoms):
            element = self.atom_elements[atom_index]
            atom_name = self.atom_ids[atom_index]
            atom_entry = etree.SubElement(atom_array, 'atom', id=atom_name, elementType=element)
            (ideal_x, ideal_y, ideal_z) = self.replace_none_with_0_0_0(self.ideal_xyz[atom_index])
            atom_entry.set('x3', str(ideal_x))
            atom_entry.set('y3', str(ideal_y))
            atom_entry.set('z3', str(ideal_z))
        bond_array = etree.SubElement(mol, 'bondArray')
        for bond_index in range(self.number_bonds):
            atom_1 = self.bond_atom_name_1[bond_index]
            atom_2 = self.bond_atom_name_2[bond_index]
            bond_order = Chem.rdchem.BondType(self.bond_order[bond_index])
            # print bond_order
            bond_entry = etree.SubElement(bond_array, 'bond')
            bond_entry.set('atomsRefs2', atom_1+' '+atom_2)
            bond_entry.set('order', str(bond_order))
        # tree = etree.ElementTree(top)
        cml_string = etree.tostring(top, pretty_print=True)
        if file_name is None:
            return cml_string
        else:
            with open(file_name, 'wb') as cml_file:
                cml_file.write(cml_string)
                cml_file.close()
            return None

    def xyz_file_or_string(self, file_name=None, ideal=True):
        """
        write a xyz format file or return a string containing the molecule in xyz format

        Args:
            file_name (str): optional filename
            ideal (bool): write the ideal coordinates (True) or model coordinates (False)? Default True: ideal.

        Returns:
            None or a string containing the molecule converted to xyz
        """
        ret_str = '{}\n'.format(self.number_atoms)
        ret_str += '{}\n'.format(self.chem_comp_id)
        for atom_index in range(self.number_atoms):
            element = self.atom_elements[atom_index]
            if ideal:
                (x, y, z) = self.replace_none_with_0_0_0(self.ideal_xyz[atom_index])
            else:
                (x, y, z) = self.replace_none_with_0_0_0(self.model_xyz[atom_index])
            ret_str += '{: <2}{:9.4f} {:9.4f} {:9.4f}\n'.format(element, x, y, z)
        if file_name is None:
            return ret_str
        else:
            with open(file_name, 'w') as file_out:
                file_out.write(ret_str)
                file_out.close()
            return None

    def image_file_or_string(self, file_name=None, wedge=True, atom_labels=True, hydrogen=False,
                             pixels_x=400, pixels_y=200, highlight_bonds=None, black=False,
                             raise_exception=False):
        """
        produces a svg image of the molecule to a string or file using RDKit.

        Args:
            file_name (str): the name of the output svg file or None to return the svg as a string.
            wedge (bool):  wedge the bonds in the image
            atom_labels (bool): include atom labels in the image
            hydrogen (bool): include hydrogen atoms in the image.
            pixels_x (int): size of image in pixels
            pixels_y (int): size of image in pixels
            highlight_bonds: an ordered dictionary of bonds to highlight key (atom_index_0, atom_index_1) to (r,g,b)
            black (bool): wipe out atom colors and make the molecular diagram black (highlight_bonds are not affected).
            raise_exception (bool): raise an exception on RDKit problems. Defaults to silently returning

        Returns:
            None or a string containing the svg string of the molecule.
        """
        if self.rwmol_cleaned is None:
            logging.warning('using original rwmol for {} image generation as clean up/sanitize failed'.
                            format(self.chem_comp_id))
            if hydrogen:
                mol_to_draw = self.rwmol_original
            else:
                mol_to_draw = self.rwmol_original_remove_h
        else:
            if hydrogen:
                mol_to_draw = self.rwmol_cleaned
            else:
                mol_to_draw = self.rwmol_cleaned_remove_h


        copy_to_draw = Chem.RWMol(mol_to_draw, True)  # boolean quickCopy means no conformers.
        # load 2D coords
        if self.__load_2d_coordinates_in_rwmol(copy_to_draw) is None:
            logging.error('Problem for {} cannot produce 2D coords.'.format(self.chem_comp_id))
        drawer = rdMolDraw2D.MolDraw2DSVG(pixels_x, pixels_y)
        # noinspection PyArgumentList
        opts = drawer.drawOptions()
        if atom_labels:
            # noinspection PyArgumentList
            for atom in copy_to_draw.GetAtoms():
                atom_index = atom.GetIdx()
                atom_name = self.atom_ids[atom_index]
                opts.atomLabels[atom_index] = atom_name
        try:
            copy_to_draw = rdMolDraw2D.PrepareMolForDrawing(copy_to_draw, wedgeBonds=wedge, addChiralHs=False)
        except Exception as e_mess:
            if raise_exception:
                raise  # re-raise
            else:
                logging.error('Problem for {} in PrepareMolForDrawing: {}'.format(self.chem_comp_id, e_mess))
                return
        if highlight_bonds is None:
            drawer.DrawMolecule(copy_to_draw)
        else:
            # highlight the atoms on each end of the bond in the colour
            highlight_atoms_colours = {}
            for key in highlight_bonds:
                highlight_atoms_colours[key[0]] = highlight_bonds[key]
                highlight_atoms_colours[key[1]] = highlight_bonds[key]
            # find the bond in the rdkit molecule
            highlight_bonds_colours = {}
            for key in highlight_bonds:
                at_index_0 = key[0]
                at_index_1 = key[1]
                for bond in copy_to_draw.GetBonds():
                    bond_index_0 = bond.GetBeginAtomIdx()
                    bond_index_1 = bond.GetEndAtomIdx()
                    if (at_index_0 == bond_index_0 and at_index_1 == bond_index_1) or \
                       (at_index_0 == bond_index_1 and at_index_1 == bond_index_0):
                        highlight_bonds_colours[bond.GetIdx()] = highlight_bonds[key]
                        break
            drawer.DrawMolecule(copy_to_draw,
                                highlightAtoms=highlight_atoms_colours.keys(),
                                highlightAtomColors=highlight_atoms_colours,
                                highlightBonds=highlight_bonds_colours.keys(),
                                highlightBondColors=highlight_bonds_colours)
        # noinspection PyArgumentList
        drawer.FinishDrawing()
        # noinspection PyArgumentList
        svg = drawer.GetDrawingText()
        # next line might be needed to get svg to display in Jupyter notebook?
        # svg = svg.replace('svg:', '')
        if black:
            svg = svg.replace('#FF0000', '#000000')  # get rid of red oxygen
            svg = svg.replace('#0000FF', '#000000')  # get rid of blue nitrogen
            svg = svg.replace('#FF7F00', '#000000')  # get rid of orange phosphorous
            svg = svg.replace('#CCCC00', '#000000')  # get rid of yellow sulphur
        if file_name is None:
            return svg
        else:
            with open(file_name, 'w') as img_file:
                img_file.write(svg)
                img_file.close()
        return None

    @staticmethod
    def replace_none_with_0_0_0(xyz_coordinates_tuple):
        """
        deals with missing coordinates - replacing None values with (0., 0., 0.)

        Args:
            xyz_coordinates_tuple: tuple (x, y, z) or None

        Returns:
            tuple (x, y, z)

        """
        if xyz_coordinates_tuple is None:
            return 0., 0., 0.
        else:
            return xyz_coordinates_tuple
