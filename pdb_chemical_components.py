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
from collections import namedtuple


class PdbChemicalComponents(object):
    """
    deals with parsing the PDB chemical chemical component cif file
    currently limited to a single cif item in the dictionary.
    """

    def __init__(self, file_name=None, cif_parser='auto'):
        """

        Args:
            file_name (str): filename
            cif_parser (str): the cif parser to use. One of 'auto' or 'mmcifIO'(EBI) or'CifFile'
        """
        self.chem_comp_id = None
        self.chem_comp_name = None
        self.pdbx_release_status = None
        self.inchikey = None
        self.Atom = namedtuple('Atom', 'atom_id pdbx_stereo_config xyz_ideal')
        self.atoms = []
        self.Bond = namedtuple('Bond', 'atom_id_1 atom_id_2 value_order pdbx_aromatic_flag pdbx_stereo_config')
        self.bonds = []
        self.bond_atom_index_1 = []
        """list of int: one for each of self.bonds the index of the matching atom_id_1 in self.atoms"""
        self.bond_atom_index_2 = []
        """list of int: one for each of self.bonds the index of the matching atom_id_2 in self.atoms"""
        self.bond_order = []
        """list of int: one for each of self.bonds the bond order for the bond got from self.bonds value_order"""
        self.bond_aromatic = []
        """list of bool: one for each of self.bonds boolean conversion of pdbx_aromatic_flag (Y or N)"""
        self.cif_parser = cif_parser
        if file_name is not None:
            self.read_ccd_from_cif_file(file_name)
            self.setup_bond_lists()

    @property
    def atom_ids(self):
        """
        tuple of the atom_id's (aka atom names) in the chem_comp

        Returns:
            (str): the atom_id's
        """
        atom_ids = []
        for atom in self.atoms:
            atom_ids.append(atom.atom_id)
        return tuple(atom_ids)

    @property
    def number_atoms(self):
        """
        The number of atoms in the chem_comp

        Returns:
            int: the number of atoms
        """
        return len(self.atoms)

    @property
    def number_bonds(self):
        """
        The number of bonds in the chem_comp

        Returns:
            int: the number of bonds
        """
        return len(self.bonds)

    def load_carbon_monoxide_hard_coded(self):
        """
        stub to produce a hard coded carbon monoxide ccd object for development idea/testing
        without file parsing

        Returns:
            None
        """
        # _chem_comp.id                                    CMO
        self.chem_comp_id = 'CMO'
        # _chem_comp.name                                  "CARBON MONOXIDE"
        self.chem_comp_name = 'CARBON MONOXIDE'
        # _chem_comp.pdbx_release_status                   REL
        self.pdbx_release_status = 'REL'
        #
        # loop_
        # _pdbx_chem_comp_descriptor.comp_id 
        # _pdbx_chem_comp_descriptor.type 
        # _pdbx_chem_comp_descriptor.program 
        # _pdbx_chem_comp_descriptor.program_version 
        # _pdbx_chem_comp_descriptor.descriptor 
        # CMO SMILES           ACDLabs              10.04 "[O+]#[C-]"                 
        # CMO SMILES_CANONICAL CACTVS               3.341 "[C-]#[O+]"                 
        # CMO SMILES           CACTVS               3.341 "[C-]#[O+]"                 
        # CMO SMILES_CANONICAL "OpenEye OEToolkits" 1.5.0 "[C-]#[O+]"                 
        # CMO SMILES           "OpenEye OEToolkits" 1.5.0 "[C-]#[O+]"                 
        # CMO InChI            InChI                1.03  InChI=1S/CO/c1-2            
        # CMO InChIKey         InChI                1.03  UGFAIRIUMAVXCW-UHFFFAOYSA-N 
        self.inchikey = 'UGFAIRIUMAVXCW-UHFFFAOYSA-N'
        # CMO C C C -1 1 N N N -0.296 8.526 17.112 0.607  0.000 0.000 C CMO 1
        # CMO O O O 1  1 N N N 0.023  7.997 18.053 -0.600 0.000 0.000 O CMO 2
        this_atom = self.Atom(atom_id='C', pdbx_stereo_config='N', xyz_ideal=(0.607, 0.000, 0.000))
        self.atoms.append(this_atom)
        this_atom = self.Atom(atom_id='O', pdbx_stereo_config='N', xyz_ideal=(-0.600, 0.000, 0.000))
        self.atoms.append(this_atom)
        # _chem_comp_bond.comp_id              CMO
        # _chem_comp_bond.atom_id_1            C
        # _chem_comp_bond.atom_id_2            O
        # _chem_comp_bond.value_order          TRIP
        # _chem_comp_bond.pdbx_aromatic_flag   N
        # _chem_comp_bond.pdbx_stereo_config   N
        # _chem_comp_bond.pdbx_ordinal         1
        this_bond = self.Bond(atom_id_1='C', atom_id_2='O', value_order='TRIP', 
                              pdbx_aromatic_flag='N', pdbx_stereo_config='N')
        self.bonds.append(this_bond)
        self.setup_bond_lists()

    def setup_bond_lists(self):
        self.bond_atom_index_1 = []
        self.bond_atom_index_2 = []
        self.bond_order = []
        self.bond_aromatic = []
        for bond in self.bonds:
            atom_id_1 = bond.atom_id_1
            index_atom_1 = self.find_atom_index(atom_id_1)
            self.bond_atom_index_1.append(index_atom_1)
            atom_id_2 = bond.atom_id_2
            index_atom_2 = self.find_atom_index(atom_id_2)
            self.bond_atom_index_2.append(index_atom_2)
            bond_order = self.map_value_order_to_int(bond.value_order)
            if bond_order == -1:
                raise RuntimeError('problem with bond order for bond {}'.format(bond))
            self.bond_order.append(bond_order)

    def find_atom_index(self, atom_id):
        for index in range(len(self.atoms)):
            this_atom = self.atoms[index]
            if atom_id == this_atom.atom_id:
                return index
        return -1

    @staticmethod
    def map_value_order_to_int(value_order):
        if value_order == 'SING':
            return 1
        elif value_order == 'DOUB':
            return 2
        elif value_order == 'TRIP':
            return 3
        else:
            return -1

    def read_ccd_from_cif_file(self, file_name):
        """
        reads the ccd from a cif file

        Args:
            file_name (str): the filename

        Returns:
            None
        """
        if not os.path.isfile(file_name):
            raise ValueError('cannot read chemical compoents from %s as file not found' % file_name)
        if self.cif_parser == 'auto':
            try:
                self.read_ccd_from_file_mmcifio(file_name)
            except ImportError:
                self.read_ccd_from_file_ciffile(file_name)
        elif self.cif_parser == 'mmcifIO':
            self.read_ccd_from_file_mmcifio(file_name)
        elif self.cif_parser == 'CifFile':
            self.read_ccd_from_file_ciffile(file_name)
        else:
            raise RuntimeError('unrecognized cif_parser {}'.format(self.cif_parser))

    def read_ccd_from_file_mmcifio(self, file_name):
        """
        reads the chemical component from file file_name using the mmcifIO parser
        https://github.com/glenveegee/PDBeCIF.git

        Args:
            file_name (str): the filename

        Returns:
            None

        Raises:
            ImportError: if the parser cannot be loaded.
        """
        import mmCif.mmcifIO as mmcifIO
        cif_parser = mmcifIO.CifFileReader(input='data', preserve_order=True)
        cif_obj = cif_parser.read(file_name, output='cif_wrapper')
        data_block = list(cif_obj.values())[0]
        chem_comp = data_block._chem_comp
        self.chem_comp_id = chem_comp['id'][0]
        self.chem_comp_name = chem_comp['name'][0]
        self.pdbx_release_status = chem_comp['pdbx_release_status'][0]
        self.atoms = []
        chem_comp_atom = data_block._chem_comp_atom
        for atom in chem_comp_atom:
            atom_id = atom['atom_id']
            pdbx_stereo_config = atom['pdbx_stereo_config']
            ideal_x = float(atom['pdbx_model_Cartn_x_ideal'])
            ideal_y = float(atom['pdbx_model_Cartn_y_ideal'])
            ideal_z = float(atom['pdbx_model_Cartn_z_ideal'])
            this_atom = self.Atom(atom_id=atom_id,
                                  pdbx_stereo_config=pdbx_stereo_config,
                                  xyz_ideal=(ideal_x, ideal_y, ideal_z))
            self.atoms.append(this_atom)
        self.bonds = []
        chem_comp_bond = data_block._chem_comp_bond
        for bond in chem_comp_bond:
            atom_id_1 = bond['atom_id_1']
            atom_id_2 = bond['atom_id_2']
            value_order = bond['value_order']
            pdbx_aromatic_flag = bond['pdbx_aromatic_flag']
            pdbx_stereo_config = bond['pdbx_stereo_config']
            this_bond = self.Bond(atom_id_1=atom_id_1, atom_id_2=atom_id_2, value_order=value_order,
                                  pdbx_aromatic_flag=pdbx_aromatic_flag, pdbx_stereo_config=pdbx_stereo_config)
            self.bonds.append(this_bond)
        pdbx_chem_comp_descriptor = data_block._pdbx_chem_comp_descriptor
        for descriptor in pdbx_chem_comp_descriptor:
            if descriptor['type'] == 'InChIKey':
                self.inchikey = descriptor['descriptor']

    def read_ccd_from_file_ciffile(self, file_name):
        """
        reads the chemical component from file file_name using the pdbx_v2.core.CifFile parser

        Args:
            file_name (str): the filename

        Returns:
            None

        Raises:
            ImportError: if CifFile parser cannot be loaded.
        """
        from pdbx_v2.core.CifFile import CifFile
        # method based on calls made by
        # https://svn-dev.wwpdb.org/svn-wwpdb/py-validation/trunk/src/python/pdboi/pdbdata/mmcifapiconnector.py
        cif_file = CifFile(file_name, parseLogFileName=None).getCifFile()
        first_data_block = cif_file.GetBlock(cif_file.GetFirstBlockName())
        table_chem_comp = first_data_block.GetTable('chem_comp')
        self.chem_comp_id = table_chem_comp(0, 'id')
        self.chem_comp_name = table_chem_comp(0, 'name')
        self.pdbx_release_status = table_chem_comp(0, 'pdbx_release_status')
        self.atoms = []
        table_chem_comp_atom = first_data_block.GetTable('chem_comp_atom')
        number_atoms = table_chem_comp_atom.GetNumRows()
        for row_num in range(number_atoms):
            atom_id = table_chem_comp_atom(row_num, 'atom_id')
            pdbx_stereo_config = table_chem_comp_atom(row_num, 'pdbx_stereo_config')
            ideal_x = float(table_chem_comp_atom(row_num, 'pdbx_model_Cartn_x_ideal'))
            ideal_y = float(table_chem_comp_atom(row_num, 'pdbx_model_Cartn_y_ideal'))
            ideal_z = float(table_chem_comp_atom(row_num, 'pdbx_model_Cartn_z_ideal'))
            this_atom = self.Atom(atom_id=atom_id,
                                  pdbx_stereo_config=pdbx_stereo_config,
                                  xyz_ideal=(ideal_x, ideal_y, ideal_z))
            self.atoms.append(this_atom)
        table_pdbx_chem_comp_descriptor = first_data_block.GetTable('pdbx_chem_comp_descriptor')
        for row_num in range(table_pdbx_chem_comp_descriptor.GetNumRows()):
            if table_pdbx_chem_comp_descriptor(row_num, 'type') == 'InChIKey':
                self.inchikey = table_pdbx_chem_comp_descriptor(row_num, 'descriptor')
