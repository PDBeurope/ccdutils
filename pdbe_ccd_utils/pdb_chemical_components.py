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
import os
from collections import namedtuple, OrderedDict


class PdbChemicalComponents(object):
    """
    deals with parsing the PDB chemical chemical component cif file.

    Notes
    Curently limited to parsing in the first cif definition in a file. Does not (yet) deal with multiple.
    """
    _chem_comp_atom_items = ('comp_id',
                             'atom_id',
                             'alt_atom_id',
                             'type_symbol',
                             'charge',
                             'pdbx_align',
                             'pdbx_aromatic_flag',
                             'pdbx_leaving_atom_flag',
                             'pdbx_stereo_config',
                             'model_Cartn_x',
                             'model_Cartn_y',
                             'model_Cartn_z',
                             'pdbx_model_Cartn_x_ideal',
                             'pdbx_model_Cartn_y_ideal',
                             'pdbx_model_Cartn_z_ideal',
                             'pdbx_component_atom_id',
                             'pdbx_component_comp_id',
                             'pdbx_ordinal')
    """list of the items used in _chem_comp_atom"""

    Bond = namedtuple('Bond', 'atom_id_1 atom_id_2 value_order pdbx_aromatic_flag pdbx_stereo_config')

    def __init__(self, file_name=None, cif_dictionary=None, cif_parser='auto'):
        """
        initializer - creates a PdbChemicalComponents object normally from a cif file

        Args:
            file_name (str): the filename of the PDB CCD file
            cif_dictionary: an ordered dictionary produced by PDBeCIF parser
            cif_parser (str): the cif parser to use. One of 'auto' or 'PDBeCIF' or'CifFile' or 'test_hard_code_cmo'
        """
        self.chem_comp_id = None
        self.chem_comp_name = None
        self.chem_comp_formula = None
        self.chem_comp_pdbx_release_status = None
        self.smiles_acdlabs = None
        self.smiles_canonical_cactvs = None
        self.smiles_cactvs = None
        self.smiles_canonical_openeye = None
        self.smiles_openeye = None
        self.inchi = None
        self.inchikey = None
        self.systematic_name_openeye = None
        self.systematic_name_acdlabs = None
        self._atoms = []
        """list of ordered dictionary"""
        self.__atom_ids = None
        self.__elements = None
        self.__stereo_configs = None
        self.__charges = None
        self.__ideal_xyz = None
        self.__model_xyz = None
        self.__pdbx_align = None
        self.bonds = []
        self.bond_atom_index_1 = []
        """list of int: one for each of self.bonds the index of the matching atom_id_1 in self.atoms"""
        self.bond_atom_index_2 = []
        """list of int: one for each of self.bonds the index of the matching atom_id_2 in self.atoms"""
        self.bond_atom_name_1 = []
        self.bond_atom_name_2 = []
        self.bond_order = []
        """list of int: one for each of self.bonds the bond order for the bond got from self.bonds value_order"""
        self.bond_aromatic = []
        """list of bool: one for each of self.bonds boolean conversion of pdbx_aromatic_flag (Y or N)"""
        self.cif_parser = cif_parser
        self.pdbecif_cif_obj = None
        if cif_parser == 'test_hard_code_cmo':
            self.__load_carbon_monoxide_hard_coded()
        elif file_name is not None:
            self.read_ccd_from_cif_file(file_name)
            self.setup_bond_lists()
        elif cif_dictionary is not None:
            self.read_ccd_from_pdbecif_cif_dictionary(cif_dictionary)
            self.setup_bond_lists()

    @staticmethod
    def empty_chem_comp_atom():
        """
        supply an empty chem_comp_atom - all items set to None

        Returns:
            OrderedDict: of items found in _chem_comp_atom_items
        """
        return OrderedDict([(k, None) for k in PdbChemicalComponents._chem_comp_atom_items])

    @property
    def number_atoms(self):
        """
        The number of atoms in the chem_comp

        Returns:
            int: the number of atoms
        """
        return len(self._atoms)

    @property
    def atom_ids(self):
        """
        tuple of the atom_id's (aka atom names) in the chem_comp

        Returns:
            (str): the atom_id's
        """
        if self.__atom_ids is None:
            self.__atom_ids = []
            for atom in self._atoms:
                self.__atom_ids.append(atom['atom_id'])
            self.__atom_ids = tuple(self.__atom_ids)
        return self.__atom_ids

    @property
    def atom_elements(self):
        """
        the elements for the atoms in the chem_comp_atom list

        Returns:
            (str): the elements for each atom
        """
        if self.__elements is None:
            self.__elements = []
            for atom in self._atoms:
                type_symbol = atom['type_symbol']
                if type_symbol is None or len(type_symbol) == 0:
                    raise RuntimeError('chem_comp_atom invalid type_symbol={}'.format(type_symbol))
                element = type_symbol[0]
                if len(type_symbol) > 1:
                    element += type_symbol[1].lower()
                self.__elements.append(element)
            self.__elements = tuple(self.__elements)
        return self.__elements

    @property
    def atom_stereo_configs(self):
        """
        the pdbx_stereo_config for the atoms in the chem_comp_atom list

        Returns:
            (str): the pdbx_stereo_config for each atom

        Raises:
            TODO custom exception if one of atom's pdbx_stereo_config is not 'N', 'R' or 'S'
        """
        if self.__stereo_configs is None:
            self.__stereo_configs = []
            for atom in self._atoms:
                self.__stereo_configs.append(atom['pdbx_stereo_config'])
            self.__stereo_configs = tuple(self.__stereo_configs)
        return self.__stereo_configs

    @property
    def atom_charges(self):
        """
        the formal charges for the atoms in the chem_comp_atom list

        Returns:
            tuple[int]: the chem_comp.charge value for each atom or None if there is a conversion error.

        """
        if self.__charges is None:
            self.__charges = []
            for atom in self._atoms:
                try:
                    charge = int(atom['charge'])
                except ValueError:
                    charge = None
                self.__charges.append(charge)
            self.__charges = tuple(self.__charges)
        return self.__charges

    @property
    def atom_charges_missing_values(self):
        if None in self.atom_charges:
            return True
        else:
            return False

    @property
    def atom_pdbx_align(self):
        """
        the pdbx_align for the atoms in the chem_comp_atom list

        Returns:
            (int): the chem_comp.pdbx_align for each atom
        """
        if self.__pdbx_align is None:
            self.__pdbx_align = []
            for atom in self._atoms:
                pdbx_align = atom['pdbx_align']
                self.__pdbx_align.append(pdbx_align)
        self.__pdbx_align = tuple(self.__pdbx_align)
        return self.__pdbx_align

    @property
    def number_bonds(self):
        """
        The number of bonds in the chem_comp

        Returns:
            int: the number of bonds
        """
        return len(self.bonds)

    @property
    def ideal_xyz(self):
        """
        The ideal coordinates from chem_comp.pdbx_model_Cartn_x_ideal, chem_comp.pdbx_model_Cartn_y _ideal,
        chem_comp.pdbx_model_Cartn_z_ideal,

        Returns:
            tuple of tuple( x, y, z) for each atom. x, y, z are floats
        """
        if self.__ideal_xyz is None:
            self.__ideal_xyz = self.__supply_model_or_ideal_coords(ideal=True)
        return self.__ideal_xyz

    @property
    def ideal_xyz_has_missing_values(self):
        if None in self.ideal_xyz:
            return True
        else:
            return False

    @property
    def model_xyz(self):
        """
        The model coordinates from chem_comp.model_Cartn_x, chem_comp.model_Cartn_y,
        chem_comp.model_Cartn_z,

        Returns:
            tuple of tuple( x, y, z) for each atom. x, y, z are floats

        Notes:
            if there is a conversion error returns coordinates None for the atom in place of x, y, z
        """
        if self.__model_xyz is None:
            self.__model_xyz = self.__supply_model_or_ideal_coords(ideal=False)
        return self.__model_xyz

    @property
    def model_xyz_has_missing_values(self):
        if None in self.model_xyz:
            return True
        else:
            return False

    def __supply_model_or_ideal_coords(self, ideal=True):
        """
        to avoid code duplicationin ideal_xyz and model_xyz

        Args:
            ideal (bool): ideal or model?

        Returns:
            tuple of tuple( x, y, z) for each atom. x, y, z are floats

        Notes:
            if there is a conversion error returns coordinates None for the atom in place of x, y, z
        """
        this_xyz = []
        if ideal:
            xyz_fields = ('pdbx_model_Cartn_x_ideal', 'pdbx_model_Cartn_y_ideal', 'pdbx_model_Cartn_z_ideal')
        else:
            xyz_fields = ('model_Cartn_x', 'model_Cartn_y', 'model_Cartn_z')
        for atom in self._atoms:
            try:
                float_x = float(atom[xyz_fields[0]])
                float_y = float(atom[xyz_fields[1]])
                float_z = float(atom[xyz_fields[2]])
                this_xyz.append((float_x, float_y, float_z))
            except ValueError:
                this_xyz.append(None)
        this_xyz = tuple(this_xyz)
        return this_xyz

    def __load_carbon_monoxide_hard_coded(self):
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
        # _chem_comp.chem_comp_pdbx_release_status                   REL
        self.chem_comp_pdbx_release_status = 'REL'
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
        my_chem_comp_atom = self.empty_chem_comp_atom()
        my_chem_comp_atom['atom_id'] = 'C'
        my_chem_comp_atom['type_symbol'] = 'C'
        my_chem_comp_atom['pdbx_stereo_config'] = 'N'
        my_chem_comp_atom['pdbx_model_Cartn_x_ideal'] = '0.607'
        my_chem_comp_atom['pdbx_model_Cartn_y_ideal'] = '0.000'
        my_chem_comp_atom['pdbx_model_Cartn_z_ideal'] = '0.000'
        my_chem_comp_atom['model_Cartn_x'] = '-0.296'
        my_chem_comp_atom['model_Cartn_y'] = '8.526'
        my_chem_comp_atom['model_Cartn_z'] = '17.112'
        my_chem_comp_atom['charge'] = '-1'
        self._atoms.append(my_chem_comp_atom)
        my_chem_comp_atom = self.empty_chem_comp_atom()
        my_chem_comp_atom['atom_id'] = 'O'
        my_chem_comp_atom['type_symbol'] = 'O'
        my_chem_comp_atom['pdbx_stereo_config'] = 'N'
        my_chem_comp_atom['pdbx_model_Cartn_x_ideal'] = '-0.600'
        my_chem_comp_atom['pdbx_model_Cartn_y_ideal'] = '0.000'
        my_chem_comp_atom['pdbx_model_Cartn_z_ideal'] = '0.000'
        my_chem_comp_atom['model_Cartn_x'] = '0.023'
        my_chem_comp_atom['model_Cartn_y'] = '7.997'
        my_chem_comp_atom['model_Cartn_z'] = '18.053'
        my_chem_comp_atom['charge'] = '1'
        self._atoms.append(my_chem_comp_atom)
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
        self.bond_atom_name_1 = []
        self.bond_atom_name_2 = []
        self.bond_order = []
        self.bond_aromatic = []
        for bond in self.bonds:
            atom_id_1 = bond.atom_id_1
            index_atom_1 = self.find_atom_index(atom_id_1)
            self.bond_atom_index_1.append(index_atom_1)
            self.bond_atom_name_1.append(atom_id_1)
            atom_id_2 = bond.atom_id_2
            index_atom_2 = self.find_atom_index(atom_id_2)
            self.bond_atom_index_2.append(index_atom_2)
            self.bond_atom_name_2.append(atom_id_2)
            bond_order = self.map_value_order_to_int(bond.value_order)
            if bond_order == -1:
                raise RuntimeError('problem with bond order for bond {}'.format(bond))
            self.bond_order.append(bond_order)
            if bond.pdbx_aromatic_flag == 'Y':
                bond_aromatic = True
            else:
                bond_aromatic = False
            self.bond_aromatic.append(bond_aromatic)

    def find_atom_index(self, atom_id):
        for index in range(len(self._atoms)):
            this_atom = self._atoms[index]
            if atom_id == this_atom['atom_id']:
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
            raise ValueError('cannot read PDB chemical components from {} as file not found'.format(file_name))
        if self.cif_parser == 'auto':
            try:
                self.read_ccd_from_file_pdbecif(file_name)
            except ImportError:
                self.read_ccd_from_file_ciffile(file_name)
        elif self.cif_parser == 'PDBeCIF':
            self.read_ccd_from_file_pdbecif(file_name)
        elif self.cif_parser == 'CifFile':
            self.read_ccd_from_file_ciffile(file_name)
        else:
            raise RuntimeError('unrecognized cif_parser {}'.format(self.cif_parser))

    def read_ccd_from_file_pdbecif(self, file_name):
        """
        reads the chemical component from file file_name using the pdbecif parser
        https://github.com/glenveegee/PDBeCIF.git

        Args:
            file_name (str): the filename

        Returns:
            None

        Raises:
            ImportError: if the parser cannot be loaded.
            RuntimeError: if a new unrecognized item has appeared
        """
        import mmCif.mmcifIO as mmcifIO
        cif_parser = mmcifIO.CifFileReader(input='data', preserve_order=True)
        self.pdbecif_cif_obj = cif_parser.read(file_name, output='cif_wrapper')
        self._read_ccd_from_pdbecif_cif_obj()

    def read_ccd_from_pdbecif_cif_dictionary(self, mmcif_dict):
        """
        reads the chemical component from a 'cif_dictionary'

        Args:
            mmcif_dict: PDBeCIF cif dictionary - normally an ordered dictionary

        Returns:
            None
        """
        from mmCif import CIFWrapper
        token_ordering = True
        # next line taken from CifFileReader method read
        self.pdbecif_cif_obj = \
            dict(((block_id, CIFWrapper(block_data, data_id=block_id, preserve_token_order=token_ordering))
                  for block_id, block_data in list(mmcif_dict.items())))
        self._read_ccd_from_pdbecif_cif_obj()

    def _read_ccd_from_pdbecif_cif_obj(self):
        """
        reads the chemical component from a PDBeCIF cif object.

        Returns:
            None
        """
        data_block = list(self.pdbecif_cif_obj.values())[0]
        # noinspection PyProtectedMember
        chem_comp = data_block._chem_comp
        for thing in 'id', 'name', 'formula', 'pdbx_release_status':
            value = chem_comp[thing][0]
            setattr(self, "chem_comp_" + thing, value)
        self._atoms = []
        # noinspection PyProtectedMember
        chem_comp_atom = data_block._chem_comp_atom
        empty_atom = self.empty_chem_comp_atom()
        if chem_comp_atom is not None:
            for atom in chem_comp_atom:
                self._atoms.append(atom)
                # check the no new attributes have been set
                for key in atom:
                    if key not in empty_atom:
                        raise RuntimeError('unrecognized item "{}" in chem_comp_atom'.format(key))
        self.bonds = []
        # noinspection PyProtectedMember
        chem_comp_bond = data_block._chem_comp_bond
        if chem_comp_bond is not None:
            for bond in chem_comp_bond:
                atom_id_1 = bond['atom_id_1']
                atom_id_2 = bond['atom_id_2']
                value_order = bond['value_order']
                pdbx_aromatic_flag = bond['pdbx_aromatic_flag']
                pdbx_stereo_config = bond['pdbx_stereo_config']
                this_bond = self.Bond(atom_id_1=atom_id_1, atom_id_2=atom_id_2, value_order=value_order,
                                      pdbx_aromatic_flag=pdbx_aromatic_flag, pdbx_stereo_config=pdbx_stereo_config)
                self.bonds.append(this_bond)
        self._pdbecif_parse_pdbx_chem_comp_descriptor(data_block)
        self._pdbecif_parse_pdbx_chem_comp_identifier(data_block)

    def _pdbecif_parse_pdbx_chem_comp_descriptor(self, data_block):
        """
        parses contents of _pdbx_chem_comp_descriptor block, SMILES strings &  inchi stuff

        Args:
            data_block: PDBeCIF datablock obtained from PDB-CCD

        Returns:
            None

        """
        # noinspection PyProtectedMember
        pdbx_chem_comp_descriptor = data_block._pdbx_chem_comp_descriptor
        if pdbx_chem_comp_descriptor is None:
            pass
        else:
            for descriptor in pdbx_chem_comp_descriptor:
                if descriptor['type'] == 'SMILES' and descriptor['program'] == 'ACDLabs':
                    self.smiles_acdlabs = descriptor['descriptor']
                elif descriptor['type'] == 'SMILES_CANONICAL' and descriptor['program'] == 'CACTVS':
                    self.smiles_canonical_cactvs = descriptor['descriptor']
                elif descriptor['type'] == 'SMILES' and descriptor['program'] == 'CACTVS':
                    self.smiles_cactvs = descriptor['descriptor']
                elif descriptor['type'] == 'SMILES_CANONICAL' and 'OpenEye' in descriptor['program']:
                    self.smiles_canonical_openeye = descriptor['descriptor']
                elif descriptor['type'] == 'SMILES' and 'OpenEye' in descriptor['program']:
                    self.smiles_openeye = descriptor['descriptor']
                elif descriptor['type'] == 'InChI':
                    self.inchi = descriptor['descriptor']
                elif descriptor['type'] == 'InChIKey':
                    self.inchikey = descriptor['descriptor']
                else:
                    logging.warn('unrecognized pdbx_chem_comp_descriptor {}'.format(descriptor))

    def _pdbecif_parse_pdbx_chem_comp_identifier(self, data_block):
        """
        parses contents of _pdbx_chem_comp_identifier block to obtain systematic names

        Args:
            data_block: PDBeCIF datablock obtained from PDB-CCD

        Returns:

        """
        # noinspection PyProtectedMember
        pdbx_chem_comp_identifier = data_block._pdbx_chem_comp_identifier
        if pdbx_chem_comp_identifier is None:
            pass
        else:
            for identifier in pdbx_chem_comp_identifier:
                if identifier['type'] == 'SYSTEMATIC NAME' and identifier['program'] == 'ACDLabs':
                    self.systematic_name_acdlabs = identifier['identifier']
                elif identifier['type'] == 'SYSTEMATIC NAME' and 'OpenEye' in identifier['program']:
                    self.systematic_name_openeye = identifier['identifier']
                else:
                    logging.warn('unrecognized chem_comp_identifier {} '.format(identifier))

    def write_ccd_cif(self, output_ccd_cif_file_name):
        """
        writes out the read in pdb_ccd as a cif file.

        Args:
            output_ccd_cif_file_name (str): file name for the output

        Returns:
            None

        Notes:
            Currently limited to PDBeCIF parser. Currently output is same as input.
        """
        if self.pdbecif_cif_obj is None:
            raise NotImplementedError('cannot write_ccd_cif currently only supported if using PDBeCIF parser')
        import mmCif.mmcifIO as mmcifIO
        cif_file_writer = mmcifIO.CifFileWriter(output_ccd_cif_file_name, preserve_order=True)
        cif_wrapper = list(self.pdbecif_cif_obj.values())[0]
        cif_file_writer.write(cif_wrapper)

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
        # noinspection PyUnresolvedReferences
        from pdbx_v2.core.CifFile import CifFile
        # method based on calls made by
        # https://svn-dev.wwpdb.org/svn-wwpdb/py-validation/trunk/src/python/pdboi/pdbdata/mmcifapiconnector.py
        cif_file = CifFile(file_name, parseLogFileName=None).getCifFile()
        first_data_block = cif_file.GetBlock(cif_file.GetFirstBlockName())
        table_chem_comp = first_data_block.GetTable('chem_comp')
        for thing in 'id', 'name', 'pdbx_release_status':
            value = table_chem_comp(0, thing)
            setattr(self, "chem_comp_" + thing, value)
        self._atoms = []
        table_chem_comp_atom = first_data_block.GetTable('chem_comp_atom')
        number_atoms = table_chem_comp_atom.GetNumRows()
        for row_num in range(number_atoms):
            my_chem_comp_atom = self.empty_chem_comp_atom()
            for key in my_chem_comp_atom:
                my_chem_comp_atom[key] = table_chem_comp_atom(row_num, key)
            self._atoms.append(my_chem_comp_atom)
        table_pdbx_chem_comp_descriptor = first_data_block.GetTable('pdbx_chem_comp_descriptor')
        for row_num in range(table_pdbx_chem_comp_descriptor.GetNumRows()):
            if table_pdbx_chem_comp_descriptor(row_num, 'type') == 'InChIKey':
                self.inchikey = table_pdbx_chem_comp_descriptor(row_num, 'descriptor')
        self.bonds = []
        table_chem_comp_bond = first_data_block.GetTable('chem_comp_bond')
        number_bonds = table_chem_comp_bond.GetNumRows()
        for row_num in range(number_bonds):
            atom_id_1 = table_chem_comp_bond(row_num, 'atom_id_1')
            atom_id_2 = table_chem_comp_bond(row_num, 'atom_id_2')
            value_order = table_chem_comp_bond(row_num, 'value_order')
            pdbx_aromatic_flag = table_chem_comp_bond(row_num, 'pdbx_aromatic_flag')
            pdbx_stereo_config = table_chem_comp_bond(row_num, 'pdbx_stereo_config')
            this_bond = self.Bond(atom_id_1=atom_id_1, atom_id_2=atom_id_2, value_order=value_order,
                                  pdbx_aromatic_flag=pdbx_aromatic_flag, pdbx_stereo_config=pdbx_stereo_config)
            self.bonds.append(this_bond)

    def __eq__(self, other):
        if type(other) is not type(self):
            return False
        for key, value in self.__dict__.items():
            if key == 'pdbecif_cif_obj':
                pass
            else:
                other_value = other.__dict__[key]
                if value != other_value:
                    return False
        return True

    def __ne__(self, other):
        return not self.__eq__(other)
