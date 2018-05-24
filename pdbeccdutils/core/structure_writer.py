#!/usr/bin/env python
# software from PDBe: Protein Data Bank in Europe; https://pdbe.org
#
# Copyright 2018 EMBL - European Bioinformatics Institute
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

import copy
import math
import xml.etree.ElementTree as ET
from collections import OrderedDict
from xml.dom import minidom

import mmCif.mmcifIO as mmcif
import rdkit
from rdkit import Chem, rdBase

import pdbeccdutils
from pdbeccdutils.core import Component, ConformerType, ReleaseStatus


def write_molecule(path, component, remove_hs=True, conf_type=ConformerType.Ideal):
    """Export molecule in a specified format. Presently supported formats
    are: PDB CCD CIF (*.cif); Mol file (*.sdf); Chemical Markup language
    (*.cml); PDB file (*.pdb); XYZ file (*.xyz); XML (*.xml).
    ConformerType.AllConformers is presently supported only for PDB.

    Args:
        path (str): Path to the file. Extension determines format to be
            used.
        component (pdbeccdutils.core.Component): Component to be exported
        remove_hs (bool, optional): Defaults to True.
        conf_type (pdbeccdutils.core.ConformerType, optional):
            Defaults to ConformerType.Ideal.

    Raises:
        ValueError: For unsupported format
    """
    extension = path.split('.')[-1].lower()
    str_representation = ''

    if extension == 'sdf':
        str_representation = to_sdf_str(component, remove_hs, conf_type)
    elif extension == 'pdb':
        str_representation = to_pdb_str(component, remove_hs, conf_type)
    elif extension in ('mmcif', 'cif'):
        to_pdb_ccd_cif_file(path, component, remove_hs)
        return
    elif extension == 'cml':
        str_representation = to_cml_str(component, remove_hs, conf_type)
    elif extension == 'xml':
        str_representation = to_xml_str(component, remove_hs, conf_type)
    elif extension == 'xyz':
        str_representation = to_xyz_str(component, remove_hs, conf_type)
    else:
        raise ValueError('Unsupported file format: {}'.format(extension))

    with open(path, 'w') as f:
        f.write(str_representation)


def to_pdb_str(component, remove_hs=True, conf_type=ConformerType.Ideal):
    """Converts structure to the PDB format.

    Args:
        component (pdbeccdutils.core.Component): Component to be
            exported.
        remove_hs (bool, optional): Defaults to True.
        conf_type (pdbeccdutils.core.ConformerType, optional):
            Defaults to ConformerType.Ideal.

    Returns:
        str: String representation of the component in the PDB format.
    """
    (mol_to_save, conf_id, conf_type) = _prepate_structure(component, remove_hs, conf_type)

    info = Chem.rdchem.AtomPDBResidueInfo()
    info.SetResidueName(component.id)
    info.SetTempFactor(20.0)
    info.SetOccupancy(1.0)
    info.SetChainId('A')
    info.SetResidueNumber(1)
    info.SetIsHeteroAtom(True)

    for atom in mol_to_save.GetAtoms():
        if atom.HasProp('name'):
            flag = atom.GetProp('name')
            atom_name = '{:<4}'.format(flag)  # make sure it is 4 characters
            info.SetName(atom_name)
        atom.SetMonomerInfo(info)

    pdb_title = 'HEADER    {} coordinates'.format(conf_type.name)
    pdb_title += ' for PDB-CCD {}\n'.format(component.id)
    pdb_title += 'COMPND    {}\n'.format(component.id)
    pdb_title += 'AUTHOR    pdbccdutils {}\n'.format(pdbeccdutils.__version__)
    pdb_title += 'AUTHOR    RDKit {}\n'.format(rdkit.__version__)
    pdb_body = ''

    try:
        pdb_body = Chem.MolToPDBBlock(mol_to_save, conf_id)
    except Exception:
        pdb_body = _to_pdb_str_fallback(mol_to_save, conf_id, info)

    pdb_string = pdb_title + pdb_body

    return pdb_string


def to_sdf_str(component, remove_hs=True, conf_type=ConformerType.Ideal):
    """Converts structure to the SDF format.

    Args:
        component (pdbeccdutils.core.Component): Component to be
            exported.
        remove_hs (bool, optional): Defaults to True.
        conf_type (pdbeccdutils.core.ConformerType, optional):
            Defaults to ConformerType.Ideal.

    Returns:
        str: String representation of the component in the SDF format
    """
    (mol_to_save, conf_id, conf_type) = _prepate_structure(component, remove_hs, conf_type)

    mol_block = []
    mappings = []
    if conf_type == ConformerType.AllConformers:
        mappings = [ConformerType.Model, ConformerType.Ideal, ConformerType.Computed]
    else:
        mappings = [conf_type]

    try:
        for conf in mappings:
            try:
                s = '{} - {} conformer'.format(component.id, conf.name)
                s += Chem.MolToMolBlock(mol_to_save, confId=component.conformers_mapping[conf])
                s += '$$$$'
                mol_block.append(s)
            except ValueError:
                pass
    except Exception:
        mappings = {m.name: component.conformers_mapping[m] for m in mappings}
        mol_block = _to_sdf_str_fallback(mol_to_save, component.id, mappings)

    return "\n".join(mol_block)


def to_xyz_str(component, remove_hs=True, conf_type=ConformerType.Ideal):
    """Converts structure to the XYZ format. Does not yet support
    ConformerType.AllConformers.

    Args:
        component (pdbeccdutils.core.Component): Component to be
            exported.
        remove_hs (bool, optional): Defaults to True.
        conf_type (pdbeccdutils.core.ConformerType, optional):
            Defaults to ConformerType.Ideal.

    Returns:
        str: String representation of the component in the XYZ format
    """
    (mol_to_save, conf_id, conf_type) = _prepate_structure(component, remove_hs, conf_type)
    conformer = mol_to_save.GetConformer(id=conf_id)

    result = list()
    result.append(str(mol_to_save.GetNumAtoms()))
    result.append(component.id)

    for atom in mol_to_save.GetAtoms():
        coords = conformer.GetAtomPosition(atom.GetIdx())
        result.append('{0:<4}{1: f}  {2: f}  {3: f}'.
                      format(atom.GetSymbol(), coords.x, coords.y, coords.z))

    return '\n'.join(result)


def to_xml_xml(component, remove_hs=True, conf_type=ConformerType.Ideal):
    """Converts structure to the XML format and returns its XML repr.

    Args:
        component (pdbeccdutils.core.Component): Component to be
            exported.
        remove_hs (bool, optional): Defaults to True.
        conf_type (pdbeccdutils.core.ConformerType, optional):
            Defaults to ConformerType.Ideal.

    Returns:
        xml.etree.ElementTree.Element: XML object
    """
    root = ET.Element('chemComp')

    id_e = ET.SubElement(root, 'id')
    name_e = ET.SubElement(root, 'name')
    formula_e = ET.SubElement(root, 'formula')
    sys_name_e = ET.SubElement(root, 'systematicName')
    s_smiles_e = ET.SubElement(root, 'stereoSmiles')
    n_smiles_e = ET.SubElement(root, 'nonStereoSmiles')
    inchi_e = ET.SubElement(root, 'inchi')

    name_e.text = component.name
    id_e.text = component.id
    formula_e.text = component.formula
    sys_name_e.text = next((x.value for x in component.descriptors
                            if x.type == 'SYSTEMATIC NAME' and x.program == 'ACDLabs'), '')
    s_smiles_e.text = next((x.value for x in component.descriptors
                            if x.type == 'SMILES_CANONICAL' and x.program == 'CACTVS'), '')
    n_smiles_e.text = next((x.value for x in component.descriptors
                            if x.type == 'SMILES' and x.program == 'CACTVS'), '')
    inchi_e.text = component.inchi

    return root


def to_xml_str(component, remove_hs=True, conf_type=ConformerType.Ideal):
    """Converts structure to the XML format. Presently just molecule
    metadata are serialized without any coordinates, which is in
    accordance with the content of the PDBeChem area.

    Args:
        component (pdbeccdutils.core.Component): Component to be
            exported.
        remove_hs (bool, optional): Defaults to True.
        conf_type (pdbeccdutils.core.ConformerType, optional):
            Defaults to ConformerType.Ideal.

    Returns:
        str: String representation of the component in CML format.
    """
    root = to_xml_xml(component, remove_hs, conf_type)

    xml = ET.tostring(root, encoding='utf-8', method='xml')
    pretty = minidom.parseString(xml)

    return pretty.toprettyxml(indent="  ")


def to_pdb_ccd_cif_file(path, component, remove_hs=True):
    """Converts structure to the PDB CIF format. Both model and ideal
    coordinates are stored. In case ideal coordinates are missing, rdkit
    attempts to generate 3D coordinates of the conformer.

    Args:
        path (str): Path to save cif file.
        component (pdbeccdutils.core.Component): Component to be
            exported.
        remove_hs (bool, optional): Defaults to True.
    """
    if type(component.ccd_cif_dict) is not dict:
        component.ccd_cif_dict = _to_pdb_ccd_cif_dict(component)

    cif_copy = copy.deepcopy(component.ccd_cif_dict)

    if remove_hs:
        h_indices = [i for i, x in enumerate(cif_copy['_chem_comp_atom']['type_symbol']) if x == "H"]
        h_names = [cif_copy['_chem_comp_atom']['atom_id'][i] for i in h_indices]

        hb_indices = []
        for key in ('atom_id_1', 'atom_id_2'):
            indices = [i for i, k in enumerate(cif_copy['_chem_comp_bond'][key]) if k in h_names]
            hb_indices += indices

        hb_indices = list(set(hb_indices))

        # scrap hydrogen atoms
        for key in cif_copy['_chem_comp_atom']:
            cif_copy['_chem_comp_atom'][key] = (
                [k for i, k in enumerate(cif_copy['_chem_comp_atom'][key]) if i not in h_indices])

        # scrap bonds to hydrogen atoms
        for key in cif_copy['_chem_comp_bond']:
            cif_copy['_chem_comp_bond'][key] = (
                [k for i, k in enumerate(cif_copy['_chem_comp_bond'][key]) if i not in hb_indices])

    cfd = mmcif.CifFileWriter(path)
    cfd.write(cif_copy)


def to_cml_str(component, remove_hs=True, conf_type=ConformerType.Ideal):
    """Converts structure to the EBI representation of the molecule in
    CML format: http://cml.sourceforge.net/schema/cmlCore.xsd

    Args:
        component (pdbeccdutils.core.Component): Component to be
            exported.
        remove_hs (bool, optional): Defaults to True.
        conf_type (pdbeccdutils.core.ConformerType, optional):
            Defaults to ConformerType.Ideal.

    Returns:
        str: String representation of the component in CML format.
    """
    (mol_to_save, conf_id, conf_type) = _prepate_structure(component, remove_hs, conf_type)

    root = ET.Element('cml')
    root.set('xsi:schemaLocation', 'http://cml.sourceforge.net/schema/cmlCore.xsd')
    root.set('xmlns:xsi', 'http://www.w3.org/2001/XMLSchema-instance')
    root.set('dictRef', 'ebiMolecule:ebiMoleculeDict.cml')
    root.set('ebiMolecule', 'http://www.ebi.ac.uk/felics/molecule')

    f_charge = sum([l.GetFormalCharge() for l in mol_to_save.GetAtoms()])
    mol = ET.SubElement(root, 'molecule', {'id': component.id, 'formalCharge': str(f_charge)})

    id_inchi = ET.SubElement(mol, 'identifier', {'dictRef': 'ebiMolecule:inchi'})
    id_inchi.text = component.inchi
    id_systematic = ET.SubElement(mol, 'identifier', {'dictRef': 'ebiMolecule:systematicName'})
    id_systematic.text = component.name
    id_formula1 = ET.SubElement(mol, 'formula', {'dictRef': 'ebiMolecule:stereoSmiles'})
    id_formula2 = ET.SubElement(mol, 'formula', {'dictRef': 'ebiMolecule:nonStereoSmiles'})
    id_formula1.text = next((x.value for x in component.descriptors
                             if x.type == 'SMILES_CANONICAL' and x.program == 'CACTVS'), '')
    id_formula2.text = next((x.value for x in component.descriptors
                             if x.type == 'SMILES' and x.program == 'CACTVS'), '')

    atom_array = ET.SubElement(mol, 'atomArray')
    conformer = mol_to_save.GetConformer(id=conf_id)

    for atom in mol_to_save.GetAtoms():
        element = atom.GetSymbol()
        a_name = _get_atom_name(atom)
        coords = conformer.GetAtomPosition(atom.GetIdx())

        a_entry = ET.SubElement(atom_array, 'atom', {'id': a_name, 'elementType': element})
        a_entry.set('x3', str(coords.x))
        a_entry.set('y3', str(coords.y))
        a_entry.set('z3', str(coords.z))

    bond_array = ET.SubElement(mol, 'bondArray')
    for bond in mol_to_save.GetBonds():
        atom_1 = _get_atom_name(bond.GetBeginAtom())
        atom_2 = _get_atom_name(bond.GetEndAtom())
        bond_order = _get_cml_bond_type(bond.GetBondType())

        bond_entry = ET.SubElement(bond_array, 'bond')
        bond_entry.set('atomsRefs2', atom_1 + ' ' + atom_2)
        bond_entry.set('order', str(bond_order))

    cml = ET.tostring(root, encoding='utf-8', method='xml')
    pretty = minidom.parseString(cml)

    return pretty.toprettyxml(indent="  ")


def _prepate_structure(component, remove_hs, conf_type):
    """Prepare structure for export based on parameters. If deemed
    conformation is missing, an exception is thrown.
    TODO: handling AllConformers for other than PDB, SDF formats.

    Args:
        component (pdbeccdutils.core.Component): Component to be
            exported.
        remove_hs (bool, optional): Defaults to True.
        conf_type (pdbeccdutils.core.ConformerType, optional):
            Defaults to ConformerType.Ideal.

    Returns:
        tuple(rdkit.Mol,int,ConformerType): mol along with properties
        to be exported.
    """
    conf_id = component.conformers_mapping[conf_type]
    mol_to_save = component._2dmol if conf_type == ConformerType.Depiction else component.mol

    if remove_hs:
        mol_to_save = component.mol_no_h

    return (mol_to_save, conf_id, conf_type)


def _to_pdb_ccd_cif_dict(component):
    """Export component to the PDB CCD CIF file.

    Args:
        component (pdbeccdutils.core.Component): Component to be
            exported.

    Returns:
        (dict of str: str)): dictionary representation of the component
            to be serialized as a cif file.
    """

    cif_dict = OrderedDict([])

    _write_pdb_ccd_cif_info(cif_dict, component)
    _write_pdb_ccd_cif_atoms(cif_dict, component)
    _write_pdb_ccd_cif_bonds(cif_dict, component)
    _write_pdb_ccd_cif_descriptor(cif_dict, component)

    return cif_dict


def _write_pdb_ccd_cif_info(cif_dict, component):
    """Writes _chem_comp namespace with general information about the
    component

    Args:
        cif_dict (dict of str: str): cif representation of the molecule
            in a dictionary.
        component (pdbeccdutils.core.Component): Component to be
            exported.
    """

    calc_formula = Chem.rdMolDescriptors.CalcMolFormula(component.mol)
    calc_weight = Chem.rdMolDescriptors.CalcExactMolWt(component.mol)
    mod_date = "{}-{:02d}-{:02d}".format(component.modified_date.year, component.modified_date.month, component.modified_date.day)

    cif_dict['pdbeccdutils'] = OrderedDict([])

    cif_dict['pdbeccdutils']['rdkit_version'] = rdkit.__version__
    cif_dict['pdbeccdutils']['core_version'] = pdbeccdutils.__version__

    label = '_chem_comp'
    cif_dict[label] = OrderedDict([])
    cif_dict[label]['id'] = component.id
    cif_dict[label]['type'] = 'NON-POLYMER'
    cif_dict[label]['pdbx_type'] = 'HETAIN'
    cif_dict[label]['formula'] = component.formula if len(component.formula) > 0 else calc_formula
    cif_dict[label]['formula_weight'] = '{:.3f}'.format(calc_weight)
    cif_dict[label]['three_letter_code'] = component.id
    cif_dict[label]['pdbx_type'] = 'HETAIN'
    cif_dict[label]['pdbx_modified_date'] = mod_date
    cif_dict[label]['pdbx_release_status'] = component.pdbx_release_status.name

    cif_dict['pdbeccdutils']['rdkit_version'] = rdkit.__version__
    cif_dict['pdbeccdutils']['core_version'] = pdbeccdutils.__version__


def _write_pdb_ccd_cif_atoms(cif_dict, component):
    """Writes the _chem_comp_atom namespace with atom details.
    Controlled dictionary:
    http://mmcif.wwpdb.org/dictionaries/mmcif_pdbx_v50.dic/Categories/chem_comp_atom.html

    Args:
        cif_dict (dict of str: str): cif representation of the molecule
            in a dictionary.
        component (pdbeccdutils.core.Component): Component to be
            exported.
    """

    if component.mol.GetNumAtoms() < 1:
        return

    # if we dont have ideal coordinates, we just try to compute them
    ideal_type = ConformerType.Ideal
    if component.has_degenerated_conformer(ConformerType.Ideal):
        ideal_type = ConformerType.Computed if component.compute_3d() else ConformerType.Ideal

    label = '_chem_comp_atom'

    cif_dict[label] = OrderedDict([])
    cif_dict[label]['comp_id'] = []
    cif_dict[label]['atom_id'] = []
    cif_dict[label]['alt_atom_id'] = []
    cif_dict[label]['type_symbol'] = []
    cif_dict[label]['charge'] = []
    cif_dict[label]['pdbx_align'] = []
    cif_dict[label]['pdbx_aromatic_flag'] = []
    cif_dict[label]['pdbx_leaving_atom_flag'] = []
    cif_dict[label]['pdbx_stereo_config'] = []
    cif_dict[label]['model_Cartn_x'] = []
    cif_dict[label]['model_Cartn_y'] = []
    cif_dict[label]['model_Cartn_z'] = []
    cif_dict[label]['pdbx_model_Cartn_x_ideal'] = []
    cif_dict[label]['pdbx_model_Cartn_y_ideal'] = []
    cif_dict[label]['pdbx_model_Cartn_z_ideal'] = []
    cif_dict[label]['pdbx_component_atom_id'] = []
    cif_dict[label]['pdbx_component_comp_id'] = []
    cif_dict[label]['pdbx_ordinal'] = []

    for atom in component.mol.GetAtoms():
        at_id = atom.GetIdx()

        model_atom = _get_atom_coord(component, at_id, ConformerType.Model)
        ideal_atom = _get_atom_coord(component, at_id, ideal_type)

        cif_dict[label]['comp_id'].append(component.id)
        cif_dict[label]['atom_id'].append(_get_atom_name(atom))
        cif_dict[label]['alt_atom_id'].append(_get_atom_name(atom))
        cif_dict[label]['type_symbol'].append(atom.GetSymbol())
        cif_dict[label]['charge'].append(str(atom.GetFormalCharge()))
        cif_dict[label]['pdbx_align'].append('?')
        cif_dict[label]['pdbx_aromatic_flag'].append('Y' if atom.GetIsAromatic() else 'N')
        cif_dict[label]['pdbx_leaving_atom_flag'].append('N')
        cif_dict[label]['pdbx_stereo_config'].append(_get_ccd_cif_chiral_type(atom))

        cif_dict[label]['model_Cartn_x'].append("{:.3f}".format(model_atom.x))
        cif_dict[label]['model_Cartn_y'].append("{:.3f}".format(model_atom.y))
        cif_dict[label]['model_Cartn_z'].append("{:.3f}".format(model_atom.z))

        cif_dict[label]['pdbx_model_Cartn_x_ideal'].append("{:.3f}".format(ideal_atom.x))
        cif_dict[label]['pdbx_model_Cartn_y_ideal'].append("{:.3f}".format(ideal_atom.y))
        cif_dict[label]['pdbx_model_Cartn_z_ideal'].append("{:.3f}".format(ideal_atom.z))

        cif_dict[label]['pdbx_component_atom_id'].append(_get_atom_name(atom))
        cif_dict[label]['pdbx_component_comp_id'].append(component.id)
        cif_dict[label]['pdbx_ordinal'].append(str(atom.GetIdx() + 1))


def _write_pdb_ccd_cif_bonds(cif_dict, component):
    """Writes the _chem_comp_bond namespace with atom details.
    Controlled dictionary:
    http://mmcif.wwpdb.org/dictionaries/mmcif_pdbx_v50.dic/Categories/chem_comp_bond.html

    Args:
        cif_dict (dict of str: str): cif representation of the molecule
            in a dictionary.
        component (pdbeccdutils.core.Component): Component to be
            exported.
    """
    if component.mol.GetNumBonds() < 1:
        return

    label = '_chem_comp_bond'
    cif_dict[label] = OrderedDict([])
    cif_dict[label]['comp_id'] = []
    cif_dict[label]['atom_id_1'] = []
    cif_dict[label]['atom_id_2'] = []
    cif_dict[label]['value_order'] = []
    cif_dict[label]['pdbx_aromatic_flag'] = []
    cif_dict[label]['pdbx_stereo_config'] = []
    cif_dict[label]['pdbx_ordinal'] = []

    for b in component.mol.GetBonds():
        atom_a = b.GetBeginAtom()
        atom_b = b.GetEndAtom()

        cif_dict[label]['comp_id'].append(component.id)
        cif_dict[label]['atom_id_1'].append(_get_atom_name(atom_a))
        cif_dict[label]['atom_id_2'].append(_get_atom_name(atom_b))
        cif_dict[label]['value_order'].append(_get_ccd_cif_bond_type(b))
        cif_dict[label]['pdbx_aromatic_flag'].append('Y' if b.GetIsAromatic() else 'N')
        cif_dict[label]['pdbx_stereo_config'].append(_get_ccd_cif_bond_stereo(b))
        cif_dict[label]['pdbx_ordinal'].append(str(b.GetIdx() + 1))


def _write_pdb_ccd_cif_descriptor(cif_dict, component):
    """Writes the _pdbx_chem_comp_descriptor namespace with details.

    Args:
        cif_dict (dict of str: str): cif representation of the molecule
            in a dictionary.
        component (pdbeccdutils.core.Component): Component to be
            exported.
    """

    if len(component.descriptors) < 1:
        return

    label = '_pdbx_chem_comp_descriptor'
    cif_dict[label] = OrderedDict([])
    cif_dict[label]['comp_id'] = []
    cif_dict[label]['type'] = []
    cif_dict[label]['program'] = []
    cif_dict[label]['descriptor'] = []

    for entry in component.descriptors:
        cif_dict[label]['comp_id'].append(component.id)
        cif_dict[label]['type'].append(entry.type)
        cif_dict[label]['program'].append(entry.program)
        cif_dict[label]['descriptor'].append(entry.value)


def _get_atom_name(atom):
    """Gets atom name. If not set ElementSymbol + Id is used.

    Args:
        atom (rdkit.Chem.rdchem.Atom): rdkit atom.

    Returns:
        str: Name of the atom.
    """
    return atom.GetProp('name') if atom.HasProp('name') else atom.GetSymbol() + str(atom.GetIdx())


def _get_cml_bond_type(bond_order):
    """Translate bond type from rdkit to CML language.

    Args:
        bond_order (rdkit.Chem.rdchem.BondType): rdkit bond type

    Returns:
        str: bond type in the CML language.
    """
    if bond_order == Chem.rdchem.BondType.SINGLE:
        return '1'
    elif bond_order == Chem.rdchem.BondType.DOUBLE:
        return '2'
    elif bond_order == Chem.rdchem.BondType.TRIPLE:
        return '3'
    elif bond_order == Chem.rdchem.BondType.AROMATIC:
        return 'A'
    else:
        return str(bond_order)


def _get_ccd_cif_bond_stereo(bond):
    """Get bond stereochemistry information to be used in CCD CIF file.
    Controlled dictionary from: http://mmcif.wwpdb.org/dictionaries/mmcif_pdbx_v50.dic/Items/_chem_comp_bond.pdbx_stereo_config.html

    Args:
        bond (rdkit.Chem.rdchem.Bond): Molecular bond.

    Returns:
        str: bond stereochemistry information for the field
            '_chem_comp_bond.pdbx_stereo_config'.
    """

    stereo = bond.GetStereo()

    if bond.GetStereo() in (Chem.rdchem.BondStereo.STEREOE, Chem.rdchem.BondStereo.STEREOCIS):
        return 'E'
    elif bond.GetStereo() in (Chem.rdchem.BondStereo.STEREOZ, Chem.rdchem.BondStereo.STEREOTRANS):
        return 'Z'
    else:
        return 'N'


def _get_ccd_cif_bond_type(bond):
    """Translate bond type from rdkit to the CCD CIF format.
    Controlled dictionary from:
    http://mmcif.wwpdb.org/dictionaries/mmcif_pdbx_v50.dic/Items/_chem_comp_bond.value_order.html

    Args:
        bond_order (rdkit.Chem.rdchem.Bond): rdkit molecular bond.

    Returns:
        str: bond type for the '_chem_comp_bond.value_order' field. If
            none of the rdkit bond types can be matched SING is returned.
    """
    bond_order = bond.GetBondType()

    if bond_order == Chem.rdchem.BondType.SINGLE:
        return 'SING'
    elif bond_order == Chem.rdchem.BondType.DOUBLE:
        return 'DOUB'
    elif bond_order == Chem.rdchem.BondType.TRIPLE:
        return 'TRIP'
    elif bond_order == Chem.rdchem.BondType.AROMATIC:
        return 'AROM'
    elif bond_order == Chem.rdchem.BondType.QUADRUPLE:
        return 'QUAD'
    else:
        return 'SING'


def _get_ccd_cif_chiral_type(atom):
    """Translate atom chiral from rdkit to the CCD CIF format.
    Controlled dictionary from:
    http://mmcif.wwpdb.org/dictionaries/mmcif_pdbx_v50.dic/Items/_chem_comp_atom.pdbx_stereo_config.html

    Args:
        atom (rdkit.Chem.rdchem.Atom): rdkit atom

    Returns:
        str: chiral type for the '_chem_comp_atom.pdbx_stereo_config'
            field. If none of the rdkit atom types can be matched 'N' is
            returned.
    """
    chiral_type = atom.GetChiralTag()

    if chiral_type == Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CW:
        return 'R'
    elif chiral_type == Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CCW:
        return 'S'
    else:
        return 'N'


def _get_atom_coord(component, at_id, conformer_type):
    """Retrieve 3D coordinates for the particular atom in a given
    conformer.

    Args:
        component (pdbeccdutils.core.Component): component to be processed.
        at_id (int): atom id to be retrieved.
        conformer_type (pdbeccdutils.core.ConformerType): conformer type

    Returns:
        rdkit.Geometry.rdGeometry.Point3D: 3D coordinates of the atom.
    """
    conf_id = component.conformers_mapping[conformer_type]

    return component.mol.GetConformer(conf_id).GetAtomPosition(at_id)

# region fallbacks


def _to_sdf_str_fallback(mol, id, mappings):
    """Fallback method to generate SDF file in case the default one in
    RDKit fails.

    Args:
        mol (rdkit.Chem.rdchem.Mol): rdkit mol to be exported
        id (str): component id.
        mappings (dict of str): Deemed mappings to be exported

    Returns:
        list of str: SDF representation of the component
    """
    content = []

    for k, v in mappings.items():
        try:
            rdkit_conformer = mol.GetConformer(v)
        except ValueError:
            continue

        atom_count = mol.GetNumAtoms()
        bond_count = mol.GetNumBonds()

        content.append('{} - {} conformer\n    RDKit   3D\n'.format(id, k))
        content.append('{:>3}{:3}  0  0  0  0  0  0  0  0999 V2000'.format(atom_count, bond_count))

        for i in range(0, atom_count):
            content.append('{:>10.4f}{:>10.4f}{:>10.4f} {:<3} 0{:>3}'.format(
                rdkit_conformer.GetAtomPosition(i).x,
                rdkit_conformer.GetAtomPosition(i).y,
                rdkit_conformer.GetAtomPosition(i).z,
                mol.GetAtomWithIdx(i).GetSymbol(),
                _charge_to_sdf(mol.GetAtomWithIdx(i).GetFormalCharge()),
            ))

        for i in range(0, bond_count):
            bond = mol.GetBondWithIdx(i)
            content.append('{:>3}{:>3}{:>3}{:>3}  0  0  0'.format(
                bond.GetBeginAtom().GetIdx() + 1,
                bond.GetEndAtom().GetIdx() + 1,
                _bond_type_to_sdf(bond),
                _bond_stereo_to_sdf(bond)
            ))
        content.append('M  END')
        content.append('$$$$')
    return content


def _to_pdb_str_fallback(mol, conf_id, info):
    """Fallback method to generate PDB file in case the default one in
    RDKit fails.

    Args:
        mol (rdkit.Chem.rdchem.Mol): Molecule to be writter.
        conf_id (int): conformer id to be writen.
        info (rdkit.Chem.rdchem.AtomPDBResidueInfo): atom metadata.

    Returns:
        str: String representation the component in the PDB format.
    """
    conformer_ids = []
    content = []

    if conf_id == -1:
        conformer_ids = [c.GetId() for c in mol.GetConformers()]
    else:
        conformer_ids = [conf_id]

    for m in conformer_ids:
        rdkit_conformer = mol.GetConformer(m)
        content.append('MODEL {:>4}'.format(m))

        for i in range(0, mol.GetNumAtoms()):
            atom = mol.GetAtomWithIdx(i)
            s = '{:<6}{:>5} {:<4} {:<3} {}{:>4}{}   {:>8.3f}{:>8.3f}{:>8.3f}{:>6.2f}{:>6.2f}          {:>2}{:>2}'\
                .format('HETATM' if info.GetIsHeteroAtom() else 'ATOM',
                        i + 1,
                        _get_atom_name(atom),
                        info.GetResidueName(),
                        info.GetChainId(),
                        info.GetResidueNumber(),
                        ' ',
                        rdkit_conformer.GetAtomPosition(i).x,
                        rdkit_conformer.GetAtomPosition(i).y,
                        rdkit_conformer.GetAtomPosition(i).z,
                        info.GetOccupancy(),
                        info.GetTempFactor(),
                        atom.GetSymbol(),
                        atom.GetFormalCharge())
            content.append(s)

        for i in range(0, mol.GetNumAtoms()):
            pivot = mol.GetAtomWithIdx(i)
            s = 'CONECT{:>5}'.format(i + 1)

            for b in pivot.GetBonds():
                end_atom = b.GetOtherAtomIdx(i)

                if end_atom < i:
                    continue

                for t in range(0, math.floor(b.GetBondTypeAsDouble())):
                    s += '{:>5}'.format(end_atom + 1)

            if len(s) > 11:
                content.append(s)

        content.append('ENDMDL')

    return "\n".join(content)


def _charge_to_sdf(charge):
    """Translate RDkit charge to the SDF language.

    Args:
        charge (int): Numerical atom charge.

    Returns:
        str: Str representation of a charge in the sdf language
    """
    if charge == -3:
        return "7"
    elif charge == -2:
        return "6"
    elif charge == -1:
        return "5"
    elif charge == 0:
        return "0"
    elif charge == 1:
        return "+1"
    elif charge == 2:
        return "+2"
    elif charge == 3:
        return "+4"
    else:
        return "0"


def _bond_stereo_to_sdf(bond):
    """Translate bond stereo information to the sdf language. Needs to
    be checked.

    Args:
        bond (rdkit.Chem.rdchem.Bond): bond to be processed.

    Returns:
        str: bond type in sdf language.
    """
    stereo = bond.GetStereo()

    if stereo == Chem.rdchem.BondStereo.STEREONONE:
        return '0'
    elif stereo == Chem.rdchem.BondStereo.STEREOANY:
        return '4'
    elif stereo in (Chem.rdchem.BondStereo.STEREOCIS, Chem.rdchem.BondStereo.STEREOTRANS):
        return '3'
    else:
        return '0'


def _bond_type_to_sdf(bond):
    """Get bond type in sdf language. Based on:
    http://www.nonlinear.com/progenesis/sdf-studio/v0.9/faq/sdf-file-format-guidance.aspx

    Args:
        bond_type (rdkit.Chem.rdchem.Bond): Bond to be processed

    Returns:
        str: String representation of the bond type.
    """
    bond_order = bond.GetBondType()

    if bond_order == Chem.rdchem.BondType.SINGLE:
        return '1'
    elif bond_order == Chem.rdchem.BondType.DOUBLE:
        return '2'
    elif bond_order == Chem.rdchem.BondType.TRIPLE:
        return '3'
    elif bond_order == Chem.rdchem.BondType.AROMATIC:
        return 'A'
    else:
        return '0'


# endregion fallbacks
