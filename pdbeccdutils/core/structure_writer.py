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

import xml.etree.ElementTree as ET
from xml.dom import minidom

from rdkit import Chem, rdBase

from pdbeccdutils.core import Component, ConformerType


def write_molecule(path, component, remove_hs=True, conf_type=ConformerType.Ideal):
    """Export molecule in a specified format.

    Args:
        path (str): Path to the file. Extension determines format to be
            used.
        component (pdbeccdutils.core.Component): Component to be exported
        remove_hs (bool, optional): Defaults to True.
        conf_type ([type], optional): Defaults to ConformerType.Ideal.

    Raises:
        NotImplementedError: For format exports under development
        ValueError: For unsupported format
    """
    extension = path.split('.')[-1].lower()
    str_representation = ''

    if extension == 'sdf':
        str_representation = to_sdf_str(component, remove_hs, conf_type)
    elif extension == 'pdb':
        str_representation = to_pdb_str(component, remove_hs, conf_type)
    elif extension in ('mmcif', 'cif'):
        str_representation = to_mmcif_str(component, remove_hs, conf_type)
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
    """Converts structure to the SDF format.

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

    pdb_title = 'TITLE     {} coordinates'.format(conf_type.name)
    pdb_title += ' for PDB CCD {}\n'.format(component.id)
    pdb_title += 'AUTHOR    ccd_utils using RDKit {}\n'.format(rdBase.rdkitVersion)
    pdb_string = pdb_title + Chem.MolToPDBBlock(mol_to_save, conf_id)

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

    return Chem.MolToMolBlock(mol_to_save, confId=conf_id)


def to_xyz_str(component, remove_hs=True, conf_type=ConformerType.Ideal):
    """Converts structure to the XYZ format.

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


def to_mmcif_str(component, remove_hs=True, conf_type=ConformerType.Ideal):
    raise NotImplementedError()


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
        bond_order = bond.GetBondType()

        bond_entry = ET.SubElement(bond_array, 'bond')
        bond_entry.set('atomsRefs2', atom_1 + ' ' + atom_2)
        bond_entry.set('order', str(bond_order))

    cml = ET.tostring(root, encoding='utf-8', method='xml')
    pretty = minidom.parseString(cml)

    return pretty.toprettyxml(indent="  ")


def _prepate_structure(component, remove_hs, conf_type):
    """Prepare structure for export based on parameters. If deemed
    conformation is missing, it is computed.
    TODO: handling AllConformers, 2D depiction

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
    conf_id = -1
    if component.has_degenerated_conformer(conf_type):
        if component.compute_3d():
            conf_id = component.conformers_mapping[ConformerType.Computed]
            conf_type = ConformerType.Computed
    else:
        conf_id = component.conformers_mapping[conf_type]

    mol_to_save = component._2dmol if conf_type == ConformerType.Depiction else component.mol

    if remove_hs:
        mol_to_save = Chem.RemoveHs(mol_to_save, sanitize=False)
        Chem.SanitizeMol(mol_to_save, catchErrors=True)

    return (mol_to_save, conf_id, conf_type)


def _get_atom_name(atom):
    """Gets atom name. If not set ElementSymbol + Id is used.

    Args:
        atom (rdkit.rdchem.Atom): rdkit Atom.

    Returns:
        str: Name of the atom.
    """
    return atom.GetProp('name') if atom.HasProp('name') else atom.GetSymbol() + str(atom.GetIdx())


def _get_bond_type(bond_order):
    """Translate bond type from rdkit to CML language.

    Args:
        bond_order (rdkit.rdchem.BondType): rdkit bond type

    Returns:
        str: bond type in CML language.
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
