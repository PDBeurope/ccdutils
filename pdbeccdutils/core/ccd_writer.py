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

"""Structure writing module. Presently the following formats are supported:
    SDF, CIF, PDB, JSON, XYZ, XML, CML.

Raises:
    CCDUtilsError: If deemed format is not supported or an unrecoverable
        error occurres.
"""
import json
import logging
import math
from pathlib import Path
from typing import List
from xml.dom import minidom
from xml.etree import ElementTree as ET
from xml.etree.ElementTree import Element, SubElement

import pdbeccdutils
import rdkit
import gemmi
from gemmi import cif
from pdbeccdutils.helpers import cif_tools
from pdbeccdutils.core.component import Component
from pdbeccdutils.core.exceptions import CCDUtilsError
from pdbeccdutils.core.models import ConformerType


def write_molecule(
    path,
    component: Component,
    remove_hs: bool = True,
    alt_names: bool = False,
    conf_type: ConformerType = ConformerType.Ideal,
):
    """Export molecule in a specified format. Presently supported formats
    are: PDB CCD CIF (*.cif); Mol file (*.sdf); Chemical Markup language
    (*.cml); PDB file (*.pdb); XYZ file (*.xyz); XML (*.xml).
    ConformerType.AllConformers is presently supported only for PDB.

    Args:
        path (str|Path): Path to the file. Suffix determines format to be
            used.
        component (Component): Component to be exported
        remove_hs (bool, optional): Defaults to True. Whether or not
            hydrogens should be removed.
        alt_names (bool, optional): Defaults to False. Whether or not
            alternate names should be exported.
        conf_type (ConformerType, optional):
            Defaults to ConformerType.Ideal. Conformer type to be
            exported.

    Raises:
        CCDUtilsError: For unsupported format
    """
    path = str(path) if isinstance(path, Path) else path

    extension = path.split(".")[-1].lower()
    str_representation = ""

    if extension in ("sdf", "mol"):
        str_representation = to_sdf_str(component, remove_hs, conf_type)
    elif extension == "pdb":
        str_representation = to_pdb_str(component, remove_hs, alt_names, conf_type)
    elif extension in ("mmcif", "cif"):
        to_pdb_ccd_cif_file(path, component, remove_hs)
        return
    elif extension == "cml":
        str_representation = to_cml_str(component, remove_hs, conf_type)
    elif extension == "xml":
        str_representation = to_xml_str(component, remove_hs, conf_type)
    elif extension == "xyz":
        str_representation = to_xyz_str(component, remove_hs, conf_type)
    elif extension == "json":
        str_representation = json.dumps(
            to_json_dict(component, remove_hs, conf_type), sort_keys=True, indent=4
        )
    else:
        raise CCDUtilsError("Unsupported file format: {}".format(extension))

    with open(path, "w") as f:
        f.write(str_representation)


def to_pdb_str(
    component: Component,
    remove_hs: bool = True,
    alt_names: bool = False,
    conf_type: ConformerType = ConformerType.Ideal,
):
    """Converts structure to the PDB format.

    Args:
        Component (Component): Component to be exported.
        remove_hs (bool, optional): Defaults to True.
        alt_names (bool, optional): Defaults to False. Whether or not
            alternate atom names should be exported.
        conf_type (ConformerType, optional): Defaults to ConformerType.Ideal.

    Returns:
        str: String representation of the component in the PDB format.
    """
    (mol_to_save, conf_id, conf_type) = _prepate_structure(
        component, remove_hs, conf_type
    )

    info = rdkit.Chem.rdchem.AtomPDBResidueInfo()
    info.SetResidueName(f"{component.id:>3}")
    info.SetTempFactor(20.0)
    info.SetOccupancy(1.0)
    info.SetChainId("A")
    info.SetResidueNumber(1)
    info.SetIsHeteroAtom(True)

    for atom in mol_to_save.GetAtoms():
        flag = _get_alt_atom_name(atom) if alt_names else _get_atom_name(atom)
        atom_name = f"{flag:<4}"  # make sure it is 4 characters
        info.SetName(atom_name)
        atom.SetMonomerInfo(info)

    pdb_title = [
        f"HEADER    {conf_type.name} coordinates",
        f" for PDB-CCD {component.id}",
        f"COMPND    {component.id}",
        f"AUTHOR    pdbccdutils {pdbeccdutils.__version__}",
        f"AUTHOR    RDKit {rdkit.__version__}",
    ]

    try:
        pdb_body = rdkit.Chem.MolToPDBBlock(mol_to_save, conf_id)
    except Exception:
        pdb_body = _to_pdb_str_fallback(
            mol_to_save, component.id, conf_id, conf_type.name
        )

    return "\n".join(pdb_title + [pdb_body])


def to_sdf_str(
    component: Component,
    remove_hs: bool = True,
    conf_type: ConformerType = ConformerType.Ideal,
):
    """Converts structure to the SDF format.

    Args:
        component (Component): Component to be exported.
        remove_hs (bool, optional): Defaults to True.
        conf_type (ConformerType, optional): Defaults to ConformerType.Ideal.

    Raises:
        CCDUtilsError: In case the structure could not be exported.

    Returns:
        str: String representation of the component in the SDF format
    """
    (mol_to_save, _, conf_type) = _prepate_structure(component, remove_hs, conf_type)

    mol_block = []

    if conf_type == ConformerType.AllConformers:
        conformers = [ConformerType.Model, ConformerType.Ideal, ConformerType.Computed]
    else:
        conformers = [conf_type]
    try:
        for c in conformers:
            try:
                conf_id = -1
                if c != ConformerType.AllConformers:
                    conf_id = component.get_conformer(c).GetId()

                block = [
                    f"{component.id} - {c.name} conformer",
                    rdkit.Chem.MolToMolBlock(mol_to_save, confId=conf_id).strip(),
                    "$$$$\n",
                ]
                mol_block += block
            except ValueError as e:
                if str(e) == "Bad Conformer Id":
                    pass
                else:
                    raise CCDUtilsError(f"Error writing SDF file - {e}")
    except Exception:
        mol_block = _to_sdf_str_fallback(mol_to_save, component.id, conformers)

    return "\n".join(mol_block)


def to_xyz_str(component, remove_hs=True, conf_type=ConformerType.Ideal):
    """Converts structure to the XYZ format. Does not yet support
    ConformerType.AllConformers.

    Args:
        component (Component): Component to be exported.
        remove_hs (bool, optional): Defaults to True.
        conf_type (ConformerType, optional): Defaults to ConformerType.Ideal.

    Returns:
        str: String representation of the component in the XYZ format
    """
    (mol_to_save, conf_id, conf_type) = _prepate_structure(
        component, remove_hs, conf_type
    )
    conformer = mol_to_save.GetConformer(id=conf_id)

    result = list()
    result.append(str(mol_to_save.GetNumAtoms()))
    result.append(component.id)

    for atom in mol_to_save.GetAtoms():
        coords = conformer.GetAtomPosition(atom.GetIdx())
        result.append(
            f"{atom.GetSymbol():<4}{coords.x: f}  {coords.y: f}  {coords.z: f}"
        )

    return "\n".join(result)


def to_xml_xml(component, remove_hs=True, conf_type=ConformerType.Ideal):
    """Converts structure to the XML format and returns its XML repr.

    Args:
        component (Component): Component to be exported.
        remove_hs (bool, optional): Defaults to True.
        conf_type (ConformerType, optional): Defaults to ConformerType.Ideal.

    Returns:
        xml.etree.ElementTree.Element: XML object
    """
    root = Element("chemComp")

    id_e = SubElement(root, "id")
    name_e = SubElement(root, "name")
    formula_e = SubElement(root, "formula")
    sys_name_e = SubElement(root, "systematicName")
    s_smiles_e = SubElement(root, "stereoSmiles")
    n_smiles_e = SubElement(root, "nonStereoSmiles")
    inchi_e = SubElement(root, "inchi")

    name_e.text = component.name
    id_e.text = component.id
    formula_e.text = component.formula
    sys_name_e.text = next(
        (
            x.value
            for x in component.descriptors
            if x.type == "SYSTEMATIC NAME" and x.program == "ACDLabs"
        ),
        "",
    )
    s_smiles_e.text = next(
        (
            x.value
            for x in component.descriptors
            if x.type == "SMILES_CANONICAL" and x.program == "CACTVS"
        ),
        "",
    )
    n_smiles_e.text = next(
        (
            x.value
            for x in component.descriptors
            if x.type == "SMILES" and x.program == "CACTVS"
        ),
        "",
    )
    inchi_e.text = component.inchi

    return root


def to_xml_str(component: Component, remove_hs=True, conf_type=ConformerType.Ideal):
    """Converts structure to the XML format. Presently just molecule
    metadata are serialized without any coordinates, which is in
    accordance with the content of the PDBeChem area.

    Args:
        component (Component): Component to be exported.
        remove_hs (bool, optional): Defaults to True.
        conf_type (ConformerType, optional): Defaults to ConformerType.Ideal.

    Returns:
        str: String representation of the component in CML format.
    """
    root = to_xml_xml(component, remove_hs, conf_type)

    xml = ET.tostring(root, encoding="utf-8", method="xml")
    xmls_tring = minidom.parseString(xml)

    return xmls_tring.toprettyxml(indent="  ")


def remove_hydrogens(cif_block_copy):
    cif_tools.preprocess_cif_category(cif_block_copy, "_chem_comp_atom.")
    cif_tools.preprocess_cif_category(cif_block_copy, "_chem_comp_bond.")

    # scrap hydrogen atoms
    h_names: List[str] = []
    atom_table = cif_block_copy.find("_chem_comp_atom.", ["type_symbol", "atom_id"])
    for i in range(len(atom_table) - 1, -1, -1):
        if atom_table[i][0] == "H":
            h_names.append(atom_table[i][1])
            del atom_table[i]

    # scrap bonds to hydrogen atoms
    if "_chem_comp_bond." not in cif_block_copy.get_mmcif_category_names():
        return

    bond_table = cif_block_copy.find("_chem_comp_bond.", ["atom_id_1", "atom_id_2"])
    for j in range(len(bond_table) - 1, -1, -1):
        if (bond_table[j][0] in h_names) or (bond_table[j][1] in h_names):
            del bond_table[j]


def to_pdb_ccd_cif_file(path, component: Component, remove_hs=True):
    """Converts structure to the PDB CIF format. Both model and ideal
    coordinates are stored. In case ideal coordinates are missing, rdkit
    attempts to generate 3D coordinates of the conformer.

    Args:
        path (str): Path to save cif file.
        component (Component): Component to be exported.
        remove_hs (bool, optional): Defaults to True.
    """

    if not isinstance(component.ccd_cif_block, gemmi.cif.Block):
        component.ccd_cif_block = _to_pdb_ccd_cif_block(component)

    temp_doc = cif.Document()
    cif_block_copy = temp_doc.add_copied_block(component.ccd_cif_block)

    _add_sw_info_cif(cif_block_copy)
    _add_2d_depiction_cif(component, cif_block_copy)
    _add_fragments_and_scaffolds_cif(component, cif_block_copy)
    _add_rdkit_properties_cif(component, cif_block_copy)
    _add_unichem_mapping_cif(component, cif_block_copy)
    _add_rdkit_conformer_cif(component, cif_block_copy, remove_hs)

    if remove_hs:
        remove_hydrogens(cif_block_copy)

    temp_doc.write_file(path, cif.Style.Pdbx)


def to_cml_str(component: Component, remove_hs=True, conf_type=ConformerType.Ideal):
    """Converts structure to the EBI representation of the molecule in
    CML format: http://cml.sourceforge.net/schema/cmlCore.xsd

    Args:
        component (Component): Component to be exported.
        remove_hs (bool, optional): Defaults to True.
        conf_type (ConformerType, optional): Defaults to ConformerType.Ideal.

    Returns:
        str: String representation of the component in CML format.
    """
    (mol_to_save, conf_id, conf_type) = _prepate_structure(
        component, remove_hs, conf_type
    )

    root = ET.Element("cml")
    root.set("xsi:schemaLocation", "http://cml.sourceforge.net/schema/cmlCore.xsd")
    root.set("xmlns:xsi", "http://www.w3.org/2001/XMLSchema-instance")
    root.set("dictRef", "ebiMolecule:ebiMoleculeDict.cml")
    root.set("ebiMolecule", "http://www.ebi.ac.uk/felics/molecule")

    f_charge = sum([a.GetFormalCharge() for a in mol_to_save.GetAtoms()])
    mol = ET.SubElement(
        root, "molecule", {"id": component.id, "formalCharge": str(f_charge)}
    )

    id_inchi = ET.SubElement(mol, "identifier", {"dictRef": "ebiMolecule:inchi"})
    id_inchi.text = component.inchi
    id_systematic = ET.SubElement(
        mol, "identifier", {"dictRef": "ebiMolecule:systematicName"}
    )
    id_systematic.text = component.name
    id_formula1 = ET.SubElement(mol, "formula", {"dictRef": "ebiMolecule:stereoSmiles"})
    id_formula2 = ET.SubElement(
        mol, "formula", {"dictRef": "ebiMolecule:nonStereoSmiles"}
    )
    id_formula1.text = next(
        (
            x.value
            for x in component.descriptors
            if x.type == "SMILES_CANONICAL" and x.program == "CACTVS"
        ),
        "",
    )
    id_formula2.text = next(
        (
            x.value
            for x in component.descriptors
            if x.type == "SMILES" and x.program == "CACTVS"
        ),
        "",
    )

    atom_array = ET.SubElement(mol, "atomArray")
    conformer = mol_to_save.GetConformer(id=conf_id)

    for atom in mol_to_save.GetAtoms():
        element = atom.GetSymbol()
        a_name = _get_atom_name(atom)
        coords = conformer.GetAtomPosition(atom.GetIdx())

        a_entry = ET.SubElement(
            atom_array, "atom", {"id": a_name, "elementType": element}
        )
        a_entry.set("x3", str(coords.x))
        a_entry.set("y3", str(coords.y))
        a_entry.set("z3", str(coords.z))

    bond_array = ET.SubElement(mol, "bondArray")
    for bond in mol_to_save.GetBonds():
        atom_1 = _get_atom_name(bond.GetBeginAtom())
        atom_2 = _get_atom_name(bond.GetEndAtom())
        bond_order = _get_cml_bond_type(bond.GetBondType())

        bond_entry = ET.SubElement(bond_array, "bond")
        bond_entry.set("atomsRefs2", atom_1 + " " + atom_2)
        bond_entry.set("order", str(bond_order))

    cml = ET.tostring(root, encoding="utf-8", method="xml")
    pretty = minidom.parseString(cml)

    return pretty.toprettyxml(indent="  ")


def to_json_dict(component: Component, remove_hs=True, conf_type=ConformerType.Ideal):
    """Returns component information in dictionary suitable for json
    formating

    Args:
        component (Component): Component to be
            exported.
        remove_hs (bool, optional): Defaults to True.
        conf_type (ConformerType, optional):
            Defaults to ConformerType.Ideal.

    Raises:
        AttributeError: If all conformers are requested. This feature is
        not supported not is planned.

    Returns:
        :obj:`dict` of :obj:`str`: dictionary representation of the component
    """
    if conf_type == ConformerType.AllConformers:
        raise AttributeError("All conformer export is not supported for json export")

    (mol_to_save, conf_id, conf_type) = _prepate_structure(
        component, remove_hs, conf_type
    )
    result = {"atoms": [], "bonds": []}
    conformer = mol_to_save.GetConformer(conf_id)

    for atom in mol_to_save.GetAtoms():
        atom_dict = {}

        coords = conformer.GetAtomPosition(atom.GetIdx())

        atom_dict["id"] = atom.GetProp("name")
        atom_dict["element"] = atom.GetSymbol()
        atom_dict["charge"] = atom.GetFormalCharge()

        if conformer.Is3D():
            atom_dict["coords"] = {
                "X": round(coords.x, 4),
                "Y": round(coords.y, 4),
                "Z": round(coords.z, 4),
            }
        else:
            atom_dict["coords"] = {"X": round(coords.x, 4), "Y": round(coords.y, 4)}

        result["atoms"].append(atom_dict)

    for bond in mol_to_save.GetBonds():
        bond_dict = {}
        bond_dict["from"] = bond.GetBeginAtomIdx()
        bond_dict["to"] = bond.GetEndAtomIdx()
        bond_dict["type"] = bond.GetBondType().name

        result["bonds"].append(bond_dict)

    return {component.id: result}


def to_json_str(component: Component, remove_hs=True, conf_type=ConformerType.Ideal):
    """Converts structure into JSON representation. https://www.json.org/

    Args:
        component (Component): Component to be exported.
        remove_hs (bool, optional): Defaults to True.
        conf_type (ConformerType, optional): Defaults to ConformerType.Ideal.

    Returns:
        str: json representation of the component as a string.

    """

    temp = to_json_dict(component, remove_hs, conf_type)
    temp["rdkit_version"] = rdkit.__version__
    temp["pdbeccdutils_version"] = pdbeccdutils.__version__

    return json.dumps(temp)


def _prepate_structure(component, remove_hs, conf_type):
    """Prepare structure for export based on parameters. If deemed
    conformation is missing, an exception is thrown.
    TODO: handling AllConformers for other than PDB, SDF formats.

    Args:
        component (pdbeccdutils.core.Component): Component to be
            exported.
        remove_hs (bool, optional): Defaults to True.
        conf_type (ConformerType, optional):
            Defaults to ConformerType.Ideal.

    Returns:
        tuple(rdkit.Mol,int,ConformerType): mol along with properties
        to be exported.
    """
    conf_id = -1  # this is AllConformers options

    if conf_type == ConformerType.Depiction:
        conf_id = 0
    else:
        for c in component.mol.GetConformers():
            if c.GetProp("name") == conf_type.name:
                conf_id = c.GetId()
                break

    mol_to_save = (
        component.mol2D if conf_type == ConformerType.Depiction else component.mol
    )

    if remove_hs:
        mol_to_save = rdkit.Chem.RemoveHs(mol_to_save, sanitize=False)
        mol_to_save.UpdatePropertyCache(strict=False)

    return (mol_to_save, conf_id, conf_type)


def _to_pdb_ccd_cif_block(component):
    """Export component to the PDB CCD CIF file.

    Args:
        component (pdbeccdutils.core.Component): Component to be
            exported.

    Returns:
        (dict of str: str)): dictionary representation of the component
            to be serialized as a cif file.
    """

    doc = cif.Document()
    cif_block = doc.add_new_block(component.id)

    _write_pdb_ccd_cif_info(cif_block, component)
    _write_pdb_ccd_cif_atoms(cif_block, component)
    _write_pdb_ccd_cif_bonds(cif_block, component)
    _write_pdb_ccd_cif_descriptor(cif_block, component)

    return cif_block


def _write_pdb_ccd_cif_info(cif_block, component):
    """Writes _chem_comp namespace with general information about the
    component

    Args:
        cif_block (cif.Block): mmcif Block object from gemmi.
        component (Component): Component to be exported.
    """

    calc_formula = rdkit.Chem.rdMolDescriptors.CalcMolFormula(component.mol)
    calc_weight = rdkit.Chem.rdMolDescriptors.CalcExactMolWt(component.mol)
    date = component.modified_date
    mod_date = f"{date.year}-{date.month:02d}-{date.day:02d}"

    label = "_chem_comp."
    cif_block.set_pairs(
        label,
        {
            "id": component.id,
            "type": "NON-POLYMER",
            "formula": component.formula or calc_formula,
            "formula_weight": f"{calc_weight:.3f}",
            "three_letter_code": component.id,
            "pdbx_type": "HETAIN",
            "pdbx_modified_date": mod_date,
            "pdbx_release_status": component.pdbx_release_status.name,
        },
        raw=False,
    )


def _write_pdb_ccd_cif_atoms(cif_block, component):
    """Writes the _chem_comp_atom namespace with atom details.
    Controlled dictionary:
    http://mmcif.wwpdb.org/dictionaries/mmcif_pdbx_v50.dic/Categories/chem_comp_atom.html

    Args:
        cif_block (cif.Block): mmcif Block object from gemmi.
        component (Component): Component to be exported.
    """

    if component.mol.GetNumAtoms() < 1:
        return

    # if we dont have ideal coordinates, we just try to compute them
    ideal_type = ConformerType.Ideal
    if component.has_degenerated_conformer(ConformerType.Ideal):
        ideal_type = (
            ConformerType.Computed if component.compute_3d() else ConformerType.Ideal
        )

    label = "_chem_comp_atom."
    atom_fields = [
        "comp_id",
        "atom_id",
        "alt_atom_id",
        "type_symbol",
        "charge",
        "pdbx_align",
        "pdbx_aromatic_flag",
        "pdbx_leaving_atom_flag",
        "pdbx_stereo_config",
        "model_Cartn_x",
        "model_Cartn_y",
        "model_Cartn_z",
        "pdbx_model_Cartn_x_ideal",
        "pdbx_model_Cartn_y_ideal",
        "pdbx_model_Cartn_z_ideal",
        "pdbx_component_atom_id",
        "pdbx_component_comp_id",
        "pdbx_ordinal",
    ]
    atom_loop = cif_block.init_loop(label, atom_fields)

    for atom in component.mol.GetAtoms():
        at_id = atom.GetIdx()
        model_atom = _get_atom_coord(component, at_id, ConformerType.Model)
        ideal_atom = _get_atom_coord(component, at_id, ideal_type)

        new_row = [
            component.id,
            _get_atom_name(atom),
            _get_alt_atom_name(atom),
            atom.GetSymbol(),
            str(atom.GetFormalCharge()),
            None,
            "Y" if atom.GetIsAromatic() else "N",
            "N",
            _get_ccd_cif_chiral_type(atom),
            f"{model_atom.x:.3f}",
            f"{model_atom.y:.3f}",
            f"{model_atom.z:.3f}",
            f"{ideal_atom.x:.3f}",
            f"{ideal_atom.y:.3f}",
            f"{ideal_atom.z:.3f}",
            _get_atom_name(atom),
            component.id,
            str(atom.GetIdx() + 1),
        ]

        atom_loop.add_row(cif.quote_list(new_row))


def _write_pdb_ccd_cif_bonds(cif_block, component):
    """Writes the _chem_comp_bond namespace with atom details.
    Controlled dictionary:
    http://mmcif.wwpdb.org/dictionaries/mmcif_pdbx_v50.dic/Categories/chem_comp_bond.html

    Args:
        cif_block (cif.Block): mmcif Block object from gemmi.
        component (Component): Component to be exported.
    """
    if component.mol.GetNumBonds() < 1:
        return

    label = "_chem_comp_bond."
    bond_fields = [
        "comp_id",
        "atom_id_1",
        "atom_id_2",
        "value_order",
        "pdbx_aromatic_flag",
        "pdbx_stereo_config",
        "pdbx_ordinal",
    ]
    bond_loop = cif_block.init_loop(label, bond_fields)

    for b in component.mol.GetBonds():
        atom_a = b.GetBeginAtom()
        atom_b = b.GetEndAtom()

        new_row = [
            component.id,
            _get_atom_name(atom_a),
            _get_atom_name(atom_b),
            _get_ccd_cif_bond_type(b),
            "Y" if b.GetIsAromatic() else "N",
            _get_ccd_cif_bond_stereo(b),
            str(b.GetIdx() + 1),
        ]
        bond_loop.add_row(cif.quote_list(new_row))


def _write_pdb_ccd_cif_descriptor(cif_block, component):
    """Writes the _pdbx_chem_comp_descriptor namespace with details.

    Args:
        cif_block (cif.Block): mmcif Block object from gemmi.
        component (Component): Component to be exported.
    """

    if not component.descriptors:
        return

    label = "_pdbx_chem_comp_descriptor."

    descriptor_fields = ["comp_id", "type", "program", "program_version", "descriptor"]
    descriptor_loop = cif_block.init_loop(label, descriptor_fields)

    for entry in component.descriptors:
        new_row = [
            component.id,
            entry.type,
            entry.program,
            entry.program_version,
            entry.value,
        ]
        descriptor_loop.add_row(cif.quote_list(new_row))


def _get_atom_name(atom):
    """Gets atom name. If not set ElementSymbol + Id is used.

    Args:
        atom (rdkit.Chem.rdchem.Atom): rdkit atom.

    Returns:
        str: Name of the atom.
    """
    return (
        atom.GetProp("name")
        if atom.HasProp("name")
        else atom.GetSymbol() + str(atom.GetIdx())
    )


def _get_alt_atom_name(atom):
    """Gets alternate atom name. If not set _get_atom_name method is
    used.

    Args:
        atom (rdkit.Chem.rdchem.Atom): rdkit atom.

    Returns:
        str: Name of the atom.
    """
    return (
        atom.GetProp("alt_name") if atom.HasProp("alt_name") else _get_atom_name(atom)
    )


def _get_cml_bond_type(bond_order):
    """Translate bond type from rdkit to CML language.

    Args:
        bond_order (rdkit.Chem.rdchem.BondType): rdkit bond type

    Returns:
        str: bond type in the CML language.
    """
    if bond_order == rdkit.Chem.rdchem.BondType.SINGLE:
        return "1"

    if bond_order == rdkit.Chem.rdchem.BondType.DOUBLE:
        return "2"

    if bond_order == rdkit.Chem.rdchem.BondType.TRIPLE:
        return "3"

    if bond_order == rdkit.Chem.rdchem.BondType.AROMATIC:
        return "A"
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

    if stereo in (
        rdkit.Chem.rdchem.BondStereo.STEREOE,
        rdkit.Chem.rdchem.BondStereo.STEREOCIS,
    ):
        return "E"

    if stereo in (
        rdkit.Chem.rdchem.BondStereo.STEREOZ,
        rdkit.Chem.rdchem.BondStereo.STEREOTRANS,
    ):
        return "Z"

    return "N"


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

    if bond_order == rdkit.Chem.rdchem.BondType.SINGLE:
        return "SING"

    if bond_order == rdkit.Chem.rdchem.BondType.DOUBLE:
        return "DOUB"

    if bond_order == rdkit.Chem.rdchem.BondType.TRIPLE:
        return "TRIP"

    if bond_order == rdkit.Chem.rdchem.BondType.AROMATIC:
        return "AROM"

    if bond_order == rdkit.Chem.rdchem.BondType.QUADRUPLE:
        return "QUAD"

    return "SING"


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

    if chiral_type == rdkit.Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CW:
        return "R"

    if chiral_type == rdkit.Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CCW:
        return "S"

    return "N"


def _get_atom_coord(component, at_id, conformer_type):
    """Retrieve 3D coordinates for the particular atom in a given
    conformer.

    Args:
        component (Component): component to be processed.
        at_id (int): atom id to be retrieved.
        conformer_type (ConformerType): conformer type

    Returns:
        rdkit.Geometry.rdGeometry.Point3D: 3D coordinates of the atom.
    """
    conformer = None
    for c in component.mol.GetConformers():
        if c.GetProp("name") == conformer_type.name:
            conformer = c
            break

    return conformer.GetAtomPosition(at_id)


# region fallbacks


def _to_sdf_str_fallback(mol, ccd_id, conformers):
    """Fallback method to generate SDF file in case the default one in
    RDKit fails.

    Args:
        mol (rdkit.Chem.rdchem.Mol): rdkit mol to be exported
        id (str): component id.
        conformers (list of of ConformerType): List of conformer types
            to be exported.

    Returns:
        list of str: SDF representation of the component
    """
    content = []

    for c in conformers:
        rdkit_conformer = None

        if c == ConformerType.AllConformers:
            rdkit_conformer = mol.GetConformer()
        else:
            for conf in mol.GetConformers():
                if conf.GetProp("name") == c.name:
                    rdkit_conformer = conf
                    break

        atom_count = mol.GetNumAtoms()
        bond_count = mol.GetNumBonds()

        content += [
            f"{ccd_id} - {c.name} conformer",
            "    RDKit   3D",
            "\n" f"{atom_count:>3}{bond_count:3}  0  0  0  0  0  0  0  0999 V2000",
        ]

        for i in range(0, atom_count):
            pos = rdkit_conformer.GetAtomPosition(i)
            atom = mol.GetAtomWithIdx(i)
            sdf_charge = __charge_to_sdf(atom.GetFormalCharge())
            content.append(
                f"{pos.x:>10.4f}{pos.y:>10.4f}{pos.z:>10.4f} {atom.GetSymbol():<3} 0{sdf_charge:>3}"
            )

        for i in range(0, bond_count):
            bond = mol.GetBondWithIdx(i)
            content.append(
                "{:>3}{:>3}{:>3}{:>3}  0  0  0".format(
                    bond.GetBeginAtom().GetIdx() + 1,
                    bond.GetEndAtom().GetIdx() + 1,
                    __bond_type_to_sdf(bond),
                    __bond_stereo_to_sdf(bond),
                )
            )
        content.append("M  END")
        content.append("$$$$\n")
    return content


def _to_pdb_str_fallback(mol, component_id, conf_id, conf_name="Model"):
    """Fallback method to generate PDB file in case the default one in
    RDKit fails.

    Args:
        mol (rdkit.Chem.rdchem.Mol): Molecule to be written.
        component_id (str): Component id.
        conf_id (int): conformer id to be written.
        conf_name (str): conformer name to be written.

    Returns:
        str: String representation the component in the PDB format.
    """
    conformer_ids = []
    content = [
        f"HEADER    {conf_name} coordinates for PDB-CCD {component_id}",
        f"COMPND    {component_id}",
        f"AUTHOR    pdbccdutils {pdbeccdutils.__version__}",
        f"AUTHOR    RDKit {rdkit.__version__}",
    ]

    if conf_id == -1:
        conformer_ids = [c.GetId() for c in mol.GetConformers()]
    else:
        conformer_ids = [conf_id]

    for m in conformer_ids:
        rdkit_conformer = mol.GetConformer(m)

        for i in range(0, mol.GetNumAtoms()):
            atom = mol.GetAtomWithIdx(i)

            s = "{:<6}{:>5} {:<4} {:>3} {}{:>4}{}   {:>8.3f}{:>8.3f}{:>8.3f}{:>6.2f}{:>6.2f}          {:>2}{:>2}".format(
                "HETATM",
                i + 1,
                atom.GetProp("name"),
                component_id,
                "A",
                1,
                " ",
                rdkit_conformer.GetAtomPosition(i).x,
                rdkit_conformer.GetAtomPosition(i).y,
                rdkit_conformer.GetAtomPosition(i).z,
                1,
                20,
                atom.GetSymbol(),
                atom.GetFormalCharge(),
            )
            content.append(s)

        for i in range(0, mol.GetNumAtoms()):
            pivot = mol.GetAtomWithIdx(i)
            s = "CONECT{:>5}".format(i + 1)

            for b in pivot.GetBonds():
                end_atom = b.GetOtherAtomIdx(i)

                if end_atom < i:
                    continue

                for t in range(0, math.floor(b.GetBondTypeAsDouble())):
                    s += "{:>5}".format(end_atom + 1)

            if len(s) > 11:
                content.append(s)

        content.append("END")

    return "\n".join(content)


# endregion fallbacks


# region helpers
def __charge_to_sdf(charge):
    """Translate RDkit charge to the SDF language.

    Args:
        charge (int): Numerical atom charge.

    Returns:
        str: Str representation of a charge in the sdf language
    """
    if charge == -3:
        return "7"
    if charge == -2:
        return "6"
    if charge == -1:
        return "5"
    if charge == 0:
        return "0"
    if charge == 1:
        return "+1"
    if charge == 2:
        return "+2"
    if charge == 3:
        return "+4"

    return "0"


def __bond_stereo_to_sdf(bond):
    """Translate bond stereo information to the sdf language. Needs to
    be checked.

    Args:
        bond (rdkit.Chem.rdchem.Bond): bond to be processed.

    Returns:
        str: bond type in sdf language.
    """
    stereo = bond.GetStereo()

    if stereo == rdkit.Chem.rdchem.BondStereo.STEREONONE:
        return "0"
    if stereo == rdkit.Chem.rdchem.BondStereo.STEREOANY:
        return "4"
    if stereo in (
        rdkit.Chem.rdchem.BondStereo.STEREOCIS,
        rdkit.Chem.rdchem.BondStereo.STEREOTRANS,
    ):
        return "3"
    return "0"


def __bond_type_to_sdf(bond):
    """Get bond type in sdf language. Based on:
    http://www.nonlinear.com/progenesis/sdf-studio/v0.9/faq/sdf-file-format-guidance.aspx

    Args:
        bond_type (rdkit.Chem.rdchem.Bond): Bond to be processed

    Returns:
        str: String representation of the bond type.
    """
    bond_order = bond.GetBondType()

    if bond_order == rdkit.Chem.rdchem.BondType.SINGLE:
        return "1"
    if bond_order == rdkit.Chem.rdchem.BondType.DOUBLE:
        return "2"
    if bond_order == rdkit.Chem.rdchem.BondType.TRIPLE:
        return "3"
    if bond_order == rdkit.Chem.rdchem.BondType.AROMATIC:
        return "A"

    return "0"


# endregion helpers


# region cif-export addition


def _add_sw_info_cif(cif_block_copy):
    """Add information to the cif file about software versions used
    to generate this information

    Args:
        cif_block_copy (Block: gemmi.cif.Block): Block representation of the molecule
        from gemmi.
    """
    category = "_software."
    sw_fields = ["name", "version", "description"]
    sw_loop = cif_block_copy.init_loop(category, sw_fields)
    sw_loop.add_row(cif.quote_list(["rdkit", rdkit.__version__, "Core functionality."]))
    sw_loop.add_row(
        cif.quote_list(
            [
                "pdbeccdutils",
                pdbeccdutils.__version__,
                "Wrapper to provide 2D templates and molecular fragments.",
            ]
        )
    )


def _add_2d_depiction_cif(component, cif_block_copy):
    """Add 2D coordinates of the component depiction

    Args:
        component (Component): pdbeccdutils component.
        cif_block_copy (Block: gemmi.cif.Block): Block representation of the molecule
        from gemmi.
    """
    if component.mol2D is None:
        return

    __add_rdkit_2d_atoms_cif(component, cif_block_copy)
    __add_rdkit_2d_bonds_cif(component, cif_block_copy)


def _add_fragments_and_scaffolds_cif(component, cif_block_copy):
    """Add fragments and scaffolds information to the CIF export

    Args:
        component (Component): Component to be exported.
        cif_block_copy (Block: gemmi.cif.Block): Block representation of the molecule
        from gemmi.
    """
    substructure_category = "_pdbe_chem_comp_substructure."
    substructure_fields = [
        "comp_id",
        "substructure_name",
        "id",
        "substructure_type",
        "substructure_smiles",
        "substructure_inchis",
        "substructure_inchikeys",
    ]
    substructure_loop = cif_block_copy.init_loop(
        substructure_category, substructure_fields
    )

    for i, scaffold in enumerate(component.scaffolds):
        mol = rdkit.Chem.MolFromSmiles(scaffold.smiles)
        inchi = rdkit.Chem.MolToInchi(mol)
        inchikey = rdkit.Chem.MolToInchiKey(mol)
        new_row = [
            component.id,
            scaffold.name,
            f"S{i+1}",
            "scaffold",
            scaffold.smiles,
            inchi or None,
            inchikey or None,
        ]
        substructure_loop.add_row(cif.quote_list(new_row))

    for j, fragment in enumerate(component.fragments):
        mol = rdkit.Chem.MolFromSmiles(fragment.smiles)
        inchi = rdkit.Chem.MolToInchi(mol)
        inchikey = rdkit.Chem.MolToInchiKey(mol)
        new_row = [
            component.id,
            fragment.name,
            f"F{j+1}",
            "fragment",
            fragment.smiles,
            inchi or None,
            inchikey or None,
        ]
        substructure_loop.add_row(cif.quote_list(new_row))

    mapping_category = "_pdbe_chem_comp_substructure_mapping."
    mapping_fields = ["comp_id", "atom_id", "substructure_id", "substructure_ordinal"]
    mapping_loop = cif_block_copy.init_loop(mapping_category, mapping_fields)

    for i, scaffold in enumerate(component.scaffolds):
        for a, mapping in enumerate(scaffold.mappings):
            for atom_name in mapping:
                new_row = [component.id, atom_name, f"S{i+1}", f"{a+1}"]
                mapping_loop.add_row(cif.quote_list(new_row))

    for i, fragment in enumerate(component.fragments):
        for a, mapping in enumerate(fragment.mappings):
            for atom_name in mapping:
                new_row = [component.id, atom_name, f"F{i+1}", f"{a+1}"]
                mapping_loop.add_row(cif.quote_list(new_row))


def _add_rdkit_properties_cif(component, cif_block_copy):
    """Add properties calculated by the rdkit to the dictionary

    Args:
        component (Component): pdbeccdutils component.
        cif_block (Block: gemmi.cif.Block): Block representation of the molecule
        from gemmi.
    """
    category = "_pdbe_chem_comp_rdkit_properties."
    cif_block_copy.set_pairs(category, {"comp_id": component.id})

    for k, v in component.physchem_properties.items():
        cif_block_copy.set_pairs(
            category, {k: f"{v:.0f}" if v.is_integer() else f"{v:.3f}"}
        )


def __add_rdkit_2d_atoms_cif(component, cif_block_copy):
    """Add atom coordinates to the mmCIF file.

    Args:
        component (Component): pdbeccdutils component.
        cif_block_copy (Block: gemmi.cif.Block): Block representation of the molecule
        from gemmi.
    """
    category = "_pdbe_chem_comp_atom_depiction."
    atom_depiction_fields = [
        "comp_id",
        "atom_id",
        "element",
        "model_Cartn_x",
        "model_Cartn_y",
        "pdbx_ordinal",
    ]
    conformer = component.mol2D.GetConformer()
    atom_depiction_loop = cif_block_copy.init_loop(category, atom_depiction_fields)

    for i, atom in enumerate(component.mol2D.GetAtoms()):
        new_row = [
            component.id,
            atom.GetProp("name"),
            atom.GetSymbol(),
            f"{conformer.GetAtomPosition(i).x:.3f}",
            f"{conformer.GetAtomPosition(i).y:.3f}",
            str(i + 1),
        ]
        atom_depiction_loop.add_row(cif.quote_list(new_row))


def __add_rdkit_2d_bonds_cif(component, cif_block_copy):
    """Add bond information to the mmCIF file.

    Args:
        component (Component): pdbeccdutils component.
        cif_block (Block: gemmi.cif.Block): Block representation of the molecule
        from gemmi.
    """
    category = "_pdbe_chem_comp_bond_depiction."
    fields = [
        "comp_id",
        "atom_id_1",
        "atom_id_2",
        "value_order",
        "bond_dir",
        "pdbx_ordinal",
    ]
    bonds_depiction_loop = cif_block_copy.init_loop(category, fields)

    try:
        copy = rdkit.Chem.Draw.rdMolDraw2D.PrepareMolForDrawing(
            component.mol2D, wedgeBonds=True, kekulize=True, addChiralHs=True
        )
    except (RuntimeError, ValueError):
        try:
            copy = rdkit.Chem.Draw.rdMolDraw2D.PrepareMolForDrawing(
                component.mol2D, wedgeBonds=False, kekulize=True, addChiralHs=True
            )
        except (RuntimeError, ValueError):
            copy = rdkit.Chem.Draw.rdMolDraw2D.PrepareMolForDrawing(
                component.mol2D, wedgeBonds=False, kekulize=True, addChiralHs=False
            )

    bond_depiction_ordinal = 0
    for i, b in enumerate(copy.GetBonds()):
        if b.GetEndAtom().GetSymbol() != "H":
            bond_depiction_ordinal += 1
            new_row = [
                component.id,
                b.GetBeginAtom().GetProp("name"),
                b.GetEndAtom().GetProp("name"),
                b.GetBondType().name,
                b.GetBondDir().name,
                str(bond_depiction_ordinal),
            ]
            bonds_depiction_loop.add_row(cif.quote_list(new_row))


def _add_unichem_mapping_cif(component, cif_block_copy):
    """Add UniChem mapping to CCD CIF.

    Args:
        component (Component): pdbeccdutils component.
        cif_block_copy (Block: gemmi.cif.Block): Block representation of the molecule
        from gemmi.
    """

    category = "_pdbe_chem_comp_external_mappings."
    fields = ["comp_id", "source", "resource", "resource_id"]
    external_mapping_loop = cif_block_copy.init_loop(category, fields)
    for mapping in component.external_mappings:
        new_row = [component.id, "UniChem", mapping[0], mapping[1]]
        external_mapping_loop.add_row(cif.quote_list(new_row))


def _add_rdkit_conformer_cif(component, cif_block_copy, remove_hs):
    """Add 3D coordinates generated by RDKit.

    component (Component): pdbeccdutils component.
    cif_block (Block: gemmi.cif.Block): Block representation of the molecule
        from gemmi.
    remove_hs (boolean): Whether or not hydrogen atoms should be written
    """

    try:
        conformer = component.get_conformer(ConformerType.Computed)
    except ValueError:
        logging.warning(f"{component.id} | Computed conformer does not exist.")
        return  # no conformer nothing to write, we quit

    category = "_pdbe_chem_comp_rdkit_conformer."
    rdkit_conformer_fields = [
        "comp_id",
        "atom_id",
        "Cartn_x_rdkit",
        "Cartn_y_rdkit",
        "Cartn_z_rdkit",
        "rdkit_method",
        "rdkit_ordinal",
    ]
    rdkit_conformer_loop = cif_block_copy.init_loop(category, rdkit_conformer_fields)
    method = conformer.GetProp("coord_generation")

    # we need to get all the allowed incides
    if remove_hs:
        atom_indices = [
            i
            for i in range(0, component.mol.GetNumAtoms())
            if component.mol.GetAtomWithIdx(i).GetSymbol() != "H"
        ]
    else:
        atom_indices = list(range(0, component.mol.GetNumAtoms()))

    for i, atom_index in enumerate(atom_indices):
        new_row = [
            component.id,
            component.mol.GetAtomWithIdx(atom_index).GetProp("name"),
            f"{conformer.GetAtomPosition(atom_index).x:.3f}",
            f"{conformer.GetAtomPosition(atom_index).y:.3f}",
            f"{conformer.GetAtomPosition(atom_index).z:.3f}",
            method,
            str(i + 1),
        ]
        rdkit_conformer_loop.add_row(cif.quote_list(new_row))


# endregion
