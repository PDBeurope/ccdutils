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
import copy
import json
import math
import xml.etree.ElementTree as ET
from collections import OrderedDict
from typing import List
from xml.dom import minidom

import pdbecif.mmcif_io as mmcif
import rdkit

import pdbeccdutils
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
        path (str): Path to the file. Extension determines format to be
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
        component (Component): Component to be exported.
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
    info.SetResidueName(f'{component.id:>3}')
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
    (mol_to_save, conf_id, conf_type) = _prepate_structure(
        component, remove_hs, conf_type
    )

    mol_block = []
    mappings = {}
    if conf_type == ConformerType.AllConformers:
        conformers = [ConformerType.Model, ConformerType.Ideal, ConformerType.Computed]
    else:
        conformers = [conf_type]
    try:
        for conf in conformers:
            try:
                block = [
                    f"{component.id} - {conf.name} conformer",
                    rdkit.Chem.MolToMolBlock(
                        mol_to_save, confId=component.conformers_mapping[conf]
                    ).strip(),
                    "$$$$\n",
                ]
                mol_block += block
            except ValueError as e:
                if str(e) == "Bad Conformer Id":
                    pass
                else:
                    raise CCDUtilsError(f"Error writing SDF file - {e}")
    except Exception:
        mappings = {m.name: component.conformers_mapping[m] for m in conformers}
        mol_block = _to_sdf_str_fallback(mol_to_save, component.id, mappings)

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
    root = ET.Element("chemComp")

    id_e = ET.SubElement(root, "id")
    name_e = ET.SubElement(root, "name")
    formula_e = ET.SubElement(root, "formula")
    sys_name_e = ET.SubElement(root, "systematicName")
    s_smiles_e = ET.SubElement(root, "stereoSmiles")
    n_smiles_e = ET.SubElement(root, "nonStereoSmiles")
    inchi_e = ET.SubElement(root, "inchi")

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
    pretty = minidom.parseString(xml)

    return pretty.toprettyxml(indent="  ")


def to_pdb_ccd_cif_file(path, component: Component, remove_hs=True):
    """Converts structure to the PDB CIF format. Both model and ideal
    coordinates are stored. In case ideal coordinates are missing, rdkit
    attempts to generate 3D coordinates of the conformer.

    Args:
        path (str): Path to save cif file.
        component (Component): Component to be exported.
        remove_hs (bool, optional): Defaults to True.
    """
    if not isinstance(component.ccd_cif_dict, dict):
        component.ccd_cif_dict = _to_pdb_ccd_cif_dict(component)

    cif_copy = copy.deepcopy(component.ccd_cif_dict)

    _add_sw_info_cif(cif_copy)
    _add_2d_depiction_cif(component, cif_copy)
    _add_fragments_and_scaffolds_cif(component, cif_copy)
    _add_rdkit_properties_cif(component, cif_copy)
    _add_unichem_mapping_cif(component, cif_copy)
    _add_rdkit_conformer_cif(component, cif_copy, remove_hs)

    if remove_hs:
        h_indices: List[int] = [
            i
            for i, x in enumerate(cif_copy["_chem_comp_atom"]["type_symbol"])
            if x == "H"
        ]
        h_names: List[str] = [
            cif_copy["_chem_comp_atom"]["atom_id"][i] for i in h_indices
        ]

        hb_indices = []
        for key in ("atom_id_1", "atom_id_2"):
            indices = [
                i
                for i, k in enumerate(cif_copy["_chem_comp_bond"][key])
                if k in h_names
            ]
            hb_indices += indices

        hb_indices = list(set(hb_indices))

        # scrap hydrogen atoms
        for key in cif_copy["_chem_comp_atom"]:
            cif_copy["_chem_comp_atom"][key] = [
                k
                for i, k in enumerate(cif_copy["_chem_comp_atom"][key])
                if i not in h_indices
            ]

        # scrap bonds to hydrogen atoms
        for key in cif_copy["_chem_comp_bond"]:
            cif_copy["_chem_comp_bond"][key] = [
                k
                for i, k in enumerate(cif_copy["_chem_comp_bond"][key])
                if i not in hb_indices
            ]

    cfd = mmcif.CifFileWriter(path)
    cfd.write({component.id: cif_copy})


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

    f_charge = sum([l.GetFormalCharge() for l in mol_to_save.GetAtoms()])
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
    conf_id = (
        0
        if conf_type == ConformerType.Depiction
        else component.conformers_mapping[conf_type]
    )
    mol_to_save = (
        component.mol2D if conf_type == ConformerType.Depiction else component.mol
    )

    if remove_hs:
        mol_to_save = rdkit.Chem.RemoveHs(mol_to_save, sanitize=False)
        mol_to_save.UpdatePropertyCache(strict=False)

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
        component (Component): Component to be
            exported.
    """

    calc_formula = rdkit.Chem.rdMolDescriptors.CalcMolFormula(component.mol)
    calc_weight = rdkit.Chem.rdMolDescriptors.CalcExactMolWt(component.mol)
    date = component.modified_date
    mod_date = f"{date.year}-{date.month:02d}-{date.day:02d}"

    label = "_chem_comp"
    cif_dict[label] = OrderedDict([])
    cif_dict[label]["id"] = component.id
    cif_dict[label]["type"] = "NON-POLYMER"
    cif_dict[label]["pdbx_type"] = "HETAIN"
    cif_dict[label]["formula"] = (
        component.formula if component.formula else calc_formula
    )
    cif_dict[label]["formula_weight"] = f"{calc_weight:.3f}"
    cif_dict[label]["three_letter_code"] = component.id
    cif_dict[label]["pdbx_type"] = "HETAIN"
    cif_dict[label]["pdbx_modified_date"] = mod_date
    cif_dict[label]["pdbx_release_status"] = component.pdbx_release_status.name


def _write_pdb_ccd_cif_atoms(cif_dict, component):
    """Writes the _chem_comp_atom namespace with atom details.
    Controlled dictionary:
    http://mmcif.wwpdb.org/dictionaries/mmcif_pdbx_v50.dic/Categories/chem_comp_atom.html

    Args:
        cif_dict (dict of str: str): cif representation of the molecule
            in a dictionary.
        component (Component): Component to be
            exported.
    """

    if component.mol.GetNumAtoms() < 1:
        return

    # if we dont have ideal coordinates, we just try to compute them
    ideal_type = ConformerType.Ideal
    if component.has_degenerated_conformer(ConformerType.Ideal):
        ideal_type = (
            ConformerType.Computed if component.compute_3d() else ConformerType.Ideal
        )

    label = "_chem_comp_atom"

    cif_dict[label] = OrderedDict([])
    cif_dict[label]["comp_id"] = []
    cif_dict[label]["atom_id"] = []
    cif_dict[label]["alt_atom_id"] = []
    cif_dict[label]["type_symbol"] = []
    cif_dict[label]["charge"] = []
    cif_dict[label]["pdbx_align"] = []
    cif_dict[label]["pdbx_aromatic_flag"] = []
    cif_dict[label]["pdbx_leaving_atom_flag"] = []
    cif_dict[label]["pdbx_stereo_config"] = []
    cif_dict[label]["model_Cartn_x"] = []
    cif_dict[label]["model_Cartn_y"] = []
    cif_dict[label]["model_Cartn_z"] = []
    cif_dict[label]["pdbx_model_Cartn_x_ideal"] = []
    cif_dict[label]["pdbx_model_Cartn_y_ideal"] = []
    cif_dict[label]["pdbx_model_Cartn_z_ideal"] = []
    cif_dict[label]["pdbx_component_atom_id"] = []
    cif_dict[label]["pdbx_component_comp_id"] = []
    cif_dict[label]["pdbx_ordinal"] = []

    for atom in component.mol.GetAtoms():
        at_id = atom.GetIdx()

        model_atom = _get_atom_coord(component, at_id, ConformerType.Model)
        ideal_atom = _get_atom_coord(component, at_id, ideal_type)

        cif_dict[label]["comp_id"].append(component.id)
        cif_dict[label]["atom_id"].append(_get_atom_name(atom))
        cif_dict[label]["alt_atom_id"].append(_get_alt_atom_name(atom))
        cif_dict[label]["type_symbol"].append(atom.GetSymbol())
        cif_dict[label]["charge"].append(str(atom.GetFormalCharge()))
        cif_dict[label]["pdbx_align"].append("?")
        cif_dict[label]["pdbx_aromatic_flag"].append(
            "Y" if atom.GetIsAromatic() else "N"
        )
        cif_dict[label]["pdbx_leaving_atom_flag"].append("N")
        cif_dict[label]["pdbx_stereo_config"].append(_get_ccd_cif_chiral_type(atom))

        cif_dict[label]["model_Cartn_x"].append(f"{model_atom.x:.3f}")
        cif_dict[label]["model_Cartn_y"].append(f"{model_atom.y:.3f}")
        cif_dict[label]["model_Cartn_z"].append(f"{model_atom.z:.3f}")

        cif_dict[label]["pdbx_model_Cartn_x_ideal"].append(f"{ideal_atom.x:.3f}")
        cif_dict[label]["pdbx_model_Cartn_y_ideal"].append(f"{ideal_atom.y:.3f}")
        cif_dict[label]["pdbx_model_Cartn_z_ideal"].append(f"{ideal_atom.z:.3f}")

        cif_dict[label]["pdbx_component_atom_id"].append(_get_atom_name(atom))
        cif_dict[label]["pdbx_component_comp_id"].append(component.id)
        cif_dict[label]["pdbx_ordinal"].append(str(atom.GetIdx() + 1))


def _write_pdb_ccd_cif_bonds(cif_dict, component):
    """Writes the _chem_comp_bond namespace with atom details.
    Controlled dictionary:
    http://mmcif.wwpdb.org/dictionaries/mmcif_pdbx_v50.dic/Categories/chem_comp_bond.html

    Args:
        cif_dict (dict of str: str): cif representation of the molecule
            in a dictionary.
        component (Component): Component to be
            exported.
    """
    if component.mol.GetNumBonds() < 1:
        return

    label = "_chem_comp_bond"
    cif_dict[label] = OrderedDict([])
    cif_dict[label]["comp_id"] = []
    cif_dict[label]["atom_id_1"] = []
    cif_dict[label]["atom_id_2"] = []
    cif_dict[label]["value_order"] = []
    cif_dict[label]["pdbx_aromatic_flag"] = []
    cif_dict[label]["pdbx_stereo_config"] = []
    cif_dict[label]["pdbx_ordinal"] = []

    for b in component.mol.GetBonds():
        atom_a = b.GetBeginAtom()
        atom_b = b.GetEndAtom()

        cif_dict[label]["comp_id"].append(component.id)
        cif_dict[label]["atom_id_1"].append(_get_atom_name(atom_a))
        cif_dict[label]["atom_id_2"].append(_get_atom_name(atom_b))
        cif_dict[label]["value_order"].append(_get_ccd_cif_bond_type(b))
        cif_dict[label]["pdbx_aromatic_flag"].append("Y" if b.GetIsAromatic() else "N")
        cif_dict[label]["pdbx_stereo_config"].append(_get_ccd_cif_bond_stereo(b))
        cif_dict[label]["pdbx_ordinal"].append(str(b.GetIdx() + 1))


def _write_pdb_ccd_cif_descriptor(cif_dict, component):
    """Writes the _pdbx_chem_comp_descriptor namespace with details.

    Args:
        cif_dict (dict of str: str): cif representation of the molecule
            in a dictionary.
        component (Component): Component to be
            exported.
    """

    if not component.descriptors:
        return

    label = "_pdbx_chem_comp_descriptor"
    cif_dict[label] = OrderedDict([])
    cif_dict[label]["comp_id"] = []
    cif_dict[label]["type"] = []
    cif_dict[label]["program"] = []
    cif_dict[label]["descriptor"] = []

    for entry in component.descriptors:
        cif_dict[label]["comp_id"].append(component.id)
        cif_dict[label]["type"].append(entry.type)
        cif_dict[label]["program"].append(entry.program)
        cif_dict[label]["descriptor"].append(entry.value)


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
    elif bond_order == rdkit.Chem.rdchem.BondType.DOUBLE:
        return "2"
    elif bond_order == rdkit.Chem.rdchem.BondType.TRIPLE:
        return "3"
    elif bond_order == rdkit.Chem.rdchem.BondType.AROMATIC:
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
    elif stereo in (
        rdkit.Chem.rdchem.BondStereo.STEREOZ,
        rdkit.Chem.rdchem.BondStereo.STEREOTRANS,
    ):
        return "Z"
    else:
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
    elif bond_order == rdkit.Chem.rdchem.BondType.DOUBLE:
        return "DOUB"
    elif bond_order == rdkit.Chem.rdchem.BondType.TRIPLE:
        return "TRIP"
    elif bond_order == rdkit.Chem.rdchem.BondType.AROMATIC:
        return "AROM"
    elif bond_order == rdkit.Chem.rdchem.BondType.QUADRUPLE:
        return "QUAD"
    else:
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
    elif chiral_type == rdkit.Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CCW:
        return "S"
    else:
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
    conf_id = component.conformers_mapping[conformer_type]

    return component.mol.GetConformer(conf_id).GetAtomPosition(at_id)


# region fallbacks


def _to_sdf_str_fallback(mol, ccd_id, mappings):
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

        content += [
            f"{ccd_id} - {k} conformer",
            "    RDKit   3D",
            "\n"
            f"{atom_count:>3}{bond_count:3}  0  0  0  0  0  0  0  0999 V2000",
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
    elif bond_order == rdkit.Chem.rdchem.BondType.DOUBLE:
        return "2"
    elif bond_order == rdkit.Chem.rdchem.BondType.TRIPLE:
        return "3"
    elif bond_order == rdkit.Chem.rdchem.BondType.AROMATIC:
        return "A"
    else:
        return "0"


def __post_process_cif_category(cif_copy, category_name):
    """Single value category needs to be string rather than array
    with a single value. Also if the category contains nothing, it should
    be removed.

    Args:
        cif_copy (dict of str: dict): Dictionary like structure of
            the CIF file.
        category_name (str): Category name to be accounted for
    """
    if not cif_copy[category_name]:  # nothing in the category => should be removed
        cif_copy.pop(category_name)
        return

    for k, v in cif_copy[category_name].items():
        if isinstance(v, list):
            if len(v) == 1:
                cif_copy[category_name][k] = v[0]

            if not v:
                cif_copy.pop(category_name)
                return


# endregion helpers


# region cif-export addition
def _add_sw_info_cif(cif_copy):
    """Add information to the cif file about software versions used
    to generate this information

    Args:
        cif_copy (dict of str: dict): Dictionary like structure of
            the CIF file.
    """
    category = "_software"

    cif_copy[category] = {}
    cif_copy[category]["name"] = ["rdkit", "pdbeccdutils"]
    cif_copy[category]["version"] = [rdkit.__version__, pdbeccdutils.__version__]
    cif_copy[category]["description"] = [
        "Core functionality.",
        "Wrapper to provide 2D templates and molecular fragments.",
    ]


def _add_2d_depiction_cif(component, cif_copy):
    """Add 2D coordinates of the component depiction

    Args:
        component (Component): pdbeccdutils component.
        cif_copy (dict of str: dict): Dictionary like structure of
            the CIF file.
    """
    if component.mol2D is None:
        return

    __add_rdkit_2d_atoms_cif(component, cif_copy)
    __add_rdkit_2d_bonds_cif(component, cif_copy)


def _add_fragments_and_scaffolds_cif(component, cif_copy):
    """Add fragments and scaffolds information to the CIF export

    Args:
        component (Component): Component to be exported.
        cif_copy (dict of str: dict): Dictionary like structure of
            the CIF file.
    """
    substructure_category = "_pdbe_chem_comp_substructure"
    mapping_category = "_pdbe_chem_comp_substructure_mapping"

    ids = [f"S{i+1}" for i in range(0, len(component.scaffolds))]
    ids += [f"F{i+1}" for i in range(0, len(component.fragments))]

    types = [f"scaffold" for i in range(0, len(component.scaffolds))]
    types += [f"fragment" for i in range(0, len(component.fragments))]

    names = [i.name for i in component.scaffolds]
    names += [i.name for i in component.fragments]

    smiles = [i.smiles for i in component.scaffolds]
    smiles += [i.smiles for i in component.fragments]

    mols = [rdkit.Chem.MolFromSmiles(i.smiles) for i in component.scaffolds]
    mols += [rdkit.Chem.MolFromSmiles(i.smiles) for i in component.fragments]

    inchis = [rdkit.Chem.MolToInchi(i) for i in mols]

    inchikeys = [rdkit.Chem.MolToInchiKey(i) for i in mols]

    # general information about all substructures
    cif_copy[substructure_category] = {}
    cif_copy[mapping_category] = {}

    cif_copy[substructure_category]["comp_id"] = [component.id] * len(ids)
    cif_copy[substructure_category]["substructure_name"] = names
    cif_copy[substructure_category]["id"] = ids
    cif_copy[substructure_category]["substructure_type"] = types
    cif_copy[substructure_category]["substructure_smiles"] = smiles
    cif_copy[substructure_category]["substructure_inchis"] = inchis
    cif_copy[substructure_category]["substructure_inchikeys"] = inchikeys

    __post_process_cif_category(cif_copy, substructure_category)

    cif_copy[mapping_category]["comp_id"] = []
    cif_copy[mapping_category]["atom_id"] = []
    cif_copy[mapping_category]["substructure_id"] = []
    cif_copy[mapping_category]["substructure_ordinal"] = []

    _add_scaffold_cif(component, cif_copy[mapping_category])
    _add_fragments_cif(component, cif_copy[mapping_category])

    __post_process_cif_category(cif_copy, mapping_category)


def _add_fragments_cif(component, cif_copy):
    """Add fragment information to the CIF export.

    Args:
        component (Component): pdbeccdutils component.
        cif_copy (dict of str: dict): Dictionary like structure of
            the _pdbe_chem_comp_substructure_mapping category.
    """
    for i, scaffold in enumerate(component.scaffolds):
        for a, mapping in enumerate(scaffold.mappings):
            for atom_name in mapping:
                cif_copy["comp_id"].append(component.id)
                cif_copy["atom_id"].append(atom_name)
                cif_copy["substructure_id"].append(f"S{i+1}")
                cif_copy["substructure_ordinal"].append(f"{a+1}")


def _add_scaffold_cif(component, cif_copy):
    """Add scaffold information to the cif export.

    Args:
        component (Component): pdbeccdutils component.
        cif_copy (dict of str: dict): Dictionary like structure of
            the CIF file.
    """
    for i, fragment in enumerate(component.fragments):
        for a, mapping in enumerate(fragment.mappings):
            for atom_name in mapping:
                cif_copy["comp_id"].append(component.id)
                cif_copy["atom_id"].append(atom_name)
                cif_copy["substructure_id"].append(f"F{i+1}")
                cif_copy["substructure_ordinal"].append(f"{a+1}")


def _add_rdkit_properties_cif(component, cif_copy):
    """Add properties calculated by the rdkit to the dictionary

    Args:
        component (Component): pdbeccdutils component.
        cif_copy (dict of str: dict): Dictionary like structure of
            the CIF file.
    """
    category = "_pdbe_chem_comp_rdkit_properties"

    cif_copy[category] = {}
    cif_copy[category]["comp_id"] = component.id

    for k, v in component.physchem_properties.items():
        cif_copy[category][k] = f"{v:.0f}" if v.is_integer() else f"{v:.3f}"

    __post_process_cif_category(cif_copy, category)


def __add_rdkit_2d_atoms_cif(component, cif_copy):
    """Add atom coordinates to the mmCIF file.

    Args:
        component (Component): pdbeccdutils component.
        cif_copy (dict of str: dict): Dictionary like structure of
            the CIF file.
    """
    category = "_pdbe_chem_comp_atom_depiction"

    cif_copy[category] = {}
    conformer = component.mol2D.GetConformer()

    cif_copy[category]["comp_id"] = [
        component.id for i in range(0, conformer.GetNumAtoms())
    ]
    cif_copy[category]["atom_id"] = [
        atom.GetProp("name") for atom in component.mol2D.GetAtoms()
    ]
    cif_copy[category]["element"] = [
        atom.GetSymbol() for atom in component.mol2D.GetAtoms()
    ]
    cif_copy[category]["model_Cartn_x"] = [
        f"{conformer.GetAtomPosition(i).x:.3f}"
        for i in range(0, conformer.GetNumAtoms())
    ]
    cif_copy[category]["model_Cartn_y"] = [
        f"{conformer.GetAtomPosition(i).y:.3f}"
        for i in range(0, conformer.GetNumAtoms())
    ]
    cif_copy[category]["pdbx_ordinal"] = [
        i + 1 for i in range(0, conformer.GetNumAtoms())
    ]

    __post_process_cif_category(cif_copy, category)


def __add_rdkit_2d_bonds_cif(component, cif_copy):
    """Add bond information to the mmCIF file.

    Args:
        component (Component): pdbeccdutils component.
        cif_copy (dict of str: dict): Dictionary like structure of
            the CIF file.
    """
    category = "_pdbe_chem_comp_bond_depiction"
    cif_copy[category] = {}

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

    bonds = [b for b in copy.GetBonds() if b.GetEndAtom().GetSymbol() != "H"]

    cif_copy[category]["comp_id"] = [component.id] * len(bonds)
    cif_copy[category]["atom_id_1"] = [b.GetBeginAtom().GetProp("name") for b in bonds]
    cif_copy[category]["atom_id_2"] = [b.GetEndAtom().GetProp("name") for b in bonds]
    cif_copy[category]["value_order"] = [b.GetBondType().name for b in bonds]
    cif_copy[category]["bond_dir"] = [b.GetBondDir().name for b in bonds]
    cif_copy[category]["pdbx_ordinal"] = list(range(1, len(bonds) + 1))

    __post_process_cif_category(cif_copy, category)


def _add_unichem_mapping_cif(component, cif_copy):
    """Add UniChem mapping to CCD CIF.

    Args:
        component (Component): pdbeccdutils component.
        cif_copy (dict of str: dict): Dictionary like structure of
            the CIF file.
    """
    category = "_pdbe_chem_comp_external_mappings"
    cif_copy[category] = {}

    cif_copy[category]["comp_id"] = []
    cif_copy[category]["source"] = []
    cif_copy[category]["resource"] = []
    cif_copy[category]["resource_id"] = []

    for i in component.external_mappings:
        cif_copy[category]["comp_id"].append(component.id)
        cif_copy[category]["source"].append("UniChem")
        cif_copy[category]["resource"].append(i[0])
        cif_copy[category]["resource_id"].append(i[1])

    __post_process_cif_category(cif_copy, category)


def _add_rdkit_conformer_cif(component, cif_copy, remove_hs):
    """Add 3D coordinates generated by RDKit.

    component (Component): pdbeccdutils component.
    cif_copy (dict of str: dict): Dictionary like structure of
        the CIF file.
    remove_hs (boolean): Whether or not hydrogen atoms should be written
    """
    try:
        conformer = component.mol.GetConformer(
            component.conformers_mapping[ConformerType.Computed]
        )
    except ValueError:
        return  # no conformer nothing to write, we quit

    category = "_pdbe_chem_comp_rdkit_conformer"
    cif_copy[category] = {}

    atom_count = (
        sum(1 for atom in component.mol.GetAtoms() if atom.GetSymbol() != "H")
        if remove_hs
        else len(component.mol.GetAtoms())
    )

    # we need to get all the allowed incides
    if remove_hs:
        atom_indices = [
            i
            for i in range(0, component.mol.GetNumAtoms())
            if component.mol.GetAtomWithIdx(i).GetSymbol() != "H"
        ]
    else:
        atom_indices = list(range(0, component.mol.GetNumAtoms()))

    cif_copy[category]["comp_id"] = [component.id] * atom_count
    cif_copy[category]["atom_id"] = [
        component.mol.GetAtomWithIdx(i).GetProp("name") for i in atom_indices
    ]
    cif_copy[category]["Cartn_x_rdkit"] = [
        f"{conformer.GetAtomPosition(i).x:.3f}" for i in atom_indices
    ]
    cif_copy[category]["Cartn_y_rdkit"] = [
        f"{conformer.GetAtomPosition(i).y:.3f}" for i in atom_indices
    ]
    cif_copy[category]["Cartn_z_rdkit"] = [
        f"{conformer.GetAtomPosition(i).z:.3f}" for i in atom_indices
    ]
    cif_copy[category]["rdkit_method"] = ["ETKDGv2"] * atom_count
    cif_copy[category]["rdkit_ordinal"] = [i for i in range(1, atom_count + 1)]

    __post_process_cif_category(cif_copy, category)


# endregion
