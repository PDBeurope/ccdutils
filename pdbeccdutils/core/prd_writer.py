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
import math
from pathlib import Path
import pdbeccdutils
import rdkit
import gemmi
from gemmi import cif
from pdbeccdutils.core import ccd_writer
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
        str_representation = ccd_writer.to_sdf_str(component, remove_hs, conf_type)
    elif extension == "pdb":
        str_representation = to_pdb_str(component, remove_hs, alt_names, conf_type)
    elif extension in ("mmcif", "cif"):
        to_pdb_ccd_cif_file(path, component, remove_hs)
        return
    elif extension == "cml":
        str_representation = ccd_writer.to_cml_str(component, remove_hs, conf_type)
    elif extension == "xml":
        str_representation = ccd_writer.to_xml_str(component, remove_hs, conf_type)
    elif extension == "xyz":
        str_representation = ccd_writer.to_xyz_str(component, remove_hs, conf_type)
    elif extension == "json":
        str_representation = json.dumps(
            ccd_writer.to_json_dict(component, remove_hs, conf_type),
            sort_keys=True,
            indent=4,
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
    (mol_to_save, conf_id, conf_type) = ccd_writer._prepate_structure(
        component, remove_hs, conf_type
    )

    for atom in mol_to_save.GetAtoms():
        flag = (
            ccd_writer._get_alt_atom_name(atom)
            if alt_names
            else ccd_writer._get_atom_name(atom)
        )
        atom_name = f"{flag:<4}"  # make sure it is 4 characters
        res_info = atom.GetPDBResidueInfo()
        res_info.SetName(atom_name)
        res_info.SetTempFactor(20.0)
        res_info.SetOccupancy(1.0)
        res_info.SetChainId("A")

    pdb_title = [
        f"HEADER    {conf_type.name} coordinates",
        f" for PDB-PRD {component.id}",
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

    ccd_writer._add_sw_info_cif(cif_block_copy)
    ccd_writer._add_2d_depiction_cif(component, cif_block_copy)
    ccd_writer._add_fragments_and_scaffolds_cif(component, cif_block_copy)
    ccd_writer._add_rdkit_properties_cif(component, cif_block_copy)
    ccd_writer._add_unichem_mapping_cif(component, cif_block_copy)
    ccd_writer._add_rdkit_conformer_cif(component, cif_block_copy, remove_hs)

    if remove_hs:
        ccd_writer.remove_hydrogens(cif_block_copy)

    temp_doc.write_file(path, cif.Style.Pdbx)


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
    ccd_writer._write_pdb_ccd_cif_bonds(cif_block, component)
    ccd_writer._write_pdb_ccd_cif_descriptor(cif_block, component)

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

    label = "_chem_comp."
    cif_block.set_pairs(
        label,
        {
            "id": component.id,
            "formula": component.formula or calc_formula,
            "formula_weight": f"{calc_weight:.3f}",
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
        "pdbx_component_comp_id",
        "pdbx_residue_numbering",
        "pdbx_component_atom_id",
        "pdbx_polymer_type",
        "pdbx_ref_id",
        "pdbx_component_id",
        "pdbx_ordinal",
    ]
    atom_loop = cif_block.init_loop(label, atom_fields)

    for atom in component.mol.GetAtoms():
        at_id = atom.GetIdx()
        model_atom = ccd_writer._get_atom_coord(component, at_id, ConformerType.Model)
        ideal_atom = ccd_writer._get_atom_coord(component, at_id, ideal_type)
        res_info = atom.GetPDBResidueInfo()

        new_row = [
            component.id,
            ccd_writer._get_atom_name(atom),
            ccd_writer._get_alt_atom_name(atom),
            atom.GetSymbol(),
            str(atom.GetFormalCharge()),
            None,
            "Y" if atom.GetIsAromatic() else "N",
            "N",
            ccd_writer._get_ccd_cif_chiral_type(atom),
            f"{model_atom.x:.3f}",
            f"{model_atom.y:.3f}",
            f"{model_atom.z:.3f}",
            f"{ideal_atom.x:.3f}",
            f"{ideal_atom.y:.3f}",
            f"{ideal_atom.z:.3f}",
            res_info.GetResidueName(),
            _get_residue_id(atom),
            _get_component_atom_id(atom),
            _get_res_type(atom),
            _get_ref_id(atom),
            _get_component_id(atom),
            str(atom.GetIdx() + 1),
        ]

        atom_loop.add_row(cif.quote_list(new_row))


def _get_component_id(atom):
    """Gets component id. If not returns None.

    Args:
        atom (rdkit.Chem.rdchem.Atom): rdkit atom.

    Returns:
        str: component id the atom.
    """
    return atom.GetProp("comp_id") if atom.HasProp("comp_id") else None


def _get_ref_id(atom):
    """Gets ref id. If not returns None.

    Args:
        atom (rdkit.Chem.rdchem.Atom): rdkit atom.

    Returns:
        str: ref id of the atom.
    """
    return atom.GetProp("ref_id") if atom.HasProp("ref_id") else None


def _get_res_type(atom):
    """Gets res type. If not returns None.

    Args:
        atom (rdkit.Chem.rdchem.Atom): rdkit atom.

    Returns:
        str: res type of the atom.
    """
    return atom.GetProp("res_type") if atom.HasProp("res_type") else None


def _get_component_atom_id(atom):
    """Gets component atom id. If not returns atom name is used.

    Args:
        atom (rdkit.Chem.rdchem.Atom): rdkit atom.

    Returns:
        str: component atom id.
    """
    return (
        atom.GetProp("component_atom_id")
        if atom.HasProp("component_atom_id")
        else ccd_writer._get_atom_name(atom)
    )


def _get_residue_id(atom):
    """Returns residue id. If not returns None

    Args:
        atom (rdkit.Chem.rdchem.Atom): rdkit atom.
    """
    return atom.GetProp("residue_id") if atom.HasProp("residue_id") else None


# region fallbacks


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
            res_info = atom.GetPDBResidueInfo()

            s = "{:<6}{:>5} {:<4} {:>3} {}{:>4}{}   {:>8.3f}{:>8.3f}{:>8.3f}{:>6.2f}{:>6.2f}          {:>2}{:>2}".format(
                "HETATM",
                i + 1,
                atom.GetProp("name"),
                res_info.GetResidueName(),
                "A",
                res_info.GetResidueNumber(),
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
