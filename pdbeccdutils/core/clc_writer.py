import json
import math
from pathlib import Path

import rdkit
import gemmi
from gemmi import cif
from xml.dom import minidom
from datetime import date as Date
from xml.etree import ElementTree as ET
from xml.etree.ElementTree import Element, SubElement

import pdbeccdutils
from pdbeccdutils.core.component import Component
from pdbeccdutils.core.exceptions import CCDUtilsError
from pdbeccdutils.core.models import ConformerType
from pdbeccdutils.core import ccd_writer
from pdbeccdutils.helpers import mol_tools


def write_molecule(
    path,
    component: Component,
    remove_hs: bool = True,
    conf_type: ConformerType = ConformerType.Model,
):
    """Export molecule in a specified format. Presently supported formats
    are: PDB mmCIF (*.cif); Mol file (*.sdf); Chemical Markup language
    (*.cml); PDB file (*.pdb); XYZ file (*.xyz); XML (*.xml).
    ConformerType.AllConformers is presently supported only for PDB.

    Args:
        path (str|Path): Path to the file. Suffix determines format to be
            used.
        component (Component): Component to be exported
        remove_hs (bool, optional): Defaults to True. Whether or not
            hydrogens should be removed.
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
        str_representation = to_pdb_str(component, remove_hs, conf_type)
    elif extension in ("mmcif", "cif"):
        to_pdb_clc_cif_file(path, component, remove_hs)
        return
    elif extension == "cml":
        str_representation = to_cml_str(component, remove_hs, conf_type)
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
    conf_type: ConformerType = ConformerType.Model,
):
    """Converts structure to the PDB format.

    Args:
        Component (Component): Component to be exported.
        remove_hs (bool, optional): Defaults to True.
        conf_type (ConformerType, optional): Defaults to ConformerType.Ideal.

    Returns:
        str: String representation of the component in the PDB format.
    """
    (mol_to_save, conf_id, conf_type) = ccd_writer._prepate_structure(
        component, remove_hs, conf_type
    )
    residue_numbers = _get_residue_number(mol_to_save)
    for atom in mol_to_save.GetAtoms():
        flag = ccd_writer._get_atom_name(atom)
        atom_name = f"{flag:<4}"  # make sure it is 4 characters
        res_info = atom.GetPDBResidueInfo()
        res_num = residue_numbers[atom.GetProp("residue_id")]
        res_info.SetResidueNumber(res_num)
        res_info.SetName(atom_name)
        res_info.SetTempFactor(20.0)
        res_info.SetOccupancy(1.0)
        res_info.SetChainId("A")

    pdb_title = [
        f"HEADER    {conf_type.name} coordinates",
        f" for PDB-CLC {component.id}",
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


def to_cml_str(component: Component, remove_hs=True, conf_type=ConformerType.Model):
    """Converts structure to the EBI representation of the molecule in
    CML format: http://cml.sourceforge.net/schema/cmlCore.xsd

    Args:
        component (Component): Component to be exported.
        remove_hs (bool, optional): Defaults to True.
        conf_type (ConformerType, optional): Defaults to ConformerType.Ideal.

    Returns:
        str: String representation of the component in CML format.
    """
    (mol_to_save, conf_id, conf_type) = ccd_writer._prepate_structure(
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
    id_inchikey = ET.SubElement(mol, "identifier", {"dictRef": "ebiMolecule:inchikey"})
    id_inchikey.text = component.inchikey
    id_formula = ET.SubElement(mol, "formula", {"dictRef": "ebiMolecule:stereoSmiles"})
    id_formula.text = next(
        (
            x.value
            for x in component.descriptors
            if x.type == "SMILES" and x.program == "RDKIT"
        ),
        "",
    )

    atom_array = ET.SubElement(mol, "atomArray")
    conformer = mol_to_save.GetConformer(id=conf_id)

    for atom in mol_to_save.GetAtoms():
        element = atom.GetSymbol()
        a_name = ccd_writer._get_atom_name(atom)
        coords = conformer.GetAtomPosition(atom.GetIdx())

        a_entry = ET.SubElement(
            atom_array, "atom", {"id": a_name, "elementType": element}
        )
        a_entry.set("x3", str(coords.x))
        a_entry.set("y3", str(coords.y))
        a_entry.set("z3", str(coords.z))

    bond_array = ET.SubElement(mol, "bondArray")
    for bond in mol_to_save.GetBonds():
        atom_1 = ccd_writer._get_atom_name(bond.GetBeginAtom())
        atom_2 = ccd_writer._get_atom_name(bond.GetEndAtom())
        bond_order = ccd_writer._get_cml_bond_type(bond.GetBondType())

        bond_entry = ET.SubElement(bond_array, "bond")
        bond_entry.set("atomsRefs2", atom_1 + " " + atom_2)
        bond_entry.set("order", str(bond_order))

    cml = ET.tostring(root, encoding="utf-8", method="xml")
    pretty = minidom.parseString(cml)

    return pretty.toprettyxml(indent="  ")


def to_xml_xml(component):
    """Converts structure to the XML format and returns its XML repr.

    Args:
        component (Component): Component to be exported.

    Returns:
        xml.etree.ElementTree.Element: XML object
    """
    root = Element("chemComp")

    id_e = SubElement(root, "id")
    formula_e = SubElement(root, "formula")
    s_smiles_e = SubElement(root, "stereoSmiles")
    n_smiles_e = SubElement(root, "nonStereoSmiles")
    inchi_e = SubElement(root, "inchi")

    id_e.text = component.id
    formula_e.text = component.formula
    s_smiles_e.text = next(
        (
            x.value
            for x in component.descriptors
            if x.type == "SMILES_CANONICAL" and x.program == "RDKit"
        ),
        "",
    )
    n_smiles_e.text = next(
        (
            x.value
            for x in component.descriptors
            if x.type == "SMILES" and x.program == "RDKit"
        ),
        "",
    )
    inchi_e.text = component.inchi

    return root


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
        f"HEADER    {conf_name} coordinates for PDB-CLC {component_id}",
        f"COMPND    {component_id}",
        f"AUTHOR    pdbccdutils {pdbeccdutils.__version__}",
        f"AUTHOR    RDKit {rdkit.__version__}",
    ]

    if conf_id == -1:
        conformer_ids = [c.GetId() for c in mol.GetConformers()]
    else:
        conformer_ids = [conf_id]

    residue_numbers = _get_residue_number(mol)
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
                residue_numbers[atom.GetProp("residue_id")],
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


def to_pdb_clc_cif_file(path, component: Component, remove_hs=True):
    """Converts structure to the PDB mmCIF format.
    Args:
        path (str): Path to save cif file.
        component (Component): Component to be exported.
        remove_hs (bool, optional): Defaults to True.
    """

    if not isinstance(component.ccd_cif_block, gemmi.cif.Block):
        component.ccd_cif_block = _to_pdb_clc_cif_block(component)

    temp_doc = gemmi.cif.Document()
    cif_block_copy = temp_doc.add_copied_block(component.ccd_cif_block)

    ccd_writer._add_sw_info_cif(cif_block_copy)
    ccd_writer._add_2d_depiction_cif(component, cif_block_copy)
    ccd_writer._add_fragments_and_scaffolds_cif(component, cif_block_copy)
    ccd_writer._add_rdkit_properties_cif(component, cif_block_copy)
    ccd_writer._add_unichem_mapping_cif(component, cif_block_copy)
    ccd_writer._add_rdkit_conformer_cif(component, cif_block_copy, remove_hs)

    if remove_hs:
        ccd_writer.remove_hydrogens(cif_block_copy)

    temp_doc.write_file(path, gemmi.cif.Style.Pdbx)


def _to_pdb_clc_cif_block(component):
    """Export component to the PDB CCD CIF file.

    Args:
        component (pdbeccdutils.core.Component): Component to be
            exported.

    Returns:
        (dict of str: str)): dictionary representation of the component
            to be serialized as a cif file.
    """

    doc = gemmi.cif.Document()
    cif_block = doc.add_new_block(component.id)

    _write_pdb_clc_cif_info(cif_block, component)
    _write_pdb_clc_cif_atoms(cif_block, component)
    _write_pdb_clc_cif_bonds(cif_block, component)
    _write_pdb_clc_cif_descriptor(cif_block, component)

    return cif_block


def _write_pdb_clc_cif_info(cif_block, component):
    """Writes _chem_comp namespace with general information about the
    component

    Args:
        cif_block (cif.Block): mmcif Block object from gemmi.
        component (Component): Component to be exported.
    """

    calc_formula = rdkit.Chem.rdMolDescriptors.CalcMolFormula(component.mol)
    calc_weight = rdkit.Chem.rdMolDescriptors.CalcExactMolWt(component.mol)
    date = component.modified_date
    mod_date = (
        f"{date.year}-{date.month:02d}-{date.day:02d}"
        if isinstance(date, Date)
        else date
    )

    label = "_chem_comp."
    cif_block.set_pairs(
        label,
        {
            "id": component.id,
            "type": "NON-POLYMER",
            "pdbx_type": "HETAIN",
            "formula": component.formula or calc_formula,
            "formula_weight": f"{calc_weight:.3f}",
            "pdbx_modified_date": mod_date,
            "pdbx_release_status": component.pdbx_release_status.name,
            "pdbx_processing_site": "PDBe",
        },
        raw=False,
    )


def _write_pdb_clc_cif_atoms(cif_block, component):
    """Writes the _chem_comp_atom namespace with atom details.
    Controlled dictionary:
    http://mmcif.wwpdb.org/dictionaries/mmcif_pdbx_v50.dic/Categories/chem_comp_atom.html

    Args:
        cif_block (cif.Block): mmcif Block object from gemmi.
        component (Component): Component to be exported.
    """

    if component.mol.GetNumAtoms() < 1:
        return

    label = "_chem_comp_atom."
    atom_fields = [
        "comp_id",
        "atom_id",
        "alt_atom_id",
        "type_symbol",
        "charge",
        "pdbx_aromatic_flag",
        "pdbx_leaving_atom_flag",
        "pdbx_stereo_config",
        "model_Cartn_x",
        "model_Cartn_y",
        "model_Cartn_z",
        "pdbx_component_comp_id",
        "pdbx_residue_numbering",
        "pdbx_component_atom_id",
        "pdbx_ordinal",
    ]
    atom_loop = cif_block.init_loop(label, atom_fields)
    residue_numbers = _get_residue_number(component.mol)
    for atom in component.mol.GetAtoms():
        at_id = atom.GetIdx()
        model_atom = ccd_writer._get_atom_coord(component, at_id, ConformerType.Model)
        res_info = atom.GetPDBResidueInfo()
        new_row = [
            component.id,
            cif.as_string(ccd_writer._get_atom_name(atom)),
            cif.as_string(ccd_writer._get_atom_name(atom)),
            atom.GetSymbol(),
            str(atom.GetFormalCharge()),
            "Y" if atom.GetIsAromatic() else "N",
            "N",
            ccd_writer._get_ccd_cif_chiral_type(atom),
            f"{model_atom.x:.3f}",
            f"{model_atom.y:.3f}",
            f"{model_atom.z:.3f}",
            res_info.GetResidueName(),
            residue_numbers[atom.GetProp("residue_id")],
            mol_tools.get_component_atom_id(atom),
            str(atom.GetIdx() + 1),
        ]

        atom_loop.add_row(cif.quote_list(new_row))


def _write_pdb_clc_cif_bonds(cif_block, component):
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
            cif.as_string(ccd_writer._get_atom_name(atom_a)),
            cif.as_string(ccd_writer._get_atom_name(atom_b)),
            ccd_writer._get_ccd_cif_bond_type(b),
            "Y" if b.GetIsAromatic() else "N",
            ccd_writer._get_ccd_cif_bond_stereo(b),
            str(b.GetIdx() + 1),
        ]
        bond_loop.add_row(gemmi.cif.quote_list(new_row))


def _write_pdb_clc_cif_descriptor(cif_block, component):
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


def _get_residue_number(mol):
    """Maps residue ids of chemcial components from protein
    to bound-molecule

    Args:
        component: Component

    Returns:
        A dictionary of mappings of residue ids
    """
    residue_id_mapping = {}
    res_num = 0
    for atom in mol.GetAtoms():
        residue_id = atom.GetProp("residue_id")
        if residue_id not in residue_id_mapping:
            res_num += 1
            residue_id_mapping[residue_id] = res_num

    return residue_id_mapping
