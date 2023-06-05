import json
import math
import rdkit
from pathlib import Path
from datetime import date
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
from pdbeccdutils.core.ccd_reader import CCDReaderResult
from pdbeccdutils.core import ccd_writer, bm_writer
from pdbeccdutils.core.component import Component
from pdbeccdutils.core.exceptions import CCDUtilsError
from pdbeccdutils.core.models import (
    CCDProperties,
    ConformerType,
    Descriptor,
    ReleaseStatus,
)
import pdbeccdutils
from pdbeccdutils.helpers import mol_tools
from gemmi import cif


def inchi_to_component(
    inchi: str, comp_id: str, comp_name: str, sanitize=True, addHs=False
) -> Component:
    """Returns Component representation from InChI

    Args:
        inchi: InChI of a molecule
    """

    warnings = []
    errors = []
    mol_result = mol_tools.mol_from_inchi(inchi)
    mol = mol_result.mol
    if not mol:
        raise CCDUtilsError(
            f"rdkit.Chem.rdchem.Mol can not be generated from inchi due to {mol_result.errors}"
        )

    if mol_result.warnings:
        warnings.append(mol_result.warnings)

    if mol_result.errors:
        errors.append(mol_result.errors)

    if addHs:
        mol = rdkit.Chem.AddHs(mol, addCoords=True)

    mol = _setup_conformer(mol)

    sanitized = False
    if sanitize:
        sanitized = mol_tools.sanitize(mol)

    for atom in mol.GetAtoms():
        atom_name = ccd_writer._get_atom_name(atom)
        atom.SetProp("name", atom_name)

    inchi_result = mol_tools.inchi_from_mol(mol)
    if inchi_result.warnings:
        warnings.append(inchi_result.warnings)
    if inchi_result.errors:
        errors.append(inchi_result.errors)

    inchikey = mol_tools.inchikey_from_inchi(inchi_result.inchi)
    descriptors = [
        Descriptor(
            type="SMILES",
            program="rdkit",
            program_version=rdkit.__version__,
            value=rdkit.Chem.MolToSmiles(mol),
        ),
        Descriptor(
            type="InChI",
            program="rdkit",
            program_version=rdkit.__version__,
            value=inchi_result.inchi,
        ),
        Descriptor(
            type="InChIKey",
            program="rdkit",
            program_version=rdkit.__version__,
            value=inchikey,
        ),
    ]
    properties = CCDProperties(
        id=comp_id,
        name=comp_name,
        formula=CalcMolFormula(mol),
        modified_date=date.today(),
        pdbx_release_status=ReleaseStatus.NOT_SET,
        weight="",
    )

    comp = Component(mol, None, properties, descriptors)

    comp.ccd_cif_block = _to_pdb_cif_block(comp)

    reader_result = CCDReaderResult(
        warnings=warnings, errors=errors, component=comp, sanitized=sanitized
    )

    return reader_result


def _setup_conformer(mol) -> rdkit.Chem.rdchem.Conformer:
    """Setup model conformer in the rdkit Mol object

    Args:
        atoms: atoms of bound-molecule

    Returns:
        RDKit generated Conformer of the component.
    """
    options = rdkit.Chem.AllChem.ETKDGv3()
    mol = rdkit.Chem.AddHs(mol)
    conf_id = rdkit.Chem.AllChem.EmbedMolecule(mol, options)
    rdkit.Chem.AllChem.UFFOptimizeMolecule(mol, confId=conf_id, maxIters=1000)
    conformer = mol.GetConformer(conf_id)
    conformer.SetProp("name", ConformerType.Ideal.name)
    conformer.SetProp("coord_generation", "ETKDGv3")
    return mol


def _to_pdb_cif_block(component):
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

    _write_pdb_cif_info(cif_block, component)
    _write_pdb_cif_atoms(cif_block, component)
    _write_pdb_cif_bonds(cif_block, component)
    _write_pdb_cif_descriptor(cif_block, component)

    return cif_block


def _write_pdb_cif_info(cif_block, component):
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
            "pdbx_type": "HETAIN",
            "formula": component.formula or calc_formula,
            "formula_weight": f"{calc_weight:.3f}",
            "pdbx_modified_date": mod_date,
            "pdbx_release_status": component.pdbx_release_status.name,
            "pdbx_ideal_coordinates_details": "RDKit",
            "pdbx_processing_site": None,
        },
        raw=False,
    )


def _write_pdb_cif_atoms(cif_block, component):
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
    ideal_type = (
        ConformerType.Ideal
        if not component.has_degenerated_conformer(ConformerType.Ideal)
        else None
    )

    label = "_chem_comp_atom."
    atom_fields = [
        "comp_id",
        "atom_id",
        "type_symbol",
        "charge",
        "pdbx_aromatic_flag",
        "pdbx_leaving_atom_flag",
        "pdbx_stereo_config",
        "pdbx_model_Cartn_x_ideal",
        "pdbx_model_Cartn_y_ideal",
        "pdbx_model_Cartn_z_ideal",
        "pdbx_ordinal",
    ]
    atom_loop = cif_block.init_loop(label, atom_fields)

    for atom in component.mol.GetAtoms():
        at_id = atom.GetIdx()
        ideal_atom = (
            ccd_writer._get_atom_coord(component, at_id, ideal_type)
            if ideal_type
            else None
        )
        new_row = [
            component.id,
            bm_writer._get_atom_name(atom),
            atom.GetSymbol(),
            str(atom.GetFormalCharge()),
            "Y" if atom.GetIsAromatic() else "N",
            "N",
            ccd_writer._get_ccd_cif_chiral_type(atom),
            f"{ideal_atom.x:.3f}" if ideal_atom else None,
            f"{ideal_atom.y:.3f}" if ideal_atom else None,
            f"{ideal_atom.z:.3f}" if ideal_atom else None,
            str(atom.GetIdx() + 1),
        ]

        atom_loop.add_row(cif.quote_list(new_row))


def _write_pdb_cif_bonds(cif_block, component):
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
            cif.as_string(bm_writer._get_atom_name(atom_a)),
            cif.as_string(bm_writer._get_atom_name(atom_b)),
            ccd_writer._get_ccd_cif_bond_type(b),
            "Y" if b.GetIsAromatic() else "N",
            ccd_writer._get_ccd_cif_bond_stereo(b),
            str(b.GetIdx() + 1),
        ]
        bond_loop.add_row(cif.quote_list(new_row))


def _write_pdb_cif_descriptor(cif_block, component):
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


def write_molecule(
    path,
    component: Component,
    remove_hs: bool = True,
    conf_type: ConformerType = ConformerType.Ideal,
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
        to_pdb_cif_file(path, component, remove_hs)
        return
    elif extension == "cml":
        str_representation = bm_writer.to_cml_str(component, remove_hs, conf_type)
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


def to_pdb_cif_file(path, component: Component, remove_hs=True):
    """Converts structure to the PDB mmCIF format.
    Args:
        path (str): Path to save cif file.
        component (Component): Component to be exported.
        remove_hs (bool, optional): Defaults to True.
    """

    if not isinstance(component.ccd_cif_block, cif.Block):
        component.ccd_cif_block = _to_pdb_cif_block(component)

    temp_doc = cif.Document()
    cif_block_copy = temp_doc.add_copied_block(component.ccd_cif_block)

    ccd_writer._add_sw_info_cif(cif_block_copy)
    ccd_writer._add_2d_depiction_cif(component, cif_block_copy)
    ccd_writer._add_fragments_and_scaffolds_cif(component, cif_block_copy)
    ccd_writer._add_rdkit_properties_cif(component, cif_block_copy)
    ccd_writer._add_unichem_mapping_cif(component, cif_block_copy)

    if remove_hs:
        ccd_writer.remove_hydrogens(cif_block_copy)

    temp_doc.write_file(path, cif.Style.Pdbx)


def to_pdb_str(
    component: Component,
    remove_hs: bool = True,
    conf_type: ConformerType = ConformerType.Computed,
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

    info = rdkit.Chem.rdchem.AtomPDBResidueInfo()
    info.SetResidueName(f"{component.id:>3}")
    info.SetTempFactor(20.0)
    info.SetOccupancy(1.0)
    info.SetChainId("A")
    info.SetResidueNumber(1)
    info.SetIsHeteroAtom(True)

    for atom in mol_to_save.GetAtoms():
        flag = bm_writer._get_atom_name(atom)
        atom_name = f"{flag:<4}"  # make sure it is 4 characters
        info.SetName(atom_name)
        atom.SetMonomerInfo(info)

    pdb_title = [
        f"HEADER    {conf_type.name} coordinates",
        f" for compound {component.id}",
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


def _to_pdb_str_fallback(mol, component_id, conf_id, conf_name="Computed"):
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
        f"HEADER    {conf_name} coordinates for compound {component_id}",
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
