import io
import os
import rdkit
import re

from rdkit.Chem.rdMolDescriptors import CalcMolFormula
from pdbeccdutils.core import models
from pdbeccdutils.core import ccd_reader
from pdbeccdutils.utils import config
from pdbeccdutils.core.component import Component
from collections import namedtuple

from pdbeccdutils.core.models import (
    CCDProperties,
    ConformerType,
    Descriptor,
    Residue,
    BoundMolecule,
)
from pdbeccdutils.helpers import cif_tools, conversions, mol_tools, helper
from gemmi import cif
from networkx import MultiDiGraph, connected_components
from contextlib import redirect_stderr


BMReaderResult = namedtuple(
    "BMReaderResult", ccd_reader.CCDReaderResult._fields + ("bound_molecule",)
)


def read_pdb_updated_cif_file(
    path_to_cif: str,
    pdb_id: str,
    to_discard: set[str] = config.DISCARDED_RESIDUES,
    sanitize: bool = True,
):
    """
    Read in single wwPDB Model CIF and create internal
    representation of its bound-molecules with multiple components.

    Args:
        path_to_cif (str): Path to the cif file
        sanitize (bool): [Defaults: True]

    Raises:
        ValueError: if file does not exist

    Returns:
        A list of CCDResult representations of each bound-molecule.
    """
    if not os.path.isfile(path_to_cif):
        raise ValueError(f"File '{path_to_cif}' does not exists")

    biomolecule_result = []
    bms = infer_bound_molecules(path_to_cif, to_discard)
    for i, bm in enumerate(bms, start=1):
        rdkit_stream = io.StringIO()
        with redirect_stderr(rdkit_stream):
            bm_id = f"{pdb_id}-bm{i}"
            reader_result = infer_multiple_chem_comp(path_to_cif, bm, bm_id, sanitize)
            if reader_result:
                rdkit_log = rdkit_stream.getvalue()
                if "WARNING" in rdkit_log:
                    start_index = re.search(r"\bWARNING\b:", rdkit_log).end()
                    rdkit_log = rdkit_log[start_index:].strip()
                    reader_result.warnings.append(rdkit_log)
                elif "ERROR" in rdkit_log:
                    start_index = re.search(r"\bERROR\b:", rdkit_log).end()
                    rdkit_log = rdkit_log[start_index:].strip()
                    reader_result.errors.append(rdkit_log)

                biomolecule_result.append(reader_result)

    return biomolecule_result


def infer_multiple_chem_comp(path_to_cif, bm, bm_id, sanitize=True):

    if bm.graph.number_of_nodes() <= 1:
        return

    cif_block = cif.read(path_to_cif).sole_block()
    (mol, warnings, errors) = _parse_pdb_mmcif(cif_block, bm.graph)
    sanitized = False
    if sanitize:
        sanitized = mol_tools.sanitize(mol)

    comp = Component(mol.GetMol(), cif_block)
    descriptors = [
        Descriptor(type="InChI", program="rdkit", value=comp.inchi_from_rdkit),
        Descriptor(type="InChIKey", program="rdkit", value=comp.inchikey_from_rdkit),
    ]
    bm_name = "_".join([residue.name for residue in bm.graph.nodes()])
    properties = CCDProperties(
        id=bm_id,
        name=bm_name,
        formula=CalcMolFormula(comp.mol),
        modified_date="",
        pdbx_release_status=models.ReleaseStatus.NOT_SET,
        weight="",
    )
    comp = Component(mol.GetMol(), cif_block, properties, descriptors)
    reader_result = BMReaderResult(
        warnings=warnings,
        errors=errors,
        component=comp,
        bound_molecule=bm,
        sanitized=sanitized,
    )

    return reader_result


def _parse_pdb_mmcif(cif_block, bm):
    """
    Create internal representation of the molecule from mmcif format.

    Args:
        cif_block (cif.Block): mmcif block object from gemmi
        sanitize (bool): Whether or not the rdkit component should
            be sanitized. Defaults to True.

    Returns:
        CCDReaderResult: internal representation with the results
            of parsing and Mol object.
    """
    warnings = []
    errors = []
    mol = rdkit.Chem.RWMol()
    preprocessable_categories = ["_atom_site.", "_chem_comp_bond."]

    for c in preprocessable_categories:
        w = cif_tools.preprocess_cif_category(cif_block, c)

        if w:
            warnings.append(w)

    bm_atoms = _get_boundmolecule_atoms(cif_block, bm)
    _parse_pdb_atoms(mol, bm_atoms)
    _parse_pdb_conformers(mol, bm_atoms)
    _parse_pdb_bonds(mol, bm, cif_block, errors)
    _add_connections(mol, bm, errors)
    _remove_disconnected_hydrogens(mol)
    return (mol, warnings, errors)


def _get_boundmolecule_atoms(cif_block, bm):
    if "_atom_site." not in cif_block.get_mmcif_category_names():
        return

    atoms = cif_block.get_mmcif_category("_atom_site.")
    bm_atoms = {key: [] for key in atoms}
    for i in range(len(atoms["id"])):
        if atoms["group_PDB"][i] == "HETATM":
            for residue in bm.nodes():
                if (
                    atoms["label_comp_id"][i] == residue.name
                    and atoms["auth_asym_id"][i] == residue.chain
                    and atoms["auth_seq_id"][i] == residue.res_id
                ):
                    for key in bm_atoms:
                        bm_atoms[key].append(atoms[key][i])

    return bm_atoms


def _parse_pdb_atoms(mol, atoms):
    for i in range(len(atoms["id"])):
        atom_id = atoms["label_atom_id"][i]
        chain = atoms["auth_asym_id"][i]
        res_name = atoms["label_comp_id"][i]
        res_id = atoms["auth_seq_id"][i]
        ins_code = (
            "" if not atoms["pdbx_PDB_ins_code"][i] else atoms["pdbx_PDB_ins_code"][i]
        )
        residue_id = f"{chain}{res_id}{ins_code}"
        element = atoms["type_symbol"][i]
        element = element if len(element) == 1 else element[0] + element[1].lower()
        isotope = None
        if element == "D":
            element = "H"
            isotope = 2
        elif element == "X":
            element = "*"

        atom = rdkit.Chem.Atom(element)
        atom.SetProp("name", atom_id)
        atom.SetProp("residue", res_name)
        atom.SetProp("residue_id", residue_id)

        if isotope is not None:
            atom.SetIsotope(isotope)

        mol.AddAtom(atom)


def _parse_pdb_conformers(mol, atoms):
    if not atoms:
        return

    model = _setup_pdb_conformer(atoms)
    mol.AddConformer(model, assignId=True)


def _setup_pdb_conformer(atoms):
    if not atoms:
        return

    conformer = rdkit.Chem.Conformer(len(atoms["id"]))
    for i in range(len(atoms["id"])):
        x = conversions.str_to_float(atoms["Cartn_x"][i])
        y = conversions.str_to_float(atoms["Cartn_y"][i])
        z = conversions.str_to_float(atoms["Cartn_z"][i])
        atom_position = rdkit.Chem.rdGeometry.Point3D(x, y, z)
        conformer.SetAtomPosition(i, atom_position)

    conformer.SetProp("name", ConformerType.Model.name)
    return conformer


def _parse_pdb_bonds(mol, bm, cif_block, errors):
    if (
        "_atom_site." not in cif_block.get_mmcif_category_names()
        or "_chem_comp_bond." not in cif_block.get_mmcif_category_names()
    ):
        return

    for residue in bm.nodes():
        resiude_bonds = get_chem_comp_bonds(cif_block, residue.name)
        for i in range(len(resiude_bonds.atom_id_1)):
            try:
                atom_1 = resiude_bonds.atom_id_1[i]
                mol_atom_1_idx = helper.find_atom_index(mol, residue.id, atom_1)
                atom_2 = resiude_bonds.atom_id_2[i]
                mol_atom_2_idx = helper.find_atom_index(mol, residue.id, atom_2)
                bond_order = helper.bond_pdb_order(resiude_bonds.value_order[i])
                if (mol_atom_1_idx is not None) and (mol_atom_2_idx is not None):
                    mol.AddBond(mol_atom_1_idx, mol_atom_2_idx, bond_order)
            except ValueError:
                errors.append(
                    f"Error perceiving {atom_1} - {atom_2} bond from _chem_comp_bond"
                )
            except RuntimeError:
                errors.append(f"Duplicit bond {atom_1} - {atom_2}")


def _add_connections(mol, bm, errors):
    for residue_1, residue_2, atoms in bm.edges(data=True):
        try:
            atom_1 = atoms["atom_id_1"]
            mol_atom_1_idx = helper.find_atom_index(mol, residue_1.id, atom_1)
            atom_2 = atoms["atom_id_2"]
            mol_atom_2_idx = helper.find_atom_index(mol, residue_2.id, atom_2)
            bond_order = helper.bond_pdb_order("SING")
            if (mol_atom_1_idx is not None) and (mol_atom_2_idx is not None):
                mol.AddBond(mol_atom_1_idx, mol_atom_2_idx, bond_order)
        except ValueError:
            errors.append(
                f"Error perceiving {atom_1} - {atom_2} bond from Boundmolecule connections"
            )
        except RuntimeError:
            errors.append(f"Duplicit bond {atom_1} - {atom_2}")


def get_chem_comp_bonds(cif_block, residue):
    if "_chem_comp_bond." not in cif_block.get_mmcif_category_names():
        return
    chem_comp_bonds = cif_block.get_mmcif_category("_chem_comp_bond.")
    ResidueBonds = namedtuple("ResidueBonds", "residue atom_id_1 atom_id_2 value_order")
    atom_id_1 = []
    atom_id_2 = []
    value_order = []
    last_comp_id = None
    residue_found = False
    for i in range(len(chem_comp_bonds["comp_id"])):
        chem_comp_id = chem_comp_bonds["comp_id"][i]
        if chem_comp_id == residue:
            residue_found = True
            atom_id_1.append(chem_comp_bonds["atom_id_1"][i])
            atom_id_2.append(chem_comp_bonds["atom_id_2"][i])
            value_order.append(chem_comp_bonds["value_order"][i])
        last_comp_id = chem_comp_id
        if last_comp_id != residue and residue_found:
            break

    residue_bonds = ResidueBonds(residue, atom_id_1, atom_id_2, value_order)
    return residue_bonds


def _remove_disconnected_hydrogens(mol):
    """
    Delete hydrogens without neighbours (degree = 0).
    RDKit works but gives warning for these ("WARNING: not removing hydrogen atom without neighbors").
    However we wouldn't want these hydrogens as they affect molecular properties.

    Args:
        mol (rdkit.Chem.rchem.Mol): Rdkit Mol object with the
            compound representation.
    """
    disconnected_hydrogen_indices = [
        atom.GetIdx()
        for atom in mol.GetAtoms()
        if atom.GetAtomicNum() == 1 and atom.GetDegree() == 0
    ]
    disconnected_hydrogen_indices.sort(reverse=True)

    for atom in disconnected_hydrogen_indices:
        mol.RemoveAtom(atom)


def infer_bound_molecules(structure, to_discard):
    """Identify bound molecules in the input protein structure.

    Args:
        structure (str): Path to the structure.
        to_discard (list of str): List of residue names to be discarded
    """

    bms = parse_bound_molecules(structure, to_discard)
    bound_molecules = []

    for bm_nodes in connected_components(bms.to_undirected()):
        subgraph = bms.subgraph(bm_nodes)
        bm = BoundMolecule(subgraph)
        bound_molecules.append(bm)

    bound_molecules = sorted(bound_molecules, key=lambda l: -len(l.graph.nodes))
    return bound_molecules


def __add_connections(struct_conn, bms):
    for i in range(len(struct_conn["id"])):
        (ptnr1, atom1) = find_pntr_entry(struct_conn, bms.nodes, 1, i)
        (ptnr2, atom2) = find_pntr_entry(struct_conn, bms.nodes, 2, i)

        # we want covalent connections among ligands only.
        if ptnr1 and ptnr2 and struct_conn["conn_type_id"][i] != "metalc":
            add_edge = True
            edge_data = bms.get_edge_data(ptnr1, ptnr2)
            if edge_data:
                for _, value in edge_data.items():
                    if value["atom_id_1"] == atom1 and value["atom_id_2"] == atom2:
                        add_edge = False
            if add_edge:
                bms.add_edge(ptnr1, ptnr2, atom_id_1=atom1, atom_id_2=atom2)


def __add_con_branch_link(entity_branch_link, branch_scheme, bms):
    entities = {}
    for i in range(len(branch_scheme["entity_id"])):
        if branch_scheme["entity_id"][i] not in entities:
            entities[branch_scheme["entity_id"][i]] = [branch_scheme["pdb_asym_id"][i]]
        elif (
            branch_scheme["pdb_asym_id"][i]
            not in entities[branch_scheme["entity_id"][i]]
        ):
            entities[branch_scheme["entity_id"][i]].append(
                branch_scheme["pdb_asym_id"][i]
            )

    for i in range(len(entity_branch_link["link_id"])):
        ent_id = entity_branch_link["entity_id"][i]
        for chain in entities[ent_id]:
            for node in bms.nodes:
                prtnr1 = False
                prtnr2 = False
                atom1 = False
                atom2 = False
                if (
                    node.name == entity_branch_link["comp_id_1"][i]
                    and node.chain == chain
                    and node.res_id == entity_branch_link["entity_branch_list_num_1"][i]
                ):
                    prtnr1 = node
                    atom1 = entity_branch_link["atom_id_1"][i]
                elif (
                    node.name == entity_branch_link["comp_id_2"][i]
                    and node.chain == chain
                    and node.res_id == entity_branch_link["entity_branch_list_num_2"][i]
                ):
                    prtnr2 = node
                    atom2 = entity_branch_link["atom_id_2"][i]
            if prtnr1 and prtnr2 and atom1 and atom2:
                add_edge = True
                edge_data = bms.get_edge_data(prtnr1, prtnr2)
                if edge_data:
                    for _, value in edge_data.items():
                        if value["atom_id_1"] == atom1 and value["atom_id_2"] == atom2:
                            add_edge = False
                if add_edge:
                    bms.add_edge(prtnr1, prtnr2, atom_id_1=atom1, atom_id_2=atom2)
                bms.add_edge(prtnr1, prtnr2, atom_id_1=atom1, atom_id_2=atom2)


def parse_bound_molecules(path, to_discard):
    """Parse information from the information about HETATMS from the
    `_pdbx_nonpoly_scheme` and connectivity among them from `_struct_conn`.

    Args:
        path (str): Path to the mmCIF structure
        to_discard (list of str): List of residue names to be discarded.

    Returns:
        MultiDiGraph: All the bound molecules in a given entry.
    """
    doc = cif.read(path)
    cif_block = doc.sole_block()

    if (
        "_pdbx_nonpoly_scheme." not in cif_block.get_mmcif_category_names()
        and "_pdbx_branch_scheme." not in cif_block.get_mmcif_category_names()
    ):
        return MultiDiGraph()

    if "_pdbx_nonpoly_scheme." in cif_block.get_mmcif_category_names():
        bms = parse_ligands_from_nonpoly_scheme(
            cif_block.get_mmcif_category("_pdbx_nonpoly_scheme."), to_discard
        )

    else:
        bms = MultiDiGraph()

    if "_pdbx_branch_scheme." in cif_block.get_mmcif_category_names():
        if bms:
            bms = parse_ligands_from_branch_scheme(
                cif_block.get_mmcif_category("_pdbx_branch_scheme."), to_discard, bms
            )
        else:
            bms = parse_ligands_from_branch_scheme(
                cif_block.get_mmcif_category("_pdbx_branch_scheme."),
                to_discard,
                MultiDiGraph(),
            )

    if "_struct_conn." in cif_block.get_mmcif_category_names():
        __add_connections(cif_block.get_mmcif_category("_struct_conn."), bms)

    if "_pdbx_branch_scheme." in cif_block.get_mmcif_category_names():
        __add_con_branch_link(
            cif_block.get_mmcif_category("_pdbx_entity_branch_link."),
            cif_block.get_mmcif_category("_pdbx_branch_scheme."),
            bms,
        )

    return bms


def find_pntr_entry(struct_conn, residue_pool, partner, i):
    """Helper method to find ligand residue in parsed ligands and check
    its connections.

    Just because it is chatty and used throughout the code quite often.

    Args:
        struct_conn (dict of str: str): struct_conn table.
        residue_pool (List of Residue): List of all the parsed residues.
        partner (str): Identification of the partner (should be 1 or 2)
        i (int): Index in the struct_conn table
    """
    atom_name = struct_conn[f"ptnr{partner}_label_atom_id"][i]
    for residue in residue_pool:
        if (
            residue.name == struct_conn[f"ptnr{partner}_label_comp_id"][i]
            and residue.chain == struct_conn[f"ptnr{partner}_auth_asym_id"][i]
            and residue.res_id == struct_conn[f"ptnr{partner}_auth_seq_id"][i]
        ):
            return (residue, atom_name)

    return (None, atom_name)


def parse_ligands_from_branch_scheme(branch_scheme, to_discard, g):
    """Parse ligands from the mmcif file.

    Args:
        nonpoly_scheme (dict of str: list of str): mmcif _nonpoly_scheme
            category.
        to_discard (list of str): List of residue names to be discarded.

    Returns:
        MultiDiGraph: Ligands and their connectivity in a PDB entry
    """

    for i in range(len(branch_scheme["asym_id"])):
        n = Residue(
            branch_scheme["pdb_mon_id"][i],  # aka label_comp_id
            branch_scheme["pdb_asym_id"][i],  # aka auth_asym_id
            branch_scheme["num"][i],  # aka auth_seq_id
            None,  # aka pdbx_PDB_ins_code
            branch_scheme["entity_id"][i],
        )

        if n.name not in to_discard:
            g.add_node(n)

    return g


def parse_ligands_from_nonpoly_scheme(nonpoly_scheme, to_discard):
    """Parse ligands from the mmcif file.

    Args:
        nonpoly_scheme (dict of str: list of str): mmcif _nonpoly_scheme
            category.
        to_discard (list of str): List of residue names to be discarded.

    Returns:
        MultiDiGraph: Ligands and their connectivity in a PDB entry
    """
    g = MultiDiGraph()

    for i in range(len(nonpoly_scheme["asym_id"])):
        n = Residue(
            nonpoly_scheme["pdb_mon_id"][i],  # aka label_comp_id
            nonpoly_scheme["pdb_strand_id"][i],  # aka auth_asym_id
            nonpoly_scheme["pdb_seq_num"][i],  # aka auth_seq_id
            nonpoly_scheme["pdb_ins_code"][i],  # aka pdbx_PDB_ins_code
            nonpoly_scheme["entity_id"][i],
        )

        if n.name not in to_discard:
            g.add_node(n)

    return g
