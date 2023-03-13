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

"""
A set of methods for reading in data of PDB model cif files and creating an array of
internal representation of bound molecules. The basic use can be as easy as this:

    from pdbeccdutils.core import bm_reader

    bound_molecules = bm_reader.read_pdb_updated_cif_file('/path/to/cif/xxx_updated.cif')
    for i in bound_molecules:
        rdkit_mol = i.component.mol
        rdkit_inchikey = i.component.inchikey_from_rdkit
"""

import os
import rdkit
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
from pdbeccdutils.core import ccd_reader
from pdbeccdutils.core.component import Component

from pdbeccdutils.core.models import (
    CCDProperties,
    ConformerType,
    Descriptor,
    Residue,
    BoundMolecule,
)
from pdbeccdutils.helpers import cif_tools, conversions, mol_tools, helper
from gemmi import cif
from networkx import DiGraph, connected_components


def read_pdb_updated_cif_file(path_to_cif: str, sanitize: bool = True):
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
    bms = infer_bound_molecules(path_to_cif, ["HOH"])
    for bm in bms:
        reader_result = infer_multiple_chem_comp(path_to_cif, bm.to_dict(), sanitize)
        if reader_result:
            biomolecule_result.append(reader_result)

    return biomolecule_result


def infer_multiple_chem_comp(path_to_cif: str, bm: dict, sanitize: bool = True):
    """
    Read in single wwPDB Model CIF and single bound-molecule to create internal
    representation of the bound molecule.

    Args:
        path_to_cif (str): Path to the cif file
        bm (dict): Dictionary representation of bound-molecule
        sanitize (bool): [Defaults: True]


    Returns:
        A dictionary of CCDResult representation of bound-molecule.
    """

    if len(bm["residues"]) <= 1:
        return

    else:
        warnings = []
        errors = []
        sanitized = False
        cif_block = cif.read(path_to_cif).sole_block()
        cif_tools.preprocess_cif_category(cif_block, "_atom_site.")
        cif_tools.preprocess_cif_category(cif_block, "_chem_comp_bond.")

        mol = rdkit.Chem.RWMol()
        index_atoms, bm_atoms_dict = _parse_pdb_atom_site(
            mol, cif_block.get_mmcif_category("_atom_site."), bm
        )
        _parse_pdb_conformers_site(mol, bm_atoms_dict, index_atoms)
        _parse_pdb_bonds_site(
            mol,
            cif_block.get_mmcif_category("_chem_comp_bond."),
            bm_atoms_dict,
            errors,
            index_atoms,
            bm,
        )

        _handle_disconnected_hydrogens(mol)
        if sanitize:
            sanitized = mol_tools.sanitize(mol)

        comp = Component(mol.GetMol(), cif_block)
        descriptors = [
            Descriptor(type="InChI", program="rdkit", value=comp.inchi_from_rdkit),
            Descriptor(
                type="InChIKey", program="rdkit", value=comp.inchikey_from_rdkit
            ),
        ]
        properties = CCDProperties(
            id="",
            name="",
            formula=CalcMolFormula(comp.mol),
            modified_date="",
            pdbx_release_status="",
            weight=round(comp.physchem_properties["exactmw"], 3),
        )
        comp = Component(mol.GetMol(), cif_block, properties, descriptors)
        reader_result = ccd_reader.CCDReaderResult(
            warnings=warnings, errors=errors, component=comp, sanitized=sanitized
        )

        return reader_result


def _handle_disconnected_hydrogens(mol):
    """
    Delete hydrogens without neighbours (degree = 0).
    RDKit works but gives warning for these ("WARNING: not removing hydrogen atom without neighbors").
    However we wouldn't want these hydrogens as they affect molecular properties.

    Args:
        mol (rdkit.Chem.rchem.Mol): Rdkit Mol object with the
            compound representation.
    """
    disconnected_hydrogens = [
        atom
        for atom in mol.GetAtoms()
        if atom.GetAtomicNum() == 1 and atom.GetDegree() == 0
    ]
    atoms_to_remove = [atom.GetIdx() for atom in disconnected_hydrogens]
    atoms_to_remove.sort(reverse=True)

    for atom in atoms_to_remove:
        mol.RemoveAtom(atom)


def _parse_pdb_atom_site(mol, atoms, bm):
    """
    Setup atoms in the bound molecule

    Args:
        mol (rdkit.Chem.rchem.Mol): RDkit Mol object with the
            compound representation.
        atoms (dict): dictionary with parsed _chem_comp_atom
            category.
        bm (dict): dictionary representation of bound-molecule
    """
    if not atoms:
        return

    exclude_alternate_atoms = True
    if exclude_alternate_atoms:
        seen = {}

    index_atoms = []
    bm_atoms_dict = {}

    for i in range(len(atoms["id"])):
        if atoms["group_PDB"][i] != "HETATM":  # only consider nonpolys
            continue
        for j in bm["residues"]:
            if (
                atoms["label_comp_id"][i] == j["label_comp_id"]
                and atoms["auth_asym_id"][i] == j["auth_asym_id"]
                and atoms["auth_seq_id"][i] == j["auth_seq_id"]
            ):
                atom_tuple = (
                    atoms["auth_seq_id"][i],
                    atoms["auth_atom_id"][i],
                    atoms["auth_comp_id"][i],
                )
                if exclude_alternate_atoms:
                    if atom_tuple in seen:
                        continue
                    seen[atom_tuple] = True

                for categories in atoms:
                    if categories not in bm_atoms_dict:
                        bm_atoms_dict[categories] = []
                    else:
                        bm_atoms_dict[categories].append(atoms[categories][i])

                element = atoms["type_symbol"][i]
                element = (
                    element if len(element) == 1 else element[0] + element[1].lower()
                )
                isotope = None

                if element == "D":
                    element = "H"
                    isotope = 2
                elif element == "X":
                    element = "*"

                atom = rdkit.Chem.Atom(element)
                atom.SetProp("name", "/".join(atom_tuple))
                index_atoms.append(
                    (
                        f"{atoms['auth_asym_id'][i]}/{atoms['auth_seq_id'][i]}",
                        atoms["auth_comp_id"][i],
                        atoms["auth_atom_id"][i],
                    )
                )

                if isotope is not None:
                    atom.SetIsotope(isotope)

                mol.AddAtom(atom)

    return index_atoms, bm_atoms_dict


def _parse_pdb_conformers_site(mol, atoms, index_atoms):
    """
    Setup model cooordinates in the rdkit Mol object.

    Args:
        mol (rdkit.Chem.rdchem.Mol): RDKit Mol object with the compound
            representation.
        atoms (dict): mmcif category with atom info category.
        index_atoms (list): List of intx atoms
    """
    if not atoms:
        return

    model = _setup_pdb_conformer_site(
        atoms, "Cartn_{}", ConformerType.Model.name, index_atoms
    )

    mol.AddConformer(model, assignId=True)


def _setup_pdb_conformer_site(atoms, label, name, index_atoms):
    """
    Setup a conformer

    Args:
        atoms (dict): mmcif category with the atom info.
        label (str): Namespace with the [x,y,z] coordinates.
        name (str): Conformer name.
        index_atoms (lsit): List of index atoms

    Returns:
        rdkit.Chem.rdchem.Conformer: Conformer of the component.
    """
    if not atoms:
        return

    conformer = rdkit.Chem.Conformer(len(index_atoms))

    j = 0
    for i in range(len(atoms["id"])):
        if atoms["group_PDB"][i] != "HETATM":
            continue
        x = conversions.str_to_float(atoms[label.format(("x"))][i])
        y = conversions.str_to_float(atoms[label.format(("y"))][i])
        z = conversions.str_to_float(atoms[label.format(("z"))][i])
        index = helper.find_element_in_list(
            index_atoms,
            (
                f"{atoms['auth_asym_id'][i]}/{atoms['auth_seq_id'][i]}",
                atoms["auth_comp_id"][i],
                atoms["auth_atom_id"][i],
            ),
        )

        atom_position = rdkit.Chem.rdGeometry.Point3D(x, y, z)
        conformer.SetAtomPosition(index, atom_position)
        j = j + 1

        conformer.SetProp("name", name)

    return conformer


def _parse_pdb_bonds_site(mol, bonds, atoms, errors, index_atoms, bm):
    """
    Setup bonds in the compound

    Args:
        mol (rdkit.Chem.rdchem.Mol): Molecule which receives bonds.
        bonds (dict): mmcif category with the bonds info. TODO not being used
        atoms (dict): mmcif category with the atom info. TODO not being used
        errors (list[str]): Issues encountered while parsing.
        index_atoms (list): List of index atoms
        bm (dict): Dictionary representation of bound-molecule
    """
    if not bonds:
        return

    bond_per_type = {}
    for i in range(len(bonds["atom_id_1"])):
        content = (
            bonds["atom_id_1"][i],
            bonds["atom_id_2"][i],
            bonds["value_order"][i],
        )
        if bonds["comp_id"][i] in bond_per_type:
            bond_per_type[bonds["comp_id"][i]].append(content)
        else:
            bond_per_type[bonds["comp_id"][i]] = [content]

    # get connection mapping
    connection_mapping = {}
    for j in bm["residues"]:
        connection_mapping[j["id"]] = j

    for j in bm["residues"]:
        lig = j["label_comp_id"]
        if lig not in bond_per_type:  # ions
            continue
        chain_res = f"{j['auth_asym_id']}/{j['auth_seq_id']}"
        try:
            for pairs in bond_per_type[lig]:
                tuple_to_find_1 = (chain_res, lig, pairs[0])
                tuple_to_find_2 = (chain_res, lig, pairs[1])
                atom_1 = helper.find_element_in_list(index_atoms, tuple_to_find_1)
                atom_2 = helper.find_element_in_list(index_atoms, tuple_to_find_2)
                bond_order = ccd_reader._bond_pdb_order(pairs[2])
                if any(a is None for a in [atom_1, atom_2, bond_order]):
                    pass
                else:
                    mol.AddBond(atom_1, atom_2, bond_order)
        except ValueError:
            errors.append(
                f"Error perceiving {atom_1} - {atom_2} bond in _chem_comp_bond"
            )
        except RuntimeError:
            errors.append(f"Duplicit bond {atom_1} - {atom_2}")

    # Add bound molecule connections
    try:
        for connection in bm["connections"]:
            first_res = connection_mapping[connection[0]]
            second_res = connection_mapping[connection[1]]
            chain_res_1 = f"{first_res['auth_asym_id']}/{first_res['auth_seq_id']}"
            chain_res_2 = f"{second_res['auth_asym_id']}/{second_res['auth_seq_id']}"
            lig_1 = first_res["label_comp_id"]
            lig_2 = second_res["label_comp_id"]
            tuple_to_find_1 = (chain_res_1, lig_1, connection[2]["atom_id_1"])
            tuple_to_find_2 = (chain_res_2, lig_2, connection[2]["atom_id_2"])
            atom_1 = helper.find_element_in_list(index_atoms, tuple_to_find_1)
            atom_2 = helper.find_element_in_list(index_atoms, tuple_to_find_2)

            if any(a is None for a in [atom_1, atom_2]):
                pass
            else:
                mol.AddBond(atom_1, atom_2, ccd_reader._bond_pdb_order("SING"))
    except ValueError:
        errors.append(f"Error perceiving {atom_1} - {atom_2} bond in _struct_conn")
    except RuntimeError:
        errors.append(f"Duplicit bond {atom_1} - {atom_2}")


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
                    atom1 = node.name
                elif (
                    node.name == entity_branch_link["comp_id_2"][i]
                    and node.chain == chain
                    and node.res_id == entity_branch_link["entity_branch_list_num_2"][i]
                ):
                    prtnr2 = node
                    atom2 = node.name
            if prtnr1 and prtnr2:
                bms.add_edge(prtnr1, prtnr2, atom_id_1=atom1, atom_id_2=atom2)


def parse_bound_molecules(path, to_discard):
    """Parse information from the information about HETATMS from the
    `_pdbx_nonpoly_scheme` and connectivity among them from `_struct_conn`.

    Args:
        path (str): Path to the mmCIF structure
        to_discard (list of str): List of residue names to be discarded.

    Returns:
        DiGraph: All the bound molecules in a given entry.
    """
    doc = cif.read(path)
    cif_block = doc.sole_block()

    if (
        "_pdbx_nonpoly_scheme." not in cif_block.get_mmcif_category_names()
        and "_pdbx_branch_scheme." not in cif_block.get_mmcif_category_names()
    ):
        return DiGraph()

    if "_pdbx_nonpoly_scheme." in cif_block.get_mmcif_category_names():
        bms = parse_ligands_from_nonpoly_scheme(
            cif_block.get_mmcif_category("_pdbx_nonpoly_scheme."), to_discard
        )

    else:
        bms = DiGraph()

    if "_pdbx_branch_scheme." in cif_block.get_mmcif_category_names():
        if bms:
            bms = parse_ligands_from_branch_scheme(
                cif_block.get_mmcif_category("_pdbx_branch_scheme."), to_discard, bms
            )
        else:
            bms = parse_ligands_from_branch_scheme(
                cif_block.get_mmcif_category("_pdbx_branch_scheme."),
                to_discard,
                DiGraph(),
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
        DiGraph: Ligands and their connectivity in a PDB entry
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
        DiGraph: Ligands and their connectivity in a PDB entry
    """
    g = DiGraph()

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
