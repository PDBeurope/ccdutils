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
A set of methods for identifying bound-molecules (covalently bonded CCDs ) from mmCIF files
of proteins and creating MultiDiGraph representation of the molecules.
"""

from gemmi import cif
from pdbeccdutils.helpers import helper
from networkx import MultiDiGraph, connected_components
from pdbeccdutils.core.models import (
    Residue,
    AssemblyResidue,
    BoundMolecule,
)


def infer_bound_molecules(structure, to_discard, assembly=False):
    """Identify bound molecules in the input protein structure.

    Args:
        structure (str): Path to the structure.
        to_discard (list of str): List of residue names to be discarded
    """

    bms = parse_bound_molecules(structure, to_discard, assembly)
    bound_molecules = []

    for bm_nodes in connected_components(bms.to_undirected()):
        subgraph = bms.subgraph(bm_nodes)
        bm = BoundMolecule(subgraph)
        bound_molecules.append(bm)

    bound_molecules = sorted(bound_molecules, key=lambda l: -len(l.graph.nodes))
    return bound_molecules


def __add_connections(struct_conn: dict[str, list[str]], bms: MultiDiGraph) -> None:
    """Adds covalent connections between chemical components in a
    bound-molecule from `_struct_conn` table

    Args:
        struct_conn: Dictionary of `_struct_conn` category from mmCIF
        bms: Graph corresponding to bound-molecule
    """
    for i in range(len(struct_conn["id"])):
        (ptnr1, atom1) = find_pntr_entry(struct_conn, bms.nodes, 1, i)
        (ptnr2, atom2) = find_pntr_entry(struct_conn, bms.nodes, 2, i)

        # we want covalent connections among ligands only.
        if ptnr1 and ptnr2 and struct_conn["conn_type_id"][i] != "metalc":
            add_edge = True
            # In MultiDiGraph two nodes can have multiple edges even if the edges share the same edge data.
            # _struct_conn has duplicate links between atoms, they are removed by the below process.
            edge_data = bms.get_edge_data(ptnr1, ptnr2)
            if edge_data:
                for _, value in edge_data.items():
                    if value["atom_id_1"] == atom1 and value["atom_id_2"] == atom2:
                        add_edge = False
            if add_edge:
                bms.add_edge(ptnr1, ptnr2, atom_id_1=atom1, atom_id_2=atom2)


def __add_con_branch_link(
    entity_branch_link: dict[str, list[str]],
    branch_scheme: dict[str, list[str]],
    bms: MultiDiGraph,
) -> None:
    """Adds connections between chemical components in a bound-molecule from
    `_pdbx_entity_branch_link` table

    Args:
        entity_branch_link: Dictionary of `_pdbx_entity_branch_link` category
        branch_scheme: Dictionary of `_pdbx_branch_scheme` category
    """
    entities = {}
    for i in range(len(branch_scheme["entity_id"])):
        if not entities.get(branch_scheme["entity_id"][i]):
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


def parse_bound_molecules(
    path: str, to_discard: list[str], assembly=False
) -> MultiDiGraph:
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
            cif_block.get_mmcif_category("_pdbx_nonpoly_scheme."), to_discard, assembly
        )

    else:
        bms = MultiDiGraph()

    if "_pdbx_branch_scheme." in cif_block.get_mmcif_category_names():
        if bms:
            bms = parse_ligands_from_branch_scheme(
                cif_block.get_mmcif_category("_pdbx_branch_scheme."),
                to_discard,
                bms,
                assembly,
            )
        else:
            bms = parse_ligands_from_branch_scheme(
                cif_block.get_mmcif_category("_pdbx_branch_scheme."),
                to_discard,
                MultiDiGraph(),
                assembly,
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


def find_pntr_entry(
    struct_conn: dict[str, list[str]], residue_pool: list[Residue], partner: int, i: int
):
    """Helper method to find ligand residue in parsed ligands and check
    its connections.

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


def parse_ligands_from_branch_scheme(
    branch_scheme: dict[str, list[str]],
    to_discard: list[str],
    g: MultiDiGraph,
    assembly=False,
):
    """Parse ligands from `_pdbx_branch_scheme` category of mmCIF file

    Args:
        branch_scheme: Dictionary of `_pdbx_branch_scheme` category
        to_discard: List of residue names to be not considered as bound-molecule
        g: A Graph object with nodes as Resiudes and their connectivity as edges

    Returns:
        A MultiDiGraph object with nodes as Resiudes and their connectivity a
    """

    if assembly:
        for i in range(len(branch_scheme["asym_id"])):
            (orig_auth_asym_id, operator) = helper.get_additional_fields(
                branch_scheme["pdb_asym_id"][i]
            )
            n = AssemblyResidue(
                branch_scheme["pdb_mon_id"][i],  # aka label_comp_id
                branch_scheme["pdb_asym_id"][i],  # aka auth_asym_id
                branch_scheme["num"][i],  # aka auth_seq_id
                None,  # aka pdbx_PDB_ins_code
                branch_scheme["entity_id"][i],
                orig_auth_asym_id,
                operator,
            )

            if n.name not in to_discard:
                g.add_node(n)
    else:
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


def parse_ligands_from_nonpoly_scheme(nonpoly_scheme, to_discard, assembly=False):
    """Parse ligands from the mmcif file.

    Args:
        nonpoly_scheme (dict of str: list of str): mmcif _nonpoly_scheme
            category.
        to_discard (list of str): List of residue names to be discarded.

    Returns:
        MultiDiGraph: Ligands and their connectivity in a PDB entry
    """
    g = MultiDiGraph()

    if assembly:
        for i in range(len(nonpoly_scheme["asym_id"])):
            (orig_auth_asym_id, operator) = helper.get_additional_fields(
                nonpoly_scheme["pdb_strand_id"][i]
            )
            n = AssemblyResidue(
                nonpoly_scheme["pdb_mon_id"][i],  # aka label_comp_id
                nonpoly_scheme["pdb_strand_id"][i],  # aka auth_asym_id
                nonpoly_scheme["pdb_seq_num"][i],  # aka auth_seq_id
                nonpoly_scheme["pdb_ins_code"][i],  # aka pdbx_PDB_ins_code
                nonpoly_scheme["entity_id"][i],
                orig_auth_asym_id,
                operator,
            )

            if n.name not in to_discard:
                g.add_node(n)
    else:
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
