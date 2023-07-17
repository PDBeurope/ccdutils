#!/usr/bin/env python
# software from PDBe: Protein Data Bank in Europe; https://pdbe.org
#
# Copyright 2021 EMBL - European Bioinformatics Institute
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
Set of methods to format data for gemmi parser
"""

import gemmi
from gemmi import cif


def preprocess_cif_category(cif_block, label):
    """
    Checks if the category is present in gemmi.cif.Block object

    Args:
        cif_block (Block): mmcif Block from gemmi.
        label (str): name of the category

    Returns:
        str: Possible error encountered
    """
    if label not in cif_block.get_mmcif_category_names():
        return f"Namespace {label} does not exist."


def fix_updated_mmcif(input_str: str, output_str: str) -> None:
    """Keeps only first model, remove alternate conformations of atoms and residues.
    Updates _pdbx_branch_scheme, _struct_conn, _pdbx_nonpoly_scheme and _atom_site

    Args:
        input_str : Path to input structure
        output_str: Path to write processed structure

    """

    cif_block = gemmi.cif.read(input_str).sole_block()
    st = gemmi.make_structure_from_block(cif_block)
    # remove all other models except one
    del st[1:]

    # remove alternate conformations of atoms and residues
    st.remove_alternative_conformations()
    groups = gemmi.MmcifOutputGroups(True)
    groups.group_pdb = True
    st_cif_block = st.make_mmcif_block(groups)

    cif_block.set_mmcif_category(
        "_atom_site.", st_cif_block.get_mmcif_category("_atom_site."), raw=True
    )
    cif_block_categories = cif_block.get_mmcif_category_names()
    if "_struct_conn." in cif_block_categories:
        cif_block.set_mmcif_category(
            "_struct_conn.", st_cif_block.get_mmcif_category("_struct_conn."), raw=True
        )

    # keep only first conformer of residues in _pdbx_branch_scheme
    if "_pdbx_branch_scheme." in cif_block_categories:
        new_pdbx_branch_scheme = _filter_pdbx_branch_scheme(cif_block, st)
        cif_block.set_mmcif_category(
            "_pdbx_branch_scheme.", new_pdbx_branch_scheme, raw=False
        )
        # Remove links of residues not in _pdbx_branch_scheme from _pdbx_entity_branch_link
        if "_pdbx_entity_branch_link." in cif_block_categories:
            new_entity_branch_link = _filter_pdbx_entity_branch_link(
                cif_block, new_pdbx_branch_scheme
            )
            cif_block.set_mmcif_category(
                "_pdbx_entity_branch_link.", new_entity_branch_link, raw=False
            )

    # keep only first conformer of residues in _pdbx_nonpoly_scheme
    if "_pdbx_nonpoly_scheme." in cif_block_categories:
        new_nonpoly_scheme = _filter_pdbx_nonpoly_scheme(cif_block, st)
        cif_block.set_mmcif_category(
            "_pdbx_nonpoly_scheme.", new_nonpoly_scheme, raw=False
        )

    cif_block.write_file(output_str, cif.Style.Pdbx)


def _filter_pdbx_branch_scheme(
    cif_block: gemmi.cif.Block, structure: gemmi.Structure
) -> dict[str, list[str]]:
    """Removes residues from _pdbx_branch_scheme that are not in _atom_site

    Args:
        cif_block: Parsed Block object from gemmi of input structure
        structure: Structure object from gemmi of input structure
            after removing alternative conformations

    Returns:
        Updated dictionary of '_pdbx_branch_scheme' category

    """
    branch_scheme = cif_block.get_mmcif_category("_pdbx_branch_scheme.")
    branch_entity_ids = set(branch_scheme["entity_id"])
    branch_chains = set(branch_scheme["pdb_asym_id"])
    # chain number from gemmi maps to _atom_site.auth_asym_id
    st_residues = [
        (res.entity_id, chain, res.name, str(res.seqid.num))
        for chain in branch_chains
        for res in structure[0][chain]
        if res.entity_id in branch_entity_ids
    ]
    new_branch_scheme = {key: [] for key in branch_scheme}
    for i in range(len(branch_scheme["entity_id"])):
        branch_scheme_residue = (
            branch_scheme["entity_id"][i],
            branch_scheme["pdb_asym_id"][
                i
            ],  # _pdbx_branch_scheme.pdb_asym_id maps to _atom_site.auth_asym_id
            branch_scheme["pdb_mon_id"][i],
            branch_scheme["num"][i],
        )
        if branch_scheme_residue in st_residues:
            for key in new_branch_scheme:
                new_branch_scheme[key].append(branch_scheme[key][i])

    return new_branch_scheme


def _filter_pdbx_entity_branch_link(
    cif_block: gemmi.cif.Block, branch_scheme: dict
) -> dict[str, list[str]]:
    """Removes residues from `_pdbx_entity_branch_link` that are not in _pdbx_branch_scheme

    Args:
        cif_block: Parsed Block object from gemmi of input structure
        branch_scheme: Dictionary of `_pdbx_branch_scheme` of input structure
            after removing alternative conformations

    Returns:
        Updated dictionary of '_pdbx_entity_branch_link' category

    """
    branch_scheme_residues = [
        (
            branch_scheme["entity_id"][i],
            branch_scheme["mon_id"][i],
            branch_scheme["num"][i],
        )
        for i in range(len(branch_scheme["entity_id"]))
    ]

    entity_branch_link = cif_block.get_mmcif_category("_pdbx_entity_branch_link.")

    new_entity_branch_link = {key: [] for key in entity_branch_link}

    for i in range(len(entity_branch_link["entity_id"])):
        comp_1 = (
            entity_branch_link["entity_id"][i],
            entity_branch_link["comp_id_1"][i],
            entity_branch_link["entity_branch_list_num_1"][i],
        )
        comp_2 = (
            entity_branch_link["entity_id"][i],
            entity_branch_link["comp_id_2"][i],
            entity_branch_link["entity_branch_list_num_2"][i],
        )
        if (comp_1 in branch_scheme_residues) and (comp_2 in branch_scheme_residues):
            for key in new_entity_branch_link:
                new_entity_branch_link[key].append(entity_branch_link[key][i])

    return new_entity_branch_link


def _filter_pdbx_nonpoly_scheme(
    cif_block: gemmi.cif.Block, structure: gemmi.Structure
) -> dict[str, list[str]]:
    """Removes residues from _pdbx_nonpoly_scheme that are not in _atom_site
    Args:
        cif_block: Parsed Block object from gemmi of input structure
        structure: Structure object from gemmi of input structure
            after removing alternative conformations

    Returns:
        Updated dictionary of '_pdbx_nonpoly_scheme' category

    """
    nonpoly_scheme = cif_block.get_mmcif_category("_pdbx_nonpoly_scheme.")
    nonpoly_entity_ids = set(nonpoly_scheme["entity_id"])
    # subchain number from gemmi maps to _atom_site.label_asym_id
    st_residues = [
        (res.entity_id, res.subchain, res.name, str(res.seqid.num))
        for chain in structure[0]
        for res in chain
        if res.entity_id in nonpoly_entity_ids
    ]
    new_nonpoly_scheme = {key: [] for key in nonpoly_scheme}
    for i in range(len(nonpoly_scheme["asym_id"])):
        nonpoly_scheme_residue = (
            nonpoly_scheme["entity_id"][i],
            nonpoly_scheme["asym_id"][
                i
            ],  # _pdbx_nonpoly_scheme.asym_id maps to _atom_site.label_asym_id
            nonpoly_scheme["pdb_mon_id"][i],
            nonpoly_scheme["pdb_seq_num"][i],
        )
        if nonpoly_scheme_residue in st_residues:
            for key in new_nonpoly_scheme:
                new_nonpoly_scheme[key].append(nonpoly_scheme[key][i])

    return new_nonpoly_scheme


def get_prd_cc_code(prd_code: str):
    """Returns PRDCC code from PRD code

    Args:
        prd_code: ID of PRD
    """
    prefix, code = prd_code.split("_")
    prdcc_code = f"{prefix}CC_{code}"
    return prdcc_code
