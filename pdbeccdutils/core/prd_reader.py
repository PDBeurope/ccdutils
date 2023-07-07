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


import logging
import os
from typing import Dict
import rdkit
from pdbeccdutils.core.component import Component
from pdbeccdutils.core import ccd_reader
from pdbeccdutils.core.exceptions import CCDUtilsError
from pdbeccdutils.core.models import (
    CCDProperties,
    ReleaseStatus,
)
from pdbeccdutils.helpers import cif_tools, conversions, mol_tools
from gemmi import cif

# categories that need to be 'fixed'
# str => list[str]
preprocessable_categories = [
    "_chem_comp_atom.",
    "_chem_comp_bond.",
    "_pdbx_chem_comp_identifier.",
    "_pdbx_chem_comp_descriptor.",
    "_chem_comp.",
]


def read_pdb_cif_file(
    path_to_cif: str, sanitize: bool = True
) -> ccd_reader.CCDReaderResult:
    """
    Read in single wwPDB CCD CIF component and create its internal
    representation.

    Args:
        path_to_cif (str): Path to the cif file
        sanitize (bool): [Defaults: True]

    Raises:
        ValueError: if file does not exist

    Returns:
        CCDReaderResult: Results of the parsing altogether with
        the internal representation of the component.
    """
    if not os.path.isfile(path_to_cif):
        raise ValueError("File '{}' does not exists".format(path_to_cif))

    doc = cif.read(path_to_cif)
    cif_block = doc.sole_block()

    return _parse_pdb_mmcif(cif_block, sanitize)


def read_pdb_components_file(
    path_to_cif: str, sanitize: bool = True
) -> Dict[str, ccd_reader.CCDReaderResult]:
    """
    Process multiple compounds stored in the wwPDB CCD
    `components.cif` file.

    Args:
        path_to_cif (str): Path to the `prdcc-all.cif` file with
            multiple ligands in it.
        sanitize (bool): Whether or not the components should be sanitized
            Defaults to True.

    Raises:
        ValueError: if the file does not exist.

    Returns:
        dict[str, CCDReaderResult]: Internal representation of all
        the components in the `components.cif` file.
    """
    if not os.path.isfile(path_to_cif):
        raise ValueError("File '{}' does not exists".format(path_to_cif))

    result_bag = {}

    for block in cif.read(path_to_cif):
        try:
            result_bag[block.name] = _parse_pdb_mmcif(block, sanitize)
        except CCDUtilsError as e:
            logging.error(
                f"ERROR: Data block {block.name} not processed. Reason: ({str(e)})."
            )

    return result_bag


# region parse mmcif


def _parse_pdb_mmcif(cif_block, sanitize=True):
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
    sanitized = False
    mol = rdkit.Chem.RWMol()

    for c in preprocessable_categories:
        w = cif_tools.preprocess_cif_category(cif_block, c)

        if w:
            warnings.append(w)

    _parse_pdb_atoms(mol, cif_block)
    ccd_reader._parse_pdb_conformers(mol, cif_block)
    ccd_reader._parse_pdb_bonds(mol, cif_block, errors)
    ccd_reader._handle_implicit_hydrogens(mol)

    if sanitize:
        sanitized = mol_tools.sanitize(mol)

    descriptors = ccd_reader._parse_pdb_descriptors(
        cif_block, "_pdbx_chem_comp_descriptor.", "descriptor"
    )
    descriptors += ccd_reader._parse_pdb_descriptors(
        cif_block, "_pdbx_chem_comp_identifier.", "identifier"
    )
    properties = _parse_pdb_properties(cif_block)

    comp = Component(mol.GetMol(), cif_block, properties, descriptors)
    reader_result = ccd_reader.CCDReaderResult(
        warnings=warnings, errors=errors, component=comp, sanitized=sanitized
    )

    return reader_result


def _parse_pdb_atoms(mol, cif_block):
    """
    Setup atoms in the component

    Args:
        mol (rdkit.Chem.rdchem.Mol): Rdkit Mol object with the
            compound representation.
        cif_block (cif.Block): mmCIF block object from gemmi.
    """
    if "_chem_comp_atom." not in cif_block.get_mmcif_category_names():
        return

    atoms = cif_block.find(
        "_chem_comp_atom.",
        [
            "atom_id",
            "type_symbol",
            "alt_atom_id",
            "pdbx_leaving_atom_flag",
            "charge",
            "pdbx_component_comp_id",
            "pdbx_residue_numbering",
            "pdbx_component_atom_id",
            "pdbx_polymer_type",
            "pdbx_ref_id",
            "pdbx_component_id",
        ],
    )
    for row in atoms:
        atom_id = cif.as_string(row["_chem_comp_atom.atom_id"])
        element = cif.as_string(row["_chem_comp_atom.type_symbol"])
        alt_atom_id = cif.as_string(row["_chem_comp_atom.alt_atom_id"])
        leaving_atom = cif.as_string(row["_chem_comp_atom.pdbx_leaving_atom_flag"])
        charge = cif.as_string(row["_chem_comp_atom.charge"])
        comp_atom_id = cif.as_string(row["_chem_comp_atom.pdbx_component_atom_id"])
        res_name = cif.as_string(row["_chem_comp_atom.pdbx_component_comp_id"])
        residue_id = cif.as_string(row["_chem_comp_atom.pdbx_residue_numbering"])
        res_type = cif.as_string(row["_chem_comp_atom.pdbx_polymer_type"])
        ref_id = cif.as_string(row["_chem_comp_atom.pdbx_ref_id"])
        comp_id = cif.as_string(row["_chem_comp_atom.pdbx_component_id"])

        element = element if len(element) == 1 else element[0] + element[1].lower()
        isotope = None

        if element == "D":
            element = "H"
            isotope = 2
        elif element == "X":
            element = "*"

        atom = rdkit.Chem.Atom(element)
        atom.SetProp("name", atom_id)
        atom.SetProp("alt_name", alt_atom_id)
        atom.SetBoolProp("leaving_atom", leaving_atom == "Y")
        atom.SetProp("component_atom_id", comp_atom_id)
        atom.SetProp("ref_id", ref_id)
        atom.SetProp("comp_id", comp_id)
        atom.SetProp("res_type", res_type)
        atom.SetProp("residue_id", residue_id)
        atom.SetFormalCharge(conversions.str_to_int(charge))

        res_info = rdkit.Chem.AtomPDBResidueInfo()
        res_info.SetResidueName(res_name)
        res_info.SetIsHeteroAtom(True)
        atom.SetMonomerInfo(res_info)

        if isotope is not None:
            atom.SetIsotope(isotope)

        mol.AddAtom(atom)


def _parse_pdb_properties(cif_block):
    """Parse useful information from _chem_comp category

    Args:
        cif_block (cif.Block): mmcif block object from gemmi

    Returns:
        Properties: dataclass with the CCD properties.
    """
    properties = None
    if "_chem_comp." in cif_block.get_mmcif_category_names():
        rel_status = ReleaseStatus.from_str(
            cif_block.find_value("_chem_comp.pdbx_release_status")
        )
        formula_weight = cif_block.find_value("_chem_comp.formula_weight")
        weight = 0.0 if cif.is_null(formula_weight) else cif.as_number(formula_weight)

        properties = CCDProperties(
            id=cif_tools.get_prd_cc_code(
                cif.as_string(cif_block.find_value("_chem_comp.id"))
            ),
            name=cif.as_string(cif_block.find_value("_chem_comp.name")),
            formula=cif.as_string(cif_block.find_value("_chem_comp.formula")),
            modified_date="",
            pdbx_release_status=rel_status,
            weight=weight,
        )
    return properties


# endregion parse mmcif
