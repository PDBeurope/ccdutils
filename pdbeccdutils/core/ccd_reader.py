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
A set of methods for reading in data and creating internal representation
of molecules. The basic use can be as easy as this:

    from pdbeccdutils.core import ccd_reader

    ccdutils_component = ccd_reader.read_pdb_cif_file('/path/to/cif/ATP.cif').component
    rdkit_mol = ccdutils_component.mol
"""

import logging
import os
from datetime import date
from typing import Dict, List, NamedTuple

import rdkit
from pdbeccdutils.core.component import Component
from pdbeccdutils.core.exceptions import CCDUtilsError
from pdbeccdutils.core.models import (
    CCDProperties,
    ConformerType,
    Descriptor,
    ReleaseStatus,
)
from pdbeccdutils.helpers import cif_tools, conversions, mol_tools, helper
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


class CCDReaderResult(NamedTuple):
    """
    NamedTuple for the result of reading an individual PDB chemical
    component definition (CCD).

    Attributes:
        component (Component): internal representation of the CCD read-in.
        errors (list[str]): A list of any errors
            found while reading the CCD. If no warnings found `errors`
            will be empty.
        warnings (list[str]): A list of any warnings
            found while reading the CCD. If no warnings found `warnings`
            will be empty.
        sanitized (bool): Whether or not the molecule was sanitized
    """

    warnings: List[str]
    errors: List[str]
    component: Component
    sanitized: bool


def read_pdb_cif_file(path_to_cif: str, sanitize: bool = True) -> CCDReaderResult:
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
) -> Dict[str, CCDReaderResult]:
    """
    Process multiple compounds stored in the wwPDB CCD
    `components.cif` file.

    Args:
        path_to_cif (str): Path to the `components.cif` file with
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
    _parse_pdb_conformers(mol, cif_block)
    _parse_pdb_bonds(mol, cif_block, errors)
    _handle_implicit_hydrogens(mol)

    if sanitize:
        sanitized = mol_tools.sanitize(mol)

    descriptors = _parse_pdb_descriptors(
        cif_block, "_pdbx_chem_comp_descriptor.", "descriptor"
    )
    descriptors += _parse_pdb_descriptors(
        cif_block, "_pdbx_chem_comp_identifier.", "identifier"
    )
    properties = _parse_pdb_properties(cif_block)

    comp = Component(mol.GetMol(), cif_block, properties, descriptors)
    reader_result = CCDReaderResult(
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
        ["atom_id", "type_symbol", "alt_atom_id", "pdbx_leaving_atom_flag", "charge"],
    )
    for row in atoms:
        atom_id = cif.as_string(row["_chem_comp_atom.atom_id"])
        element = cif.as_string(row["_chem_comp_atom.type_symbol"])
        alt_atom_id = cif.as_string(row["_chem_comp_atom.alt_atom_id"])
        leaving_atom = cif.as_string(row["_chem_comp_atom.pdbx_leaving_atom_flag"])
        charge = cif.as_string(row["_chem_comp_atom.charge"])

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
        atom.SetFormalCharge(conversions.str_to_int(charge))

        if isotope is not None:
            atom.SetIsotope(isotope)

        mol.AddAtom(atom)


def _parse_pdb_conformers(mol, cif_block):
    """Setup model and ideal cooordinates in the rdkit Mol object.

    Args:
        mol (rdkit.Chem.rdchem.Mol): RDKit Mol object with the compound
            representation.
        cif_block (cif.Block): mmCIF block object from gemmi.
    """

    if "_chem_comp_atom." not in cif_block.get_mmcif_category_names():
        return

    required_fields = [
        "model_Cartn_x",
        "model_Cartn_y",
        "model_Cartn_z",
        "pdbx_model_Cartn_x_ideal",
        "pdbx_model_Cartn_y_ideal",
        "pdbx_model_Cartn_z_ideal",
    ]
    atoms = cif_block.find("_chem_comp_atom.", required_fields)
    ideal = _setup_pdb_conformer(
        atoms, "pdbx_model_Cartn_{}_ideal", ConformerType.Ideal.name
    )
    model = _setup_pdb_conformer(atoms, "model_Cartn_{}", ConformerType.Model.name)

    mol.AddConformer(ideal, assignId=True)
    mol.AddConformer(model, assignId=True)


def _setup_pdb_conformer(atoms, label, name):
    """Setup a conformer

    Args:
        atoms (Table): mmcif table with the atom info.
        label (str): Conformer category
        name (str): Conformer name.

    Returns:
        rdkit.Chem.rdchem.Conformer: Conformer of the component.
    """
    if not atoms:
        return

    conformer = rdkit.Chem.Conformer(len(atoms))

    for row in atoms:
        x = conversions.str_to_float(
            row["_chem_comp_atom.{}".format(label.format("x"))]
        )
        y = conversions.str_to_float(
            row["_chem_comp_atom.{}".format(label.format("y"))]
        )
        z = conversions.str_to_float(
            row["_chem_comp_atom.{}".format(label.format("z"))]
        )

        atom_position = rdkit.Chem.rdGeometry.Point3D(x, y, z)
        conformer.SetAtomPosition(row.row_index, atom_position)

    conformer.SetProp("name", name)

    return conformer


def _parse_pdb_bonds(mol, cif_block, errors):
    """
    Setup bonds in the compound

    Args:
        mol (rdkit.Chem.rdchem.Mol): Molecule which receives bonds.
        cif_block (cif.Block): mmcif block Block from gemmi.
        errors (list[str]): Issues encountered while parsing.
    """
    if (
        "_chem_comp_atom." not in cif_block.get_mmcif_category_names()
        or "_chem_comp_bond." not in cif_block.get_mmcif_category_names()
    ):
        return
    atoms = cif_block.find("_chem_comp_atom.", ["atom_id"])
    bonds = cif_block.find(
        "_chem_comp_bond.", ["atom_id_1", "atom_id_2", "value_order"]
    )
    atoms_ids = list(atoms.find_column("_chem_comp_atom.atom_id"))
    for row in bonds:
        try:
            atom_1 = row["_chem_comp_bond.atom_id_1"]
            atom_1_id = atoms_ids.index(atom_1)
            atom_2 = row["_chem_comp_bond.atom_id_2"]
            atom_2_id = atoms_ids.index(atom_2)
            bond_order = helper.bond_pdb_order(row["_chem_comp_bond.value_order"])

            mol.AddBond(atom_1_id, atom_2_id, bond_order)
        except ValueError:
            errors.append(
                f"Error perceiving {atom_1} - {atom_2} bond in _chem_comp_bond"
            )
        except RuntimeError:
            errors.append(f"Duplicit bond {atom_1} - {atom_2}")


def _handle_implicit_hydrogens(mol):
    """Forbid atoms which does not have explicit Hydrogen partner to get
    implicit hydrogen.

    Args:
        mol (rkit.Chem.rdchem.Mol): Mol to be modified
    """
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 1:
            continue

        no_Hs = True
        for bond in atom.GetBonds():
            other = bond.GetOtherAtom(atom)
            if other.GetAtomicNum() == 1:
                no_Hs = False
                break

        atom.SetNoImplicit(no_Hs)


def _parse_pdb_descriptors(cif_block, cat_name, label="descriptor"):
    """Parse useful information from _pdbx_chem_comp_* category

    Args:
        cif_block (cif.Block): mmCIF Block object from gemmi
        cat_name (str): mmcif category with the
            descriptors info.
        label (str, optional): Defaults to 'descriptor'. Name of the
            category to be parsed.

    Returns:
        Descriptor: namedtuple with the property info
    """
    descriptors = []

    if cat_name not in cif_block.get_mmcif_category_names():
        return descriptors

    descriptors_block = cif_block.find(
        cat_name, [label, "type", "program", "program_version"]
    )
    for row in descriptors_block:
        d = Descriptor(
            type=cif.as_string(row[f"{cat_name}type"]),
            program=cif.as_string(row[f"{cat_name}program"]),
            program_version=cif.as_string(row[f"{cat_name}program_version"]),
            value=cif.as_string(row[f"{cat_name}{label}"]),
        )
        descriptors.append(d)

    return descriptors


def _parse_pdb_properties(cif_block):
    """Parse useful information from _chem_comp category

    Args:
        cif_block (cif.Block): mmcif block object from gemmi

    Returns:
        Properties: dataclass with the CCD properties.
    """
    properties = None
    if "_chem_comp." in cif_block.get_mmcif_category_names():
        mod_date = cif_block.find_value("_chem_comp.pdbx_modified_date")
        if cif.is_null(mod_date):
            d = date(1970, 1, 1)
        else:
            mod_date = mod_date.split("-")
            d = date(int(mod_date[0]), int(mod_date[1]), int(mod_date[2]))

        rel_status = ReleaseStatus.from_str(
            cif_block.find_value("_chem_comp.pdbx_release_status")
        )
        formula_weight = cif_block.find_value("_chem_comp.formula_weight")
        weight = 0.0 if cif.is_null(formula_weight) else cif.as_number(formula_weight)

        properties = CCDProperties(
            id=cif.as_string(cif_block.find_value("_chem_comp.id")),
            name=cif.as_string(cif_block.find_value("_chem_comp.name")),
            formula=cif.as_string(cif_block.find_value("_chem_comp.formula")),
            modified_date=d,
            pdbx_release_status=rel_status,
            weight=weight,
        )
    return properties


# endregion parse mmcif
