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
from pdbeccdutils.helpers import cif_tools, conversions, mol_tools
from pdbecif.mmcif_io import MMCIF2Dict

# categories that need to be 'fixed'
# str => list[str]
preprocessable_categories = [
    "_chem_comp_atom",
    "_chem_comp_bond",
    "_pdbx_chem_comp_identifier",
    "_pdbx_chem_comp_descriptor",
    "_chem_comp",
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

    cif_dict = list(MMCIF2Dict().parse(path_to_cif).values())[0]

    return _parse_pdb_mmcif(cif_dict, sanitize)


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

    for k, v in MMCIF2Dict().parse(path_to_cif).items():
        try:
            result_bag[k] = _parse_pdb_mmcif(v, sanitize)
        except CCDUtilsError as e:
            logging.error(f"ERROR: Data block {k} not processed. Reason: ({str(e)}).")

    return result_bag


# region parse mmcif
def _parse_pdb_mmcif(cif, sanitize=True):
    """
    Create internal representation of the molecule from mmcif format.

    Args:
        cif (dict): mmcif dictionary
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
        w = cif_tools.preprocess_cif_category(cif, c)

        if w:
            warnings.append(w)

    _parse_pdb_atoms(mol, cif)
    _parse_pdb_conformers(mol, cif)
    _parse_pdb_bonds(mol, cif, errors)
    _handle_implicit_hydrogens(mol)

    if sanitize:
        sanitized = mol_tools.sanitize(mol)

    descriptors = _parse_pdb_descriptors(
        cif, "_pdbx_chem_comp_descriptor", "descriptor"
    )
    descriptors += _parse_pdb_descriptors(
        cif, "_pdbx_chem_comp_identifier", "identifier"
    )
    properties = _parse_pdb_properties(cif["_chem_comp"])

    comp = Component(mol.GetMol(), cif, properties, descriptors)
    reader_result = CCDReaderResult(
        warnings=warnings, errors=errors, component=comp, sanitized=sanitized
    )

    return reader_result


def _parse_pdb_atoms(mol, cif):
    """
    Setup atoms in the component

    Args:
        mol (rdkit.Chem.rchem.Mol): Rdkit Mol object with the
            compound representation.
        cif (dict): mmCIF dictionary.
    """
    if "_chem_comp_atom" not in cif:
        return

    atoms = cif["_chem_comp_atom"]
    data = zip(
        atoms["atom_id"],
        atoms["type_symbol"],
        atoms["alt_atom_id"],
        atoms["pdbx_leaving_atom_flag"],
        atoms["charge"],
    )

    for atom_id, element, alt_atom_id, leaving_atom, charge in data:
        element = element if len(element) == 1 else element[0] + element[1].lower()
        isotope = None

        if element == "D":
            element = "H"
            isotope = 2
        elif element == "X":
            element = "*"

        atom = rdkit.Chem.Atom(element)
        # atom.SetChiralTag(_atom_chiral_tag(atoms['pdbx_stereo_config'][i]))
        atom.SetProp("name", atom_id)
        atom.SetProp("alt_name", alt_atom_id)
        atom.SetBoolProp("leaving_atom", leaving_atom == "Y")
        atom.SetFormalCharge(conversions.str_to_int(charge))

        if isotope is not None:
            atom.SetIsotope(isotope)

        mol.AddAtom(atom)


def _parse_pdb_conformers(mol, cif):
    """
    Setup model and ideal cooordinates in the rdkit Mol object.

    Args:
        mol (rdkit.Chem.rdchem.Mol): RDKit Mol object with the compound
            representation.
        cif (dict): mmCIF dictionary.
    """
    if "_chem_comp_atom" not in cif:
        return

    atoms = cif["_chem_comp_atom"]
    ideal = _setup_pdb_conformer(
        atoms, "pdbx_model_Cartn_{}_ideal", ConformerType.Ideal.name
    )
    model = _setup_pdb_conformer(atoms, "model_Cartn_{}", ConformerType.Model.name)

    mol.AddConformer(ideal, assignId=True)
    mol.AddConformer(model, assignId=True)


def _setup_pdb_conformer(atoms, label, name):
    """
    Setup a conformer

    Args:
        atoms (dict): mmcif category with the atom info.
        label (str): Namespace with the [x,y,z] coordinates.
        name (str): Conformer name.

    Returns:
        rdkit.Chem.rdchem.Conformer: Conformer of the component.
    """
    if not atoms:
        return

    conformer = rdkit.Chem.Conformer(len(atoms["atom_id"]))

    for i in range(len(atoms["atom_id"])):
        x = conversions.str_to_float(atoms[label.format(("x"))][i])
        y = conversions.str_to_float(atoms[label.format(("y"))][i])
        z = conversions.str_to_float(atoms[label.format(("z"))][i])

        atom_position = rdkit.Chem.rdGeometry.Point3D(x, y, z)
        conformer.SetAtomPosition(i, atom_position)

    conformer.SetProp("name", name)

    return conformer


def _parse_pdb_bonds(mol, cif, errors):
    """
    Setup bonds in the compound

    Args:
        mol (rdkit.Chem.rdchem.Mol): Molecule which receives bonds.
        cif (dict): mmcif dictionary.
        errors (list[str]): Issues encountered while parsing.
    """
    if "_chem_comp_atom" not in cif or "_chem_comp_bond" not in cif:
        return

    atoms = cif["_chem_comp_atom"]
    bonds = cif["_chem_comp_bond"]

    atoms_ids = atoms["atom_id"]
    data_struct = zip(bonds["atom_id_1"], bonds["atom_id_2"], bonds["value_order"])
    for atom_1, atom_2, order in data_struct:
        try:
            atom_1_id = atoms_ids.index(atom_1)
            atom_2_id = atoms_ids.index(atom_2)
            bond_order = _bond_pdb_order(order)

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


def _parse_pdb_descriptors(cif, cat_name, label="descriptor"):
    """
    Parse useful information from _pdbx_chem_comp_* category

    Args:
        pdbx_chem_comp_descriptors (dict): mmcif category with the
            descriptors info.
        label (str, optional): Defaults to 'descriptor'. Name of the
            category to be parsed.

    Returns:
        Descriptor: namedtuple with the property info
    """
    descriptors = []

    if cat_name not in cif:
        return descriptors

    for i in range(len(cif[cat_name][label])):
        d = Descriptor(
            type=cif[cat_name]["type"][i],
            program=cif[cat_name]["program"][i],
            value=cif[cat_name][label][i],
        )
        descriptors.append(d)

    return descriptors


def _parse_pdb_properties(chem_comp):
    """
    Parse useful information from _chem_comp category

    Args:
        chem_comp (dict): the mmcif category with _chem_comp info

    Returns:
        Properties: dataclass with the CCD properties.
    """
    properties = None
    if chem_comp:
        mod_date = chem_comp["pdbx_modified_date"][0].split("-")
        d = (
            date(1970, 1, 1)
            if mod_date[0] == "?"
            else date(int(mod_date[0]), int(mod_date[1]), int(mod_date[2]))
        )

        rel_status = chem_comp["pdbx_release_status"][0]
        rel_status = ReleaseStatus.from_str(chem_comp["pdbx_release_status"][0])
        weight = (
            0.0
            if chem_comp["formula_weight"][0] == "?"
            else float(chem_comp["formula_weight"][0])
        )

        properties = CCDProperties(
            id=chem_comp["id"][0],
            name=chem_comp["name"][0],
            formula=chem_comp["formula"][0],
            modified_date=d,
            pdbx_release_status=rel_status,
            weight=weight,
        )
    return properties


def _bond_pdb_order(value_order):
    """
    Transpils mmcif bond order into rdkit language

    Args:
        value_order (str): bond type as a str

    Returns:
        rdkit.Chem.rdchem.BondType: -- bond type
    """
    if value_order == "SING":
        return rdkit.Chem.rdchem.BondType(1)
    if value_order == "DOUB":
        return rdkit.Chem.rdchem.BondType(2)
    if value_order == "TRIP":
        return rdkit.Chem.rdchem.BondType(3)

    return None


def _atom_chiral_tag(tag):
    """Parse _chem_comp.pdbx_stereo.config from chem_comp

    Args:
        tag (str): R/S/N identification of chiral center.

    Returns:
        rdkit.Chem.ChiralType: Chiral center in RDKit language.
    """
    if tag == "N":
        return rdkit.Chem.ChiralType.CHI_UNSPECIFIED
    if tag == "S":
        return rdkit.Chem.ChiralType.CHI_TETRAHEDRAL_CCW
    if tag == "R":
        return rdkit.Chem.ChiralType.CHI_TETRAHEDRAL_CW

    return rdkit.Chem.ChiralType.CHI_UNSPECIFIED


# endregion parse mmcif
