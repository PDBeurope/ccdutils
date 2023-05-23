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
of proteins and creating Component representation of molecules.
"""

import os
import rdkit
from datetime import date
from collections import namedtuple
from gemmi import cif
from networkx import MultiDiGraph

from rdkit.Chem.rdMolDescriptors import CalcMolFormula
from pdbeccdutils.core import ccd_reader, models
from pdbeccdutils.core.boundmolecule import infer_bound_molecules
from pdbeccdutils.utils import config
from pdbeccdutils.core.component import Component
from pdbeccdutils.core.models import (
    CCDProperties,
    ConformerType,
    Descriptor,
)
from pdbeccdutils.helpers import cif_tools, conversions, mol_tools, helper


BMReaderResult = namedtuple(
    "BMReaderResult", ccd_reader.CCDReaderResult._fields + ("bound_molecule",)
)


def read_pdb_cif_file(
    path_to_cif: str,
    to_discard: set[str] = config.DISCARDED_RESIDUES,
    sanitize: bool = True,
    assembly: bool = False,
) -> list[BMReaderResult]:
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
    bms = infer_bound_molecules(path_to_cif, to_discard, assembly)
    for i, bm in enumerate(bms, start=1):
        bm_id = f"bm{i}"
        reader_result = infer_multiple_chem_comp(path_to_cif, bm, bm_id, sanitize)
        if reader_result:
            biomolecule_result.append(reader_result)

    return biomolecule_result


def infer_multiple_chem_comp(path_to_cif, bm, bm_id, sanitize=True):

    """Args:
        path_to_cif: Path to input structure
        bm: bound-molecules identified from input structure
        bm_id: ID of bound-molecule
        sanitize: True if bound-molecule need to be sanitized

    Returns:
        BMReaderResult: Namedtuple containing Component representation of bound-molecule

    """

    if bm.graph.number_of_nodes() <= 1 or bm.id != bm.orig_id:
        return

    cif_block = cif.read(path_to_cif).sole_block()
    (mol, warnings, errors) = _parse_pdb_mmcif(cif_block, bm.graph)
    sanitized = False
    if sanitize:
        sanitized = mol_tools.sanitize(mol)

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

    # define release category based on warnings or errors
    pdbx_release_status = models.ReleaseStatus.NOT_SET
    if warnings or errors:
        pdbx_release_status = models.ReleaseStatus.HOLD
    else:
        pdbx_release_status = models.ReleaseStatus.REL

    properties = CCDProperties(
        id=bm_id,
        name=bm.name,
        formula=CalcMolFormula(mol),
        modified_date=date.today().strftime("%Y-%m-%d"),
        pdbx_release_status=pdbx_release_status,
        weight="",
    )

    comp = Component(mol, cif_block, properties, descriptors)

    reader_result = BMReaderResult(
        warnings=warnings,
        errors=errors,
        component=comp,
        bound_molecule=bm,
        sanitized=sanitized,
    )

    return reader_result


def _parse_pdb_mmcif(
    cif_block: cif.Block, bm: models.BoundMolecule
) -> tuple[rdkit.Chem.rdchem.Mol, list[str], list[str]]:
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
    mol = _handle_hydrogens(mol)
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


def _parse_pdb_atoms(mol: rdkit.Chem.rdchem.Mol, atoms: dict[str, list[str]]):
    """Setup atoms of bound-molecules in the Component

    Args:
        mol: Rdkit Mol object of bound-molecule
        atoms: atoms of bound-molecules
    """

    for i in range(len(atoms["id"])):
        atom_id = atoms["label_atom_id"][i]
        chain = atoms["auth_asym_id"][i]
        res_name = atoms["label_comp_id"][i]
        res_id = int(atoms["auth_seq_id"][i])
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
        atom.SetProp("residue_id", residue_id)

        res_info = rdkit.Chem.AtomPDBResidueInfo()
        res_info.SetResidueName(res_name)
        res_info.SetResidueNumber(res_id)
        res_info.SetChainId(chain)
        atom.SetMonomerInfo(res_info)

        if isotope is not None:
            atom.SetIsotope(isotope)

        mol.AddAtom(atom)


def _parse_pdb_conformers(mol: rdkit.Chem.rdchem.Mol, atoms: dict[str, list[str]]):
    """Setup model cooordinates in the rdkit Mol object.

    Args:
        mol: RDKit Mol object of bound-molecule
        atoms: atoms of bound-molecule
    """
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


def _parse_pdb_bonds(
    mol: rdkit.Chem.rdchem.Mol,
    bm: MultiDiGraph,
    cif_block: cif.Block,
    errors: list[str],
):
    """Setup bonds in the rdkit Mol object

    Args:
        mol: RDKit Mol object of bound-molecule
        bm: bound-molecule
        errors: list of errors encountered while parsing.
    """
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


def _add_connections(
    mol: rdkit.Chem.rdchem.Mol, bm: MultiDiGraph, errors: list[str]
) -> None:
    """Add bonds between CCDs in the bound-molecule

    Args:
        mol: RDKit Mol object of bound-molecule
        bm: bound-molecule
        errors: list of errors encountered while parsing

    """
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


def get_chem_comp_bonds(cif_block: cif.Block, residue: str):
    """Returns _chem_comp_bond associated with a residue

    Args:
        cif_block: gemmi.cif.Block object of protein mmCIF file
        residue: CCD ID
    """

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
        mol (rdkit.Chem.rdchem.Mol): Rdkit Mol object with the
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


def _handle_hydrogens(mol):

    hydrogen_indices = [
        atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() == 1
    ]

    hydrogen_indices.sort(reverse=True)
    for index in hydrogen_indices:
        mol.RemoveAtom(index)

    mol.UpdatePropertyCache(strict=False)
    mol = rdkit.Chem.AddHs(mol, addCoords=True, addResidueInfo=True)
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 1:
            atom.SetProp("name", atom.GetSymbol() + str(atom.GetIdx()))
    return mol
