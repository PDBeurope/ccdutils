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

"""Generic helper functions that may be re-used
"""

import rdkit


def find_element_in_list(array, element):
    """
    Finds an element in an array. Does not crash if not found

    Args:
        array (list): list to be searched
        element (any): element to be found

    Returns:
        int: Index of the element or None if the element is not found.
    """
    try:
        index = array.index(element)
        return index
    except ValueError:
        return None


def bond_pdb_order(value_order):
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


def find_atom_index(mol, residue_id, atom_id):
    for atom in mol.GetAtoms():
        if atom.GetProp("residue_id") == residue_id and atom.GetProp("name") == atom_id:
            return atom.GetIdx()
