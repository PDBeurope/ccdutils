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


def bond_pdb_order(value_order):
    """
    Transpils mmcif bond order into rdkit language

    Args:
        value_order (str): bond type as a str

    Returns:
        rdkit.Chem.rdchem.BondType: -- bond type
    """
    if value_order.casefold() == "SING".casefold():
        return rdkit.Chem.rdchem.BondType(1)
    if value_order.casefold() == "DOUB".casefold():
        return rdkit.Chem.rdchem.BondType(2)
    if value_order.casefold() == "TRIP".casefold():
        return rdkit.Chem.rdchem.BondType(3)

    return None


def find_atom_index(mol, residue_id, atom_id):
    for atom in mol.GetAtoms():
        if atom.GetProp("residue_id") == residue_id and atom.GetProp("name") == atom_id:
            return atom.GetIdx()
