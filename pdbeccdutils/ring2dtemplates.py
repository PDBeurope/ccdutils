# software from PDBe: Protein Data Bank in Europe; http://pdbe.org
#
# Copyright 2017 EMBL - European Bioinformatics Institute
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
#
"""
module to supply an ordered dictionary of rdkit mols with good 2D coordinates for complicated ring systems.
"""
import logging
from collections import OrderedDict
# noinspection PyPackageRequirements
from rdkit import Chem
from pdbeccdutils.utilities import ring2dtemplates_file_path

ring_2d_templates = OrderedDict()


def supply_ring_2d_templates():
    """
    supplies rdkit mols with good 2D coordinates for complicated ring systems to be used as templates for drawing
    2D coordinates

    Returns:
        an ordered dictionary name -> rdkit mol with the coordinates.

    Notes:
        The mols have the bond type set to unspecified so they will match regardless of bond order so a chlorin
        ring should match the porphin template.
    """
    if len(ring_2d_templates) == 0:
        _load_ring_2d_templates()
    return ring_2d_templates


def _load_ring_2d_templates():
    logging.debug('loading ring templates from file {}'.format(ring2dtemplates_file_path))
    supplier = Chem.SDMolSupplier(ring2dtemplates_file_path)
    for rdkit_mol in supplier:
        name = rdkit_mol.GetProp("_Name")
        logging.debug('ring template {} has {} atoms'.format(name, rdkit_mol.GetNumAtoms()))
        for bond in rdkit_mol.GetBonds():
            bond.SetBondType(Chem.rdchem.BondType.UNSPECIFIED)
        ring_2d_templates[name] = rdkit_mol
