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

import csv
import os

from rdkit import Chem
from rdkit import rdBase
from rdkit.Chem import AllChem

from pdbeccdutils.utils import config
from pdbeccdutils.computations import DepictionManager
from pdbeccdutils.core import DepictionSource


class FragmentLibrary:
    """Implementation of fragment library.

    Returns:
        pdbeccdutils.core.FragmentLibrary: Instance of the object
    """

    def __init__(self, path=config.fragment_library, header=True, delimiter='\t', quotechar='"'):
        self.name = os.path.basename(path).split('.')[0]
        self.library = self._read_in_library(path, header, delimiter, quotechar)

    def _read_in_library(self, path, header, delimiter, quotechar):
        """Fragment library parser

        Args:
            path (str): Path to the fragment library.
            header (bool): Whether or not the library contains a header.
            delimiter (str): Delimiter symbol.
            quotechar (str): Quotechar symbol.
        """
        rdBase.DisableLog('rdApp.*')
        library = {}
        depiction_manager = DepictionManager()

        with open(path, 'r') as csvfile:
            library_reader = csv.reader(csvfile, delimiter=delimiter, quotechar=quotechar)
            if header:
                next(library_reader)

            for row in library_reader:
                mol = None
                if row[1] == 'SMARTS':
                    mol = Chem.MolFromSmarts(row[2])
                else:
                    mol = Chem.MolFromSmiles(row[2])
                    if row[1] == 'LIKE':
                        [bond.SetBondType(Chem.rdchem.BondType.UNSPECIFIED) for bond in mol.GetBonds()]
                Chem.SanitizeMol(mol, catchErrors=True)
                _generate_fragment_conformer(mol)
                mol = _depict_fragment(row[0], mol, depiction_manager)

                library[row[0]] = mol
        rdBase.EnableLog('rdApp.*')
        return library

    def to_image(self, path):
        """Export image with all fragments.

        Args:
            path (str): Destination of the image
        """
        img = Chem.Draw.MolsToGridImage(list(self.library.values()),
                                        legends=list(self.library.keys()), molsPerRow=10)
        img.save(path)


def _depict_fragment(name, mol, depiction_manager):
    """Get nice 2D coords of the fragments

    Args:
        name (str): Molecule name.
        mol (rdkit.Chem.rdchem.Mol): Molecule to be processed.
        depiction_manager (pdbeccdutils.computations.DepictionManager):
            DepictionManager instance to create nice 2D layout

    Returns:
        rdkit.Chem.rdChem.Mol: Mol with 2D coords.
    """
    result = depiction_manager.depict_molecule(name, mol)
    return result.mol if result.source is not DepictionSource.Failed else mol


def _generate_fragment_conformer(mol):
    """Generate 3D coordinates for the fragment so it can be depicted
    by the DepictionManager

    Args:
        mol (rdkit.Chem.rdchem.Mol): Molecule to be processed.
    """
    options = AllChem.ETKDGv2()
    options.clearConfs = False

    try:
        AllChem.EmbedMolecule(mol, options)
        AllChem.UFFOptimizeMolecule(mol)
    except Exception:
        pass  # don't care if it fails
