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

import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdCoordGen

import pdbeccdutils.utils.config as config


class FragmentLibrary:
    """Implementation of fragment library.

    Returns:
        pdbeccdutils.core.FragmentLibrary: Instance of the object
    """

    def __init__(self, path=config.fragment_library, header=True, delimiter='\t', quotechar='"'):
        self.name = os.path.basename(path).split('.')[0]
        self._read_in_library(path, header, delimiter, quotechar)

    def _read_in_library(self, path, header, delimiter, quotechar):
        """Fragment library parser

        Args:
            path (str): Path to the fragment library.
            header (bool): Whether or not the library contains a header.
            delimiter (str): Delimiter symbol.
            quotechar (str): Quotechar symbol.
        """
        rdkit.rdBase.DisableLog('rdApp.*')
        rdCoordGen.SetDefaultTemplateFileDir(config.coordgen_templates)
        rdCoordGen.CoordGenParams.coordgenScaling = 33  # to get single bond length 1.5

        self.library = {}

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
                rdCoordGen.AddCoords(mol)
                self.library[row[0]] = mol

        rdkit.rdBase.EnableLog('rdApp.*')

    def to_image(self, path):
        """Export image with all fragments.

        Args:
            path (str): Destination of the image
        """
        img = Chem.Draw.MolsToGridImage(list(self.library.values()),
                                        legends=list(self.library.keys()), molsPerRow=10)
        img.save(path)

    def generate_conformers(self):
        """Generate 3D coordinates for the fragment library.

        """
        rdkit.rdBase.DisableLog('rdApp.*')
        options = AllChem.ETKDGv2()
        options.clearConfs = False

        for k, v in self.library.items():
            try:
                AllChem.EmbedMolecule(v, options)
                AllChem.UFFOptimizeMolecule(v)
            except Exception:
                pass  # don't care if it fails
        rdkit.rdBase.EnableLog('rdApp.*')
