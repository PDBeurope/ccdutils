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
from typing import Dict

import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem, rdCoordGen

import pdbeccdutils.utils.config as config
from pdbeccdutils.core.models import FragmentEntry


class FragmentLibrary:
    """Implementation of fragment library.
    """

    def __init__(self, path: str = config.fragment_library, header: bool = True,
                 delimiter: str = '\t', quotechar: str = '"') -> None:
        self.library: Dict[str, FragmentEntry] = {}
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
                        for bond in mol.GetBonds():
                            bond.SetBondType(Chem.rdchem.BondType.UNSPECIFIED)

                Chem.SanitizeMol(mol, catchErrors=True)
                rdCoordGen.AddCoords(mol)

                self.library[row[0]] = FragmentEntry(row[0], row[6], mol)

        rdkit.rdBase.EnableLog('rdApp.*')

    def to_image(self, path, source=''):
        """Export image with all fragments.

        Args:
            path (str): Destination of the image
            source (str): Select a source which fragments are going to
                be drawn.

        """
        use_svg = path[-3:].lower() == 'svg'
        if source:
            temp = {k: v for k, v in self.library.items() if v.source == source}
        else:
            temp = self.library

        mols = [v.mol for k, v in temp.items()]
        names = [v.name for k, v in temp.items()]
        img = Chem.Draw.MolsToGridImage(mols, legends=names, molsPerRow=10, useSVG=use_svg)

        if use_svg:
            with open(path, 'w') as f:
                f.write(img)
        else:
            img.save(path)

    def generate_conformers(self):
        """Generate 3D coordinates for the fragment library.
        """
        rdkit.rdBase.DisableLog('rdApp.*')
        options = AllChem.ETKDGv2()
        options.clearConfs = False

        for v in self.library.values():
            try:
                AllChem.EmbedMolecule(v.mol, options)
                AllChem.UFFOptimizeMolecule(v.mol)
            except Exception:
                pass  # don't care if it fails
        rdkit.rdBase.EnableLog('rdApp.*')
