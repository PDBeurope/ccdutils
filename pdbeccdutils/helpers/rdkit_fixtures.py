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

import numpy as np
import rdkit


def fix_conformer(conformer):
    """In place fixing of rdkit conformer.
    In certain cases the resulting conformer (mainly from depiction process)
    can contain not valid atom coordinatesp [NaN, NaN, NaN]. This
    results in errors in downstream processes so that it is easier
    to fix it when the problem occurs

    Args:
        conformer (rdkit.Chem.rdchem.Conformer): RDKit conformer
    """    
    positions = conformer.GetPositions()

    for index, pos in enumerate(positions):
        if all(np.isnan(pos)):
            new_pos = rdkit.Chem.rdGeometry.Point3D(0, 0, 0)
            conformer.SetAtomPosition(index, new_pos)
