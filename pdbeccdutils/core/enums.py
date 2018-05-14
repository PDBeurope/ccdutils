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

from enum import IntEnum


class DepictionSource(IntEnum):
    Pubchem = 1
    Template = 2
    RdKit = 3


class ConformerType(IntEnum):
    Ideal = 1
    Model = 2
    Depiction = 3
    Computed = 4
    AllConformers = 5


class ReleaseStatus(IntEnum):
    """ an enumeration for pdbx_release_status
    allowed values include REL and HOLD, see:
    http://mmcif.wwpdb.org/dictionaries/mmcif_pdbx.dic/Items/_chem_comp.pdbx_release_status.html

    Notes
        An additional value 'NOT_SET' has been added for case where
        pdbx_release_status has not been set.
    """
    NOT_SET = 0
    DEL = 1
    HOLD = 2
    HPUB = 3
    OBS = 4
    REF_ONLY = 5
    REL = 6
