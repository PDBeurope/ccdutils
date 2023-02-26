#!/usr/bin/env python
# software from PDBe: Protein Data Bank in Europe; https://pdbe.org
#
# Copyright 2021 EMBL - European Bioinformatics Institute
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
Set of methods to format data for gemmi parser
"""


def preprocess_cif_category(cif_block, label):
    """
    Checks if the category is present in gemmi.cif.Block object

    Args:
        cif_block (Block): mmcif Block from gemmi.
        label (str): name of the category

    Returns:
        str: Possible error encountered
    """
    if label not in cif_block.get_mmcif_category_names():
        return f"Namespace {label} does not exist."
