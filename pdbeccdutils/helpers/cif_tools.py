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
Set of methods to format data for pdbecif parser
"""


def preprocess_cif_category(cif, label):
    """
    The mmcif dictionary values are either str or list(), which is a bit
    tricky to work with. This method makes list() of all of them in
    order to parse all of the in the same way.

    Args:
        vif (dict): mmcif category with the parser output.
        label (str): name of the category

    Returns:
        str: Possible error encountered
    """
    if label not in cif:
        return f"Namespace {label} does not exist."

    check_element = list(cif[label].keys())[0]
    values = (
        cif[label]
        if isinstance(cif[label][check_element], list)
        else {k: [v] for k, v in cif[label].items()}
    )
    cif[label] = values

    return None


def post_process_cif_category(cif, category_name):
    """Single value category needs to be string rather than array
    with a single value. Also if the category contains nothing, it should
    be removed.

    Args:
        cif_copy (dict of str: dict): Dictionary like structure of
            the CIF file.
        category_name (str): Category name to be accounted for
    """
    if not cif[category_name]:  # nothing in the category => should be removed
        cif.pop(category_name)
        return

    for k, v in cif[category_name].items():
        if isinstance(v, list):
            if len(v) == 1:
                cif[category_name][k] = v[0]

            if not v:
                cif.pop(category_name)
                return
