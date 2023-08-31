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
import glob
import os


def tests_dir():
    """
    directory of the tests directory where this python file is.

     Returns:
         str the name of the directory is

    """
    return os.path.dirname(os.path.abspath(__file__))


test_cif = os.path.join(tests_dir(), "ccd_mmcif_test_files", "random_sample")
prd_test_cif = os.path.join(tests_dir(), "prd_mmcif_test_files")
test_depiction = os.path.join(tests_dir(), "ccd_mmcif_test_files", "depiction_test")
test_boundmolecules = os.path.join(
    tests_dir(),
    "updated_mmcif_test_files",
)


def cif_filename(code):
    path = os.path.join(test_cif, f"{code}.cif")

    return path if os.path.isfile(path) else os.path.join(test_depiction, f"{code}.cif")


def prd_cif_filename(code):
    path = os.path.join(prd_test_cif, f"{code}.cif")
    if os.path.isfile(path):
        return path


def updated_mmcif_filename(pdb_id):
    path = os.path.join(test_boundmolecules, pdb_id, f"{pdb_id}_processed.cif.gz")
    if os.path.isfile(path):
        return path


def fixed_mmcif_filename(pdb_id):
    path = os.path.join(test_boundmolecules, pdb_id, f"{pdb_id}_processed.cif.gz")
    if os.path.isfile(path):
        return path


def bio_assembly_filename(pdb_id, au_fallback):
    if not au_fallback:
        path = os.path.join(test_boundmolecules, pdb_id, f"{pdb_id}_bio.cif.gz")
    else:
        path = os.path.join(test_boundmolecules, pdb_id, f"{pdb_id}_processed.cif.gz")
    if os.path.isfile(path):
        return path


def bms_filename(pdb_id):
    path = os.path.join(test_boundmolecules, pdb_id, "bound_molecules.json")
    if os.path.isfile(path):
        return path


def remove_discarded_ligands(category, field, discarded_ligands):
    new_category = {key: [] for key in category}
    for i in range(len(category[field])):
        if category[field][i] not in discarded_ligands:
            for key in new_category:
                new_category[key].append(category[key][i])

    return new_category


def unl_model():
    return os.path.join(tests_dir(), "ccd_mmcif_test_files", "UNL.cif")


def supply_list_of_sample_cifs():
    """
    returns the list of sample pdb ccd cifs for test.

    Args:

    Returns:
        list of filenames
    """
    return sorted(glob.glob(os.path.join(test_cif, "*.cif")))
