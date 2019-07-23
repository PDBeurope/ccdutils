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


test_cif = os.path.join(tests_dir(), 'ccd_mmcif_test_files', 'random_sample')
test_depiction = os.path.join(tests_dir(), 'ccd_mmcif_test_files', 'depiction_test')
test_cut_down_components_cif = os.path.join(tests_dir(), 'components_cif', 'cut_down_components.cif')


def cif_filename(code):
    path = os.path.join(test_cif, f'{code}.cif')

    return path if os.path.isfile(path) else os.path.join(test_depiction, f'{code}.cif')

def unl_model():
    return os.path.join(tests_dir(), 'ccd_mmcif_test_files', 'UNL.cif')

def supply_list_of_sample_cifs():
    """
    returns the list of sample pdb ccd cifs for test.

    Args:

    Returns:
        list of filenames
    """
    return sorted(glob.glob(os.path.join(test_cif, '*.cif')))
