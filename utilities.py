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
import errno
import glob
import shutil
import os


def create_directory_using_mkdir_unless_it_exists(path, clean_existing=False):
    if clean_existing:
        if os.path.isdir(path):
            shutil.rmtree(path)
    try:
        os.mkdir(path)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise


def this_script_dir():
    """
    gives the name of directory were this python file

    Returns:
        str the name of the directory is

    """
    return os.path.dirname(os.path.abspath(__file__))


test_cif_path_name = os.path.join(this_script_dir(), 'tests', 'ccd_mmcif_test_files')
test_comparison_files_path = os.path.join(this_script_dir(), 'tests', 'comparison_files')
test_components_cif_first_file_comps = os.path.join(this_script_dir(),
                                                    'tests', 'components_cif', 'components.cif.first_five_comps')


def cif_filename(code):
    return os.path.join(test_cif_path_name, code + '.cif')


def supply_list_of_sample_cifs():
    """
    returns the list of sample pdb ccd cifs for test.

    Args:

    Returns:
        list of filenames
    """
    return sorted(glob.glob(os.path.join(test_cif_path_name, '*.cif')))


def file_name_in_tsts_out(file_name, remove_existing=True):
    """
    creates the subdirectory "out" in the "tests" subdirectory (if necessary)
    and returns the file_name in this directory, If the file already exists it will remove it.
    (Cannot call method file_name_in_tests_out otherwise nose thinks it is a test)

    Args:
        file_name (str):  the name for the file
        remove_existing (bool): remove existing file?

    Returns:
        str: the filename in the subdirectory tests/out
    """
    subdir = os.path.join(this_script_dir(), 'tests', 'out')
    if not os.path.isdir(subdir):
        os.mkdir(subdir)
    out_file_name = os.path.join(subdir, file_name)
    if remove_existing and os.path.isfile(out_file_name):
        os.remove(out_file_name)
    return out_file_name
