#!/usr/bin/env python
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
#
"""
Script for PDBeChem backend infrastructure.
Processes the wwPDB Chemical Components Dictionary file components.cif
producing files for

http://ftp.ebi.ac.uk/pub/databases/msd/pdbechem/

To do this components.cif is split into individual PDB chemical component
definitions cif files, sdf files, pdb files and image files.
In addition creates chem.xml and chem_comp.list for all components.
"""
import argparse
import logging
import os
from argparse import RawTextHelpFormatter

from split_components_cif import SplitComponentsCif
from utilities import create_directory_using_mkdir_unless_it_exists


clean_existing = True  # might want an update run mode later but for now remove existing directories/files
file_subdirs = 'mmcif', 'sdf', 'sdf_nh', 'sdf_r', 'sdf_r_nh', 'pdb', 'pdb_r', 'cml', 'xyz', 'xyz_r'

def create_parser():
    """
    Sets up parse the command line options.

    Returns:
         ArgumentParser parser
    """
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=RawTextHelpFormatter)
    parser.add_argument('COMPONENTS_CIF', help='input PDB-CCD components.cif file (must be specified)')
    parser.add_argument('OUTPUT_DIR', help='directory for output (must be specified)')
    parser.add_argument('--debug', action='store_true', help='turn on debug message logging output')
    return parser


def process_components_cif(components_cif, output_dir, debug):
    """
    processes components.cif for PDBeChem type output

    Args:
        components_cif (str): file name/path for components.cif or a test version
        output_dir (str): path for the directory where output will be written
        debug (bool): produce debug type logging

    Returns:

    """
    logger = logging.getLogger(' ')
    if debug:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.WARNING)
    logger.debug('components_cif={} output_dir={}'.format(components_cif, output_dir))
    create_directory_using_mkdir_unless_it_exists(output_dir)
    files_dir = os.path.join(output_dir, 'files')
    create_directory_using_mkdir_unless_it_exists(files_dir, clean_existing)
    file_subdirs_path = {}
    for subdir in file_subdirs:
        file_subdirs_path[subdir] = os.path.join(files_dir, subdir)
        create_directory_using_mkdir_unless_it_exists(file_subdirs_path[subdir])
        logger.debug('have created files subdir {}'.format(file_subdirs_path[subdir]))
    split_cc = SplitComponentsCif(components_cif)
    logger.debug('have opened {} and it contains {} individual CCD cif definitions '.
                 format(components_cif, len(split_cc.cif_dictionary)))
    for pdb_cc_rdkit in split_cc.individual_pdb_ccd_rdkit():
        chem_comp_id = pdb_cc_rdkit.chem_comp_id
        logger.debug('chem_comp_id={}'.format(chem_comp_id))
        for subdir in file_subdirs:
            if subdir == 'mmcif':
                file_type = '.cif'
            else:
                file_type = '.' + subdir[:3]
            output_file = os.path.join(file_subdirs_path[subdir], chem_comp_id + file_type)
            if subdir == 'mmcif':
                pdb_cc_rdkit.write_ccd_cif(output_file)
            elif subdir == 'sdf':
                pdb_cc_rdkit.sdf_file_or_string(file_name=output_file, ideal=True, hydrogen=True)
            elif subdir == 'sdf_nh':
                pdb_cc_rdkit.sdf_file_or_string(file_name=output_file, ideal=True, hydrogen=False)
            elif subdir == 'sdf_r':
                pdb_cc_rdkit.sdf_file_or_string(file_name=output_file, ideal=False, hydrogen=True)
            elif subdir == 'sdf_r_nh':
                pdb_cc_rdkit.sdf_file_or_string(file_name=output_file, ideal=False, hydrogen=False)
            elif subdir == 'pdb':
                pdb_cc_rdkit.pdb_file_or_string(file_name=output_file, ideal=True)
            elif subdir == 'pdb_r':
                pdb_cc_rdkit.pdb_file_or_string(file_name=output_file, ideal=False)
            elif subdir == 'cml':
                pdb_cc_rdkit.cml_file_or_string(file_name=output_file)
            elif subdir == 'xyz':
                pdb_cc_rdkit.xyz_file_or_string(file_name=output_file, ideal=True)
            elif subdir == 'xyz_r':
                pdb_cc_rdkit.xyz_file_or_string(file_name=output_file, ideal=False)
            else:
                raise NotImplementedError('unrecognized subdir {}'.format(subdir))
            if os.path.isfile(output_file):
                logger.debug('written file {}'.format(output_file))
            else:
                logger.warn('failed to write {}'.format(output_file))


def main():
    parser = create_parser()
    args = parser.parse_args()
    process_components_cif(args.COMPONENTS_CIF, args.OUTPUT_DIR, args.debug)

if __name__ == "__main__":
    main()
