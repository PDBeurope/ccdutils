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
from argparse import RawTextHelpFormatter


def __parse_command_line_args():
    """
    Sets up and parses the command line options.

    Returns:
         the arguments name space
    """
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=RawTextHelpFormatter)
    parser.add_argument('COMPONENTS_CIF', help='input PDB-CCD components.cif file (must be specified)')
    parser.add_argument('OUTPUT_DIR', help='directory for output (must be specified)')
    parser.add_argument('--debug', action='store_true', help='turn on debug message logging output')
    return parser.parse_args()


def run():
    logger = logging.getLogger(' ')
    args = __parse_command_line_args()
    if args.debug:
        logging.basicConfig(level=logging.DEBUG)
    components_cif = str(args.COMPONENTS_CIF)
    output_dir = str(args.OUTPUT_DIR)
    logger.debug('components_cif={} output_dir ={}'.format(components_cif, output_dir))

if __name__ == "__main__":
    run()
