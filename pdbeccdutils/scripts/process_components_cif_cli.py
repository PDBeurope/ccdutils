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
In addition creates chem_comp.xml and chem_comp.list for all components.
"""
import argparse
import logging
from argparse import RawTextHelpFormatter


def create_parser():
    """
    Sets up parse the command line options.

    Returns:
         ArgumentParser parser
    """
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=RawTextHelpFormatter)
    add_arg = parser.add_argument
    add_arg('COMPONENTS_CIF', help='input PDB-CCD components.cif file (must be specified)')
    add_arg('--output_dir', '-o',
            help='create an output directory with files suitable for PDBeChem ftp directory')
    add_arg('--chem_comp_xml',
            help='write chem_comp.xml file to this file.')
    add_arg('--test_first', type=int,
            help='only process the first TEST_FIRST chemical component definitions (for testing).')
    add_arg('--library',
            help='use this fragment library in place of the one supplied with the code.')
    add_arg('--debug', action='store_true', help='turn on debug message logging output')
    return parser


def process_components_cif(args):
    """
    processes components.cif for PDBeChem type output

    Args:
        args: an argparse namespace containing the required arguments
    """
    components_cif = args.COMPONENTS_CIF
    output_dir = args.output_dir
    test_first = args.test_first
    debug = args.debug
    library = args.library
    logger = logging.getLogger(' ')
    if debug:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.WARNING)
    logger.debug('components_cif={} output_dir={}'.format(components_cif, output_dir))
    raise NotImplementedError('process_components_cif to be written')

def main():
    parser = create_parser()
    args = parser.parse_args()
    process_components_cif(args)
