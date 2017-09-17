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
proof of concept - Mogul analysis of PDB-CCD coordinates

* Jiffy to read a pdb-ccd cif and run Mogul using the CSD Python API producing a html format report of the results.
* This should have coloured 2D diagrams showing outliers for bonds, angles, torsions
  and html tables listing outliers - possibly with javascript to list all.
* use buster-report as a model.
   http://grade.globalphasing.org/tut/erice_workshop/introtutorial/buster/00_MapOnly.report/ligand/
"""
import argparse
import logging
import sys
from argparse import RawTextHelpFormatter
from pdb_ccd_mogul import PdbCCDMogul

def __parse_command_line_args():
    """
    Sets up and parses the command line options.

    Returns:
         the arguments name space
    """
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=RawTextHelpFormatter)
    parser.add_argument('CIF', help='input PDB-CCD mmcif file (must be specified)')
    parser.add_argument('HTML', help='output html report filename (must be specified)')
    parser.add_argument('--debug', action='store_true', help='turn on debug message logging output')
    return parser.parse_args()


def run():
    logger = logging.getLogger(' ')
    args = __parse_command_line_args()
    if args.debug:
        logging.basicConfig(level=logging.DEBUG, format='%(levelname)s: %(message)s',)
    cif_file = str(args.CIF)
    html_file = str(args.HTML)
    logger.debug('input PDB-CCD cif file {}'.format(cif_file))
    logger.debug('output html file {}'.format(html_file))

    try:
        pdb_ccd_mogul = PdbCCDMogul(file_name=cif_file)
        logger.debug('ideal coords {}'.format(pdb_ccd_mogul.pdb_ccd_rdkit.ideal_xyz))
        pdb_ccd_mogul.run_mogul()
    except ValueError as error_message:
        print('ERROR {}'.format(error_message))
        sys.exit(1)
    logging.debug('mogul results for {} bonds, {} angles, {} torsions and {} rings'.
                  format(len(pdb_ccd_mogul.store_bonds), len(pdb_ccd_mogul.store_angles),
                         len(pdb_ccd_mogul.store_torsions), len(pdb_ccd_mogul.store_rings)))
    pdb_ccd_mogul.prepare_file_html(html_file)
    print('have written report to {}'.format(html_file))


if __name__ == "__main__":
    run()
