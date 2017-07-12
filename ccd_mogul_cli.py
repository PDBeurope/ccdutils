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
import os
import sys
from argparse import RawTextHelpFormatter
from pdb_ccd_mogul import PdbCCDMogul
from yattag import Doc

ANGSTROM = '&Aring;'
SIGMA = '&sigma;'


def __parse_command_line_args():
    """
    Sets up and parses the command line options.

    Returns:
         the arguments name space
    """
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=RawTextHelpFormatter)
    parser.add_argument('CIF', help='input PDB-CCD mmcif file (must be specified)')
    parser.add_argument('OUT_DIR', help='output directory for the report (must be specified)')
    parser.add_argument('--debug', action='store_true', help='turn on debug message logging output')
    return parser.parse_args()


def run():
    logger = logging.getLogger(' ')
    args = __parse_command_line_args()
    if args.debug:
        logging.basicConfig(level=logging.DEBUG, format='%(levelname)s: %(message)s',)
    cif_file = str(args.CIF)
    out_dir = str(args.OUT_DIR)
    logger.debug('input PDB-CCD cif file {}'.format(cif_file))
    logger.debug('output directory {}'.format(out_dir))

    try:
        pdb_ccd_mogul = PdbCCDMogul(file_name=cif_file)
    except ValueError as error_message:
        print('ERROR {}'.format(error_message))
        sys.exit(1)
    logging.debug('mogul results for {} bonds, {} angles, {} torsions and {} rings'.
                  format(len(pdb_ccd_mogul.store_bonds), len(pdb_ccd_mogul.store_angles),
                         len(pdb_ccd_mogul.store_torsions), len(pdb_ccd_mogul.store_rings)))
    if not os.path.isdir(out_dir):
        try:
            os.mkdir(out_dir)
            logging.debug('have made output directory {}'.format(out_dir))
        except OSError as error_message:
            print('ERROR cannot mkdir {} as {}'.format(out_dir, error_message))
            sys.exit(1)

    html_file = os.path.join(out_dir, 'index.html')
    html_text = prepare_html(pdb_ccd_mogul)
    with open(html_file, "w") as text_file:
        text_file.write(html_text)
    print('have written report to {}'.format(html_file))


def prepare_html(pdb_ccd_mogul):
    doc, tag, text, line = Doc().ttl()

    chem_comp_id = pdb_ccd_mogul.pdb_ccd_rdkit.chem_comp_id
    chem_comp_name = pdb_ccd_mogul.pdb_ccd_rdkit.chem_comp_name
    svg_diagram = pdb_ccd_mogul.pdb_ccd_rdkit.image_file_or_string(atom_labels=False)
    title = 'proof of concept - Mogul analysis of PDB-CCD coordinates for {}'.format(chem_comp_id)
    bond_title, bond_rows = prepare_bond_table(pdb_ccd_mogul)
    logging.debug(bond_title)

    with tag('html'):
        with tag('head'):
            with tag('title'):
                text(title)
            with tag('style'):
                text('table, th, td {border: 2px solid black; border-collapse: collapse;}')
                text('th, td { padding: 5px; text-align: center }')
        with tag('body'):
            with tag('h1'):
                text(title)
            with tag('ul', ):
                line('li', 'chem_comp_id =' + chem_comp_id)
                line('li', 'chem_comp_name = ' + chem_comp_name)
            doc.asis(svg_diagram)
            with tag('h2'):
                text('Mogul bond results')
            if len(pdb_ccd_mogul.store_bonds) == 0:
                line('p', 'no bonds found')
            else:
                with tag('table'):
                    with tag('tr'):
                        for item in bond_title:
                            with tag('th'):
                                doc.asis(item)
                    for row in bond_rows:
                        with tag('tr'):
                            for item in row:
                                with tag('td'):
                                    text(item)
    result = doc.getvalue()
    return result


def prepare_bond_table(pdb_ccd_mogul):
    title_row = ('atoms', 'actual in ' + ANGSTROM, 'Mogul mean in ' + ANGSTROM, 'difference in ' + ANGSTROM,
                 'Mogul ' + SIGMA + ' in ' + ANGSTROM, ' Mogul # hits', '|z-value|')
    rows = []
    for bond in sorted(pdb_ccd_mogul.store_bonds, key=lambda b: b.z_score, reverse=True):
        atoms = '-'.join(bond.atoms_ids)
        actual = '{:.3f}'.format(bond.value)
        mean = '{:.3f}'.format(bond.mean)
        difference = '{:.3f}'.format(bond.value - bond.mean)
        sigma = '{:.3f}'.format(bond.standard_deviation)
        nhits = '{}'.format(bond.nhits)
        z_score = '{:.1f}'.format(bond.z_score)
        rows.append((atoms, actual, mean, difference, sigma, nhits, z_score))

    return title_row, rows


if __name__ == "__main__":
    run()
