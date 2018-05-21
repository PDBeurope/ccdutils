# software from PDBe: Protein Data Bank in Europe; http://pdbe.org
#
# Copyright 2018 EMBL - European Bioinformatics Institute
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
import os
import sys
import xml.etree.ElementTree as ET
from xml.dom import minidom

import rdkit

import pdbeccdutils
from pdbeccdutils.core import structure_writer as writer
from pdbeccdutils.core import ConformerType, DepictionSource, ccd_reader
from pdbeccdutils.utils import DepictionManager, PubChemDownloader


def create_parser():
    """
    Sets up parse the command line options.

    Returns:
         ArgumentParser parser
    """
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    add_arg = parser.add_argument
    add_arg('components_cif', help='input PDB-CCD components.cif file (must be specified)')
    add_arg('--generic_templates', default='',
            help='Path to the directory with generic templates in sdf format.')
    add_arg('--pubchem_templates', default='',
            help='Path to the directory with pubchem templates in sdf format.')
    add_arg('--output_dir', '-o',
            help='create an output directory with files suitable for PDBeChem ftp directory')
    add_arg('--test_first', type=int,
            help='only process the first TEST_FIRST chemical component definitions (for testing).')
    add_arg('--library',
            help='use this fragment library in place of the one supplied with the code.')
    add_arg('--debug', action='store_true', help='turn on debug message logging output')

    return parser


def check_args(args):
    if not os.path.isfile(args.components_cif):
        print(f'{args.components_cif} does not exist', file=sys.stderr)
        sys.exit(os.EX_NOINPUT)

    if not os.path.isdir(args.output_dir):
        os.makedirs(args.output_dir)

    if args.test_first is not None:
        if args.test_first < 1:
            print(f'Test_first mode needs to have at least 1 component.', file=sys.stderr)
            sys.exit(os.EX_NOINPUT)


def pdbechem_pipeline(args):
    """
    Processes components.cif for PDBeChem type output

    Args:
        args: an argparse namespace containing the required arguments
    """

    logger = logging.getLogger(' ')
    level = logging.DEBUG if args.debug else logging.WARNING
    format = '[%(asctime)-15s]  %(message)s'
    logging.basicConfig(level=level, format=format, datefmt='%a, %d %b %Y %H:%M:%S')
    logger.debug(f'PDBeChem pipeline using:')
    logger.debug(f'pdbeccdutils core v. {pdbeccdutils.__version__}, RDKit v. {rdkit.__version__}')

    settings = 'Settings:\n'
    settings += '{:29s}{:25s}{}\n'.format('', 'components_cif', args.components_cif)
    settings += '{:29s}{:25s}{}\n'.format('', 'output_dir', args.output_dir)
    settings += '{:29s}{:25s}{}\n'.format('', 'generic_templates', args.generic_templates)
    settings += '{:29s}{:25s}{}\n'.format('', 'pubchem_templates', args.pubchem_templates)
    settings += '{:29s}{:25s}{}\n'.format('', 'library', args.library)
    settings += '{:29s}{:25s}{}'.format('', 'DEBUG', ('ON' if args.debug else 'OFF'))

    logger.debug(settings)

    depictions = DepictionManager(args.generic_templates, args.pubchem_templates)
    pubchem_templates = PubChemDownloader(args.pubchem_templates)
    ccd_reader_result = ccd_reader.read_pdb_components_file(args.components_cif)
    counter = len(ccd_reader_result) if args.test_first is None else args.test_first

    ids = sorted(list(ccd_reader_result.keys())[:counter])

    chem_comp_xml = ET.Element('chemCompList')
    for k, v in ccd_reader_result.items():
        try:
            if len(v.warnings) > 0:
                logger.debug(f'{k} | warnings: {";".join(v.warnings)}')
            if len(v.errors) > 0:
                logger.debug(f'{k} | errors: {";".join(v.errors)}')

            component_downloaded = pubchem_templates.download_template(v.component)
            if component_downloaded:
                logger.debug(f'{k} | downloaded new pubchem template')

            process_component(k, v.component, logger, depictions, args.output_dir)

            xml_repr = writer.to_xml_xml(v.component)
            chem_comp_xml.append(xml_repr)
        except Exception:
            logger.debug(f'{k} | FAILURE.')
        counter -= 1

        if counter == 0:
            break

    # write chem_comp_xml
    xml_str = ET.tostring(chem_comp_xml, encoding='utf-8', method='xml')
    pretty = minidom.parseString(xml_str)

    with open(os.path.join(args.output_dir, "chem_comp_list.xml"), 'w') as f:
        f.write(pretty.toprettyxml(indent="  "))

    # write ligand list
    with open(os.path.join(args.output_dir, 'chem_comp.list'), 'w') as f:
        f.write("\n".join(ids))


def process_component(ccd_id, component, logger, depictions, output_dir):
    """Process the component and write all necessary files to the output
    directory.

    Args:
        ccd_id (str): Component id
        component (pdbeccdutils.core.Component): Component to be processed.
        logger (logging.loger): Log of the application.
        depictions (pdbeccdutils.utils.DepictionManager): Depiction manager.
        output_dir (str): Output directory.
    """
    parent_dir = os.path.join(output_dir, ccd_id[0], ccd_id)

    os.makedirs(parent_dir)

    if not component.sanitize():
        logger.debug(f'{ccd_id} | sanitization issue.')

    if component.inchikey_from_rdkit_matches_ccd():
        logger.debug(f'{ccd_id} | inchikey mismatch.')

    if component.has_degenerated_conformer(ConformerType.Ideal):
        logger.debug(f'{ccd_id} | has degenerated ideal coordinates.')
        result = component.compute_3d()
        if not result:
            logger.debug(f'{ccd_id} | 3D conformation could not be generated.')

    # write images
    depiction_result = component.compute_2d(depictions)

    if depiction_result.source == DepictionSource.Failed:
        logger.debug(f'{ccd_id} | failed to generate 2D image.')
    else:
        if depiction_result.score > 0.99:
            logger.debug(f'{ccd_id} | collision free image could not be generated.')

    # write images
    component.export_2d_svg(os.path.join(parent_dir, f'{ccd_id}_small.svg'), width=100)
    component.export_2d_svg(os.path.join(parent_dir, f'{ccd_id}_medium.svg'), width=300)
    component.export_2d_svg(os.path.join(parent_dir, f'{ccd_id}_large.svg'), width=450)
    component.export_2d_svg(os.path.join(parent_dir, f'{ccd_id}_large_names.svg'), width=450, names=True)

    # export structures in PDB,SDF,CIF and CML format
    writer.write_molecule(os.path.join(parent_dir, f'{ccd_id}_model.sdf'), component, False, ConformerType.Model)
    writer.write_molecule(os.path.join(parent_dir, f'{ccd_id}_ideal.sdf'), component, False, ConformerType.Ideal)
    writer.write_molecule(os.path.join(parent_dir, f'{ccd_id}_no_h_model.sdf'), component, True, ConformerType.Model)
    writer.write_molecule(os.path.join(parent_dir, f'{ccd_id}_no_h_ideal.sdf'), component, True, ConformerType.Ideal)
    writer.write_molecule(os.path.join(parent_dir, f'{ccd_id}_model.pdb'), component, False, ConformerType.Model)
    writer.write_molecule(os.path.join(parent_dir, f'{ccd_id}_ideal.pdb'), component, False, ConformerType.Ideal)
    writer.write_molecule(os.path.join(parent_dir, f'{ccd_id}_ideal.pdb'), component, False, ConformerType.Ideal)
    writer.write_molecule(os.path.join(parent_dir, f'{ccd_id}.cml'), component)
    writer.write_molecule(os.path.join(parent_dir, f'{ccd_id}.cif'), component, False)


def main():
    """Runs the PDBeChem pipeline
    """
    parser = create_parser()
    args = parser.parse_args()

    check_args(args)
    pdbechem_pipeline(args)
