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
from pdbeccdutils.core import ConformerType, DepictionSource
from pdbeccdutils.core import ccd_reader, FragmentLibrary
from pdbeccdutils.utils import DepictionManager, PubChemDownloader, config


def create_parser():
    """
    Sets up parse the command line options.

    Returns:
         ArgumentParser parser
    """
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    add_arg = parser.add_argument
    add_arg('components_cif', help='Input PDB-CCD components.cif file (must be specified)')
    add_arg('--general_templates', default=config.general_templates, type=str,
            help='Path to the directory with general templates in sdf format.')
    add_arg('--pubchem_templates', default='',
            help='Path to the directory with pubchem templates in sdf format.')
    add_arg('--output_dir', '-o', required=True,
            help='Create an output directory with files suitable for PDBeChem ftp directory')
    add_arg('--test_first', type=int,
            help='Only process the first TEST_FIRST chemical component definitions (for testing).')
    add_arg('--library', default=config.fragment_library,
            help='Use this fragment library in place of the one supplied with the code.')
    add_arg('--debug', action='store_true', help='Turn on debug message logging output')

    return parser


def check_args(args):
    """Validate suplied arguments.

    Args:
        args (ArgumentParser): an argparse namespace containing the required arguments
    """
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
        args (ArgumentParser): an argparse namespace containing the
            required arguments
    """

    # set up logging
    logger = _set_up_logger(args)
    settings = _log_settings(args)
    logger.debug(settings)

    # set up pipeline helpers
    depictions = DepictionManager(args.pubchem_templates, args.general_templates)
    pubchem_templates = PubChemDownloader(args.pubchem_templates) if os.path.isdir(args.pubchem_templates) else None
    fragment_library = FragmentLibrary(args.library)

    logger.debug(f'Reading in {args.components_cif} file...')
    ccd_reader_result = ccd_reader.read_pdb_components_file(args.components_cif)
    counter = len(ccd_reader_result) if args.test_first is None else args.test_first

    chem_comp_xml = ET.Element('chemCompList')
    ids = sorted(list(ccd_reader_result.keys())[:counter])

    for key, ccd_reader_result in ccd_reader_result.items():
        try:
            logger.debug(f'{key} | processing...')
            process_single_component(args, ccd_reader_result, fragment_library, chem_comp_xml,
                                     depictions, pubchem_templates, logger)
        except Exception as e:
            logger.debug(f'{key} | FAILURE.')
        counter -= 1

        if counter == 0:
            break

    # write components.cif wide data
    write_components_xml(args, chem_comp_xml)
    with open(os.path.join(args.output_dir, 'chem_comp.list'), 'w') as f:
        f.write("\n".join(ids))


def _set_up_logger(args):
    """Set up application level logging.

    Args:
        args (ArgumentParser): Parsed arguments.

    Returns:
        logging.logger: Application log.
    """

    logger = logging.getLogger(' ')
    level = logging.DEBUG if args.debug else logging.WARNING
    format = '[%(asctime)-15s]  %(message)s'
    logging.basicConfig(level=level, format=format, datefmt='%a, %d %b %Y %H:%M:%S')
    logger.debug(f'PDBeChem pipeline using:')
    logger.debug(f'pdbeccdutils core v. {pdbeccdutils.__version__}, RDKit v. {rdkit.__version__}')

    return logger


def _log_settings(args):
    """Compose initial message about application settions.

    Args:
        args (ArgumentParser): Application arguments.

    Returns:
        str: initial application message.
    """

    settings = 'Settings:\n'
    settings += '{:29s}{:25s}{}\n'.format('', 'components_cif', args.components_cif)
    settings += '{:29s}{:25s}{}\n'.format('', 'output_dir', args.output_dir)
    settings += '{:29s}{:25s}{}\n'.format('', 'general_templates', args.general_templates)
    settings += '{:29s}{:25s}{}\n'.format('', 'pubchem_templates', args.pubchem_templates)
    settings += '{:29s}{:25s}{}\n'.format('', 'library', args.library)
    settings += '{:29s}{:25s}{}'.format('', 'DEBUG', ('ON' if args.debug else 'OFF'))

    return settings


def process_single_component(args, ccd_reader_result, library,
                             chem_comp_xml, depictions,
                             pubchem_templates, logger):
    """Process PDB-CCD component

    Args:
        args (ArgumentParser): Application arguments
        ccd_reader_result (pdbeccdutils.core.CCDReaderResult):
            pdbeccdutils parser output.
        library (pdbeccdutils.core.FragmentLibrary): fragment library
        chem_comp_xml (xml.etree.ElementTree.Element): root of the
            PDB-CCD XML representation
        depictions (pdbeccdutils.utils.DepictionManager): Object to take
            care of the pretty depictions.
        pubchem_templates (pdbeccdutils.utils.PubChemDownloader):
            Object to take care of the puchem templates fetching.
        logger (logging.loger): Application log
    """
    ccd_id = ccd_reader_result.component.id
    parent_dir = os.path.join(args.output_dir, ccd_id[0], ccd_id)
    os.makedirs(parent_dir)
    ideal_conformer = ConformerType.Ideal
    component = ccd_reader_result.component

    # check parsing and conformer degeneration
    issues = check_component_parsing(ccd_reader_result)
    structure_check = check_component_structure(ccd_reader_result.component)
    if len(structure_check) > 1:
        ideal_conformer = ConformerType.Computed
    issues += structure_check

    # download templates if the user wants them.
    if pubchem_templates is not None:
        issues += download_template(pubchem_templates, component)

    # search fragment library
    issues += search_fragment_library(component, library)

    # write out files
    issues += generate_depictions(component, depictions, parent_dir)
    issues += export_structure_formats(component, parent_dir, ideal_conformer)

    # get xml representation
    xml_repr = writer.to_xml_xml(ccd_reader_result.component)
    chem_comp_xml.append(xml_repr)

    # write log
    [logger.debug(f'{ccd_id} | {msg}') for msg in issues]


def check_component_parsing(ccd_reader_result):
    issues = []

    if len(ccd_reader_result.warnings) > 0:
        issues.append(f'warnings: {";".join(ccd_reader_result.warnings)}')

    if len(ccd_reader_result.errors) > 0:
        issues.append(f'errors: {";".join(ccd_reader_result.errors)}')

    if not ccd_reader_result.component.sanitize():
        issues.append('sanitization issue.')

    if not ccd_reader_result.component.inchikey_from_rdkit_matches_ccd():
        issues.append('inchikey mismatch.')

    return issues


def check_component_structure(component):
    issues = []
    if component.has_degenerated_conformer(ConformerType.Ideal):
        issues.append('has degenerated ideal coordinates.')
        result = component.compute_3d()
        if not result:
            issues.append('3D conformation could not be generated.')

    return issues


def download_template(pubchem_templates, component):
    component_downloaded = pubchem_templates.download_template(component)
    if component_downloaded:
        return ['downloaded new pubchem template.']

    return []


def generate_depictions(component, depictions, parent_dir):
    issues = []
    depiction_result = component.compute_2d(depictions)
    ccd_id = component.id

    if depiction_result.source == DepictionSource.Failed:
        issues.append('failed to generate 2D image.')
    else:
        if depiction_result.score > 0.99:
            issues.append('collision free image could not be generated.')

    component.export_2d_svg(os.path.join(parent_dir, f'{ccd_id}_small.svg'), width=100)
    component.export_2d_svg(os.path.join(parent_dir, f'{ccd_id}_medium.svg'), width=300)
    component.export_2d_svg(os.path.join(parent_dir, f'{ccd_id}_large.svg'), width=450)
    component.export_2d_svg(os.path.join(parent_dir, f'{ccd_id}_large_names.svg'), width=450, names=True)

    return issues


def search_fragment_library(component, library):
    matches = component.library_search(library)

    if matches > 0:
        return [f'{matches} matches found in the library `{library.name}`.']
    else:
        return []


def export_structure_formats(component, parent_dir, ideal_conformer):
    """Writes out component in a different formats as required for the
    PDBeChem FTP area

    Args:
        component (pdbeccdutils.core.Component): Component being processed.
        parent_dir (str): Working directory.
        ideal_conformer (pdbeccdutils.core.ConformerType): ConformerType
            to be used for ideal coordinates.

    Returns:
        (list of str): Encountered issues.
    """
    issues = []

    issues += write_molecule(os.path.join(parent_dir, f'{component.id}_model.sdf'), component, False, ConformerType.Model)
    issues += write_molecule(os.path.join(parent_dir, f'{component.id}_ideal.sdf'), component, False, ideal_conformer)
    issues += write_molecule(os.path.join(parent_dir, f'{component.id}_no_h_model.sdf'), component, True, ConformerType.Model)
    issues += write_molecule(os.path.join(parent_dir, f'{component.id}_no_h_ideal.sdf'), component, True, ideal_conformer)
    issues += write_molecule(os.path.join(parent_dir, f'{component.id}_model.pdb'), component, False, ConformerType.Model)
    issues += write_molecule(os.path.join(parent_dir, f'{component.id}_ideal.pdb'), component, False, ideal_conformer)
    issues += write_molecule(os.path.join(parent_dir, f'{component.id}_ideal.pdb'), component, False, ideal_conformer)
    issues += write_molecule(os.path.join(parent_dir, f'{component.id}.cml'), component, True, ConformerType.Model)
    issues += write_molecule(os.path.join(parent_dir, f'{component.id}.cif'), component, False, ConformerType.Model)

    return issues


def write_molecule(path, component, remove_hydrogens, conformer_type):
    """Write out deemed structure.

    Args:
        path (): [description]
        component (pdbeccdutils.core.Component): Component to be written.
        remove_hydrogens (bool): Whether or not Hs will be removed.
        conformer_type (pdbeccdutils.core.Component): Conformer to be written.

    Returns:
        (list of str): encountered issues
    """
    try:
        writer.write_molecule(path, component, remove_hydrogens, conformer_type)
        return []
    except Exception:
        with open(path, 'w') as f:
            f.write('')
        return [f'{path} could not be writter.']


def write_components_xml(args, chem_comp_xml):
    """Write out XML representation of the components.cif file

    Args:
        args (ArgumentParser): Application arguments
        chem_comp_xml (xml.etree.ElementTree.Element): xml object with
            the data.
    """
    xml_str = ET.tostring(chem_comp_xml, encoding='utf-8', method='xml')
    pretty = minidom.parseString(xml_str)

    with open(os.path.join(args.output_dir, "chem_comp_list.xml"), 'w') as f:
        f.write(pretty.toprettyxml(indent="  "))


def main():
    """Runs the PDBeChem pipeline
    """
    parser = create_parser()
    args = parser.parse_args()

    check_args(args)
    pdbechem_pipeline(args)
