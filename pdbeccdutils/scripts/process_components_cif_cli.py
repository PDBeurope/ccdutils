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
import json
import logging
import os
import sys
import xml.etree.ElementTree as ET
from xml.dom import minidom
from typing import Dict, Any

import rdkit

import pdbeccdutils
from pdbeccdutils.core import ccd_reader, ccd_writer
from pdbeccdutils.core.exceptions import CCDUtilsError
from pdbeccdutils.core.component import Component
from pdbeccdutils.core.depictions import DepictionManager
from pdbeccdutils.core.fragment_library import FragmentLibrary
from pdbeccdutils.core.models import ConformerType, DepictionSource
from pdbeccdutils.utils import PubChemDownloader, config


# region pre-light tasks
def create_parser():
    """
    Sets up parse the command line options.

    Returns:
         argparse.Namespace parser
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
        args (argparse.Namespace): an argparse namespace containing the
            required arguments
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


def _set_up_logger(args):
    """Set up application level logging.

    Args:
        args (argparse.Namespace): Parsed arguments.
    """

    logger = logging.getLogger(__name__)
    level = logging.DEBUG if args.debug else logging.ERROR
    format = '[%(asctime)-15s]  %(message)s'
    logging.basicConfig(level=level, format=format, datefmt='%a, %d %b %Y %H:%M:%S')
    logger.info(f'PDBeChem pipeline using:')
    logger.info(f'pdbeccdutils core v. {pdbeccdutils.__version__}, RDKit v. {rdkit.__version__}')


def _log_settings(args):
    """Compose initial message about application settions.

    Args:
        args (argparse.Namespace): Application arguments.

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
# endregion


def pdbechem_pipeline(args):
    """
    Processes components.cif for PDBeChem type output

    Args:
        args (argparse.Namespace): an argparse namespace containing the
            required arguments
    """

    # set up logging
    _set_up_logger(args)
    logger = logging.getLogger(__name__)
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
            logger.info(f'{key} | processing...')
            process_single_component(args, ccd_reader_result, fragment_library, chem_comp_xml,
                                     depictions, pubchem_templates)
        except Exception as e:
            logger.error(f'{key} | FAILURE {str(e)}.')
        counter -= 1

        if counter == 0:
            break

    # write components.cif wide data
    write_components_xml(args, chem_comp_xml)
    with open(os.path.join(args.output_dir, 'chem_comp.list'), 'w') as f:
        f.write("\n".join(ids))


def process_single_component(args, ccd_reader_result, library,
                             chem_comp_xml, depictions,
                             pubchem_templates):
    """Process PDB-CCD component

    Args:
        args (argparse.Namespace): Application arguments
        ccd_reader_result (CCDReaderResult):
            pdbeccdutils parser output.
        library (FragmentLibrary): fragment library
        chem_comp_xml (xml.etree.ElementTree.Element): root of the
            PDB-CCD XML representation
        depictions (DepictionManager): Object to take
            care of the pretty depictions.
        pubchem_templates (PubChemDownloader):
            Object to take care of the puchem templates fetching.
    """
    ccd_id = ccd_reader_result.component.id
    parent_dir = os.path.join(args.output_dir, ccd_id[0], ccd_id)
    os.makedirs(parent_dir, exist_ok=True)
    ideal_conformer = ConformerType.Ideal
    component = ccd_reader_result.component
    json_output = {'het_code': ccd_id}

    # check parsing and conformer degeneration
    issues = check_component_parsing(ccd_reader_result)
    structure_check = check_component_structure(ccd_reader_result.component)
    if len(structure_check) == 1:
        ideal_conformer = ConformerType.Computed
    issues += structure_check

    # download templates if the user wants them.
    if pubchem_templates is not None:
        issues += download_template(pubchem_templates, component)

    # search fragment library
    issues += search_fragment_library(component, library, json_output)

    # get scaffolds
    issues += compute_component_scaffolds(component, json_output)
    # write out files
    issues += generate_depictions(component, depictions, parent_dir)
    issues += export_structure_formats(component, parent_dir, ideal_conformer)

    # get xml representation
    xml_repr = ccd_writer.to_xml_xml(ccd_reader_result.component)
    chem_comp_xml.append(xml_repr)

    with open(os.path.join(parent_dir, f'{ccd_id}_substructures.json'), 'w') as f:
        json.dump(json_output, f, sort_keys=True, indent=4)

    # write log
    logger = logging.getLogger(__name__)
    [logger.debug(f'{ccd_id} | {msg}') for msg in issues]


def check_component_parsing(ccd_reader_result):
    """Checks components parsing and highlights issues encountered with
    the molecule: errors/warnings during the parsing process,
    unrecoverable sanitization issues, inchikey mismatch between what
    was in the source file and is reproduced by rdkit.

    Args:
        ccd_reader_result (CCDReaderResult):
            Output of the parsing process.

    Returns:
        (list of str): possible issues encountered.
    """
    issues = []

    if len(ccd_reader_result.warnings) > 0:
        issues.append(f'warnings: {";".join(ccd_reader_result.warnings)}')

    if len(ccd_reader_result.errors) > 0:
        issues.append(f'errors: {";".join(ccd_reader_result.errors)}')

    if not ccd_reader_result.component.sanitized:
        issues.append('sanitization issue.')

    if not ccd_reader_result.component.inchikey_from_rdkit_matches_ccd():
        issues.append('inchikey mismatch.')

    return issues


def check_component_structure(component: Component):
    """Checks whether or not the component has degenerated ideal
    coordinates. If so, new conformer is attempted to be generated.

    Args:
        component (Component): Component to be
            processed.

    Returns:
        list of str: possible issues encountered.
    """
    issues = []
    if component.has_degenerated_conformer(ConformerType.Ideal):
        issues.append('has degenerated ideal coordinates.')
        result = component.compute_3d()
        if not result:
            issues.append('error in generating 3D conformation.')

    return issues


def download_template(pubchem_templates: PubChemDownloader, component: Component):
    """Attempts to download a pubchem template for the given component

    Args:
        pubchem_templates (PubChemDownloader):
            Pubchem downloader instance.
        component (Component): Component to be used.

    Returns:
        (list of str): information whether or not the new template has
            been downloaded.
    """
    component_downloaded = pubchem_templates.process_template(component)
    if component_downloaded:
        return ['downloaded new pubchem template.']

    return []


def generate_depictions(component: Component, depictions: DepictionManager, parent_dir: str):
    """Generate nice 2D depictions for the component. Presently depictions
    are generated in the following resolutions (100,200,300,400,500) with
    and without atom names.

    Args:
        component (Component): Component to be depicted.
        depictions (DepictionManager): Helper class
            to carry out depiction process.
        parent_dir (str): Where the depiction should be stored

    Returns:
        list of str: Possible issues encountered.
    """
    issues = []
    depiction_result = component.compute_2d(depictions)
    ccd_id = component.id

    if depiction_result.source == DepictionSource.Failed:
        issues.append('failed to generate 2D image.')
    else:
        if depiction_result.score > 0.99:
            issues.append('collision free image could not be generated.')
        issues.append(f'2D generated using {depiction_result.source.name} with score {depiction_result.score}.')

    for i in range(100, 600, 100):
        component.export_2d_svg(os.path.join(parent_dir, f'{ccd_id}_{i}.svg'), width=i)
        component.export_2d_svg(os.path.join(parent_dir, f'{ccd_id}_{i}_names.svg'), width=i, names=True)

    return issues


def search_fragment_library(component: Component, library: FragmentLibrary, json_output: Dict[str, Any]):
    """Search fragment library to find hits

    Args:
        component (Component): Component to be processed
        library (FragmentLibrary): Fragment library to be used.
        json_output (Dict[str, Any]): dictionary like structure with the
            results to be stored.

    Returns:
        list of str: info msg.
    """
    json_output['fragments'] = []
    matches = component.library_search(library)

    for k, v in component.fragments.items():
        json_output['fragments'].append({
            'name': k,
            'smiles': rdkit.Chem.MolToSmiles(library.library[k]),
            'mapping': v
        })

    if matches > 0:
        return [f'{matches} matches found in the library `{library.name}`.']
    else:
        return []


def compute_component_scaffolds(component: Component, json_output: Dict[str, Any]):
    """Compute scaffolds for a given component.

    Args:
        component (Component): Component to be processed
        json_output (Dict[str, Any]): dictionary like structure with the
            results to be stored.

    Returns:
        list of str: logging information.
    """
    try:
        scaffolds = component.get_scaffolds()
    except CCDUtilsError as e:
        return [str(e)]

    json_output['scaffolds'] = []
    for scaffold in scaffolds:
        atom_names = component.locate_fragment(scaffold)
        scaffold_atom_names = []
        for match in atom_names:
            scaffold_atom_names.append([i.GetProp('name') for i in match])

        json_output['scaffolds'].append({
            'smiles': rdkit.Chem.MolToSmiles(scaffold),
            'mapping': scaffold_atom_names
        })

    return [f'{len(scaffolds)} scaffold(s) were found.']


def export_structure_formats(component, parent_dir, ideal_conformer):
    """Writes out component in a different formats as required for the
    PDBeChem FTP area.

    Args:
        component (Component): Component being processed.
        parent_dir (str): Working directory.
        ideal_conformer (ConformerType): ConformerType
            to be used for ideal coordinates.

    Returns:
        list of str: Encountered issues.
    """
    issues = []

    issues += write_molecule(os.path.join(parent_dir, f'{component.id}_model.sdf'), component, False, ConformerType.Model)
    issues += write_molecule(os.path.join(parent_dir, f'{component.id}_ideal.sdf'), component, False, ideal_conformer)
    issues += write_molecule(os.path.join(parent_dir, f'{component.id}_ideal_alt.pdb'), component, True, ideal_conformer)
    issues += write_molecule(os.path.join(parent_dir, f'{component.id}_model_alt.pdb'), component, True, ConformerType.Model)
    issues += write_molecule(os.path.join(parent_dir, f'{component.id}_ideal.pdb'), component, False, ideal_conformer)
    issues += write_molecule(os.path.join(parent_dir, f'{component.id}_model.pdb'), component, False, ConformerType.Model)
    issues += write_molecule(os.path.join(parent_dir, f'{component.id}.cml'), component, False, ConformerType.Model)
    issues += write_molecule(os.path.join(parent_dir, f'{component.id}.cif'), component, False, ConformerType.Model)

    return issues


def write_molecule(path, component, alt_names, conformer_type):
    """Write out deemed structure.

    Args:
        path str: Path where the molecule will be stored.
        component (Component): Component to be written.
        alt_names (bool): Whether or not molecule will be written with
            alternate names.
        conformer_type (Component): Conformer to be written.

    Returns:
        list of str: encountered issues
    """
    try:
        ccd_writer.write_molecule(path, component, remove_hs=False, alt_names=alt_names,
                                  conf_type=conformer_type)
        return []
    except Exception:
        with open(path, 'w') as f:
            f.write('')
        return [f'error writing {path}.']


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
