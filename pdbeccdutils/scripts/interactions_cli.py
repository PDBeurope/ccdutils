#!/usr/bin/env python
# software from PDBe: Protein Data Bank in Europe; https://pdbe.org
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

import argparse
import gzip
import json
import logging
import os
import pathlib
import shutil
import subprocess
import sys
import xml.etree.ElementTree as ET
from collections import namedtuple

import rdkit

import pdbeccdutils
from pdbeccdutils.computations.interactions import ProtLigInteractions
from pdbeccdutils.core import ccd_reader, ccd_writer
from pdbeccdutils.core.depictions import DepictionManager
from pdbeccdutils.core.exceptions import CCDUtilsError
from pdbeccdutils.core.models import ConformerType

to_gzip = ('.json', '.pdb', '.cif')
# region logging


def _set_up_logger(args):
    """Set up application level logging.

    Args:
        args argparse.Namespace: Parsed arguments.

    Returns:
        logging.logger: Application log.
    """
    level = logging.DEBUG if args.debug else logging.WARNING
    format = '[%(asctime)-15s]  %(message)s'
    logging.basicConfig(level=level, format=format, datefmt='%a, %d %b %Y %H:%M:%S')

    logging.getLogger().disabled = True


def _log_settings(args):
    """Compose initial message about application settings.

    Args:
        args (argparse.Namespace): Application arguments.

    Returns:
        str: initial application message.
    """

    settings = 'Settings:\n'
    settings += '{:29s}{:25s}{}\n'.format('', 'generic_templates', args.config.generic_templates)
    settings += '{:29s}{:25s}{}\n'.format('', 'pubchem_templates', args.config.pubchem_templates)
    settings += '{:29s}{:25s}{}\n'.format('', 'components', args.config.components)
    settings += '{:29s}{:25s}{}\n'.format('', 'output_dir', args.config.output_dir)
    settings += '{:29s}{:25s}{}\n'.format('', 'input_dir', args.config.input_dir)
    settings += '{:29s}{:25s}{}\n'.format('', 'node', args.config.node)
    settings += '{:29s}{:25s}{}\n'.format('', 'coordinate_server', args.config.coordinate_server)
    settings += '{:29s}{:25s}{}\n'.format('', 'ChimeraX', args.config.chimerax)
    settings += '{:29s}{:25s}{}\n'.format('', 'no_header', ('ON' if args.no_header else 'OFF'))
    settings += '{:29s}{:25s}{}\n'.format('', 'gzip', ('ON' if args.gzip else 'OFF'))
    settings += '{:29s}{:25s}{}\n'.format('', 'DEBUG', ('ON' if args.debug else 'OFF'))
    settings += '{:29s}{:25s}{}\n'.format('', 'discarded_ligands', ','.join(args.config.discarded_ligands))
    settings += '{:29s}{:25s}{}'.format('', 'running on ' + str(len(args.pdbs)) + ' structures.', '')

    return settings
# endregion


def create_parser():
    """
    Sets up a parser to get command line options.

    Returns:
         argparse.Namespace parser
    """
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)

    temp_config = {
        "pubchem_templates": "path to the pubchem templates",
        "generic_templates": "path to the generic templates",
        "output_dir": "path to the output directory",
        "components": "Path to the components in the mmcif file (expected structure: mmcif/A/ATP/ATP.cif)",
        "node": "Path to the node.js binary",
        "coordinate_server": "Path to the coordinate server local.js binary",
        "input_dir": "path to the input data location (expected structure: ./tq/1tqn/assembly_generation/1tqn-assembly.xml and ./tq/1tqn/clean_mmcif/1tqn_updated.cif.gz)",
        "chimerax": "location of the ChimeraX binary",
        "discarded_ligands": "list of het codes to be discarded (e.g. HOH, SO4, etc.)"
    }

    parser.add_argument('--config', required=True,
                        help=f'Configuration file in the following JSON format:\n {json.dumps(temp_config, sort_keys=True, indent=4)}')

    selection_group = parser.add_mutually_exclusive_group(required=True)
    selection_group.add_argument('-p', '--pdbs', type=str, nargs='+', help='List of pdb ids to be processed.')
    selection_group.add_argument('-sf', '--selection-file', type=str, help='Selections as above, but listed in a file.')
    parser.add_argument('--debug', action='store_true', help='Turn on debug message logging output')
    parser.add_argument('--no_header', action='store_true', help='Turn off header information for the script.')
    parser.add_argument('--gzip', action='store_true', help='Whether or not should be all the files compressed to save space.')

    return parser


def __validate_file_exists(name, p):
    """Check if the file exists

    Args:
        name (str): Name of the parameter
        p (str): Value of the parameter
    """
    if not os.path.isfile(p):
        print(f'Supplied parameter {name} is an invalid path {p}.', file=sys.stderr)
        sys.exit(os.EX_NOINPUT)


def check_args(args):
    """Validate suplied arguments.

    Args:
        args (argparse.Namespace): an argparse namespace containing the required arguments
    """

    __validate_file_exists('Config', args.config)

    with open(args.config, 'r') as f:
        js = json.load(f)
        args.config = namedtuple('Config', js.keys())(*js.values())

    if args.selection_file is not None:
        with open(args.selection_file, 'r') as f:
            args.pdbs = list(map(lambda l: l.strip(), f.readlines()))

    __validate_file_exists('ChimeraX', args.config.chimerax)
    __validate_file_exists('coordinate_server', args.config.coordinate_server)


def interactions_pipeline(args) -> bool:
    """Given arguments compute protein-ligand interactions.

    Args:
        args (argparse.Namespace): Command line arguments

    Returns:
        bool: Whether or not errors have been encountered
    """
    logger = logging.getLogger(__name__)
    _set_up_logger(args)
    success = True

    if not args.no_header:
        logger.info(f'PDBe protein-ligand interactions pipeline using:')
        logger.info(f'pdbeccdutils core v. {pdbeccdutils.__version__} RDKit v. {rdkit.__version__}')
        settings = _log_settings(args)
        logger.info(settings)

    depictor = DepictionManager(pubchem_templates_path=args.config.pubchem_templates)

    for pdb in args.pdbs:
        logger.debug(f'Processing... {pdb}.')
        try:
            _process_single_structure(args, depictor, pdb)
        except Exception as e:
            success = False
            logger.error(f'{pdb} | FAILED')
            logger.exception(e)

    return success


def _process_single_structure(args, depictor, pdb):
    """Process a single pdb file and get composition of bound molecules
    along with the protein-ligand contacts and nice 2D depictions of
    ligands present in the bound molecules.

    TODO
    This is temporary implementation before the process is split into two

    Args:
        args (argparse.Namespace): Application configuration.
        depictor (DepictionManager): Helper class to get nice 2D images.
        pdb (str): Pdb id.
    """
    logger = logging.getLogger(__name__)
    config = args.config
    result_bag = {'entry': pdb, 'depictions': {}, 'boundMolecules': []}
    wd = os.path.join(config.output_dir, pdb[1:3], pdb)
    os.makedirs(wd, exist_ok=True)

    assembly_path = __get_cs_structure(config, wd, pdb)
    protonated_cif_path = os.path.join(wd, f'{pdb}_h.cif')
    protonated_pdb_path = os.path.join(wd, f'{pdb}_h.pdb')

    interactions = ProtLigInteractions(assembly_path, config.discarded_ligands)

    if len(interactions.bound_molecules) == 0:
        logger.debug('No bound molecules found. Skipping entry.')
        return
    else:
        logger.debug(f'{len(interactions.bound_molecules)} bound molecule(s) found.')

    __add_hydrogens(config, wd, assembly_path, protonated_cif_path, protonated_pdb_path)

    logger.debug('Initializing arpeggio.')
    i = 0
    interactions.path = protonated_cif_path
    interactions.initialize()
    for bm in interactions.bound_molecules:
        i += 1
        logger.debug(f'Running arpeggio (bm{i}) for: {str(bm)}')

        for ligand in map(lambda l: l.name, bm.residues):
            if ligand not in result_bag['depictions']:
                ligand_path = os.path.join(config.components, ligand[0], ligand, f'{ligand}.cif')
                ligand_layout = __create_ligand_layout(depictor, ligand_path)
                result_bag['depictions'][ligand] = ligand_layout[ligand]

        contacts = interactions.get_interaction(bm.to_arpeggio_format())
        contacts_filtered = list(filter(lambda l:
                                        l['interacting_entities'] in ('INTER', 'SELECTION_WATER'),
                                        contacts))

        result_bag['boundMolecules'].append({
            'id': f'bm{i}',
            'contacts': contacts_filtered,
            'composition': bm.to_dict()
        })

    with open(os.path.join(wd, 'contacts.json'), 'w') as f:
        json.dump(result_bag, f, sort_keys=True, indent=4)

    if args.gzip:
        __gzip_folder(wd)


def __get_assembly_id(path):
    """Get biological assembly id.

    Args:
        path (str): Path to the assembly configuration in XML format.

    Raises:
        AttributeError: If the configuration file does not contain
            preferred id. This should not happen.

    Returns:
        str: preferred assembly id for a given PDB entry.
    """
    if not os.path.isfile(path):
        raise AttributeError(f'Assembly configuration file {path} does not exist.')

    xml = ET.parse(path)

    for node in xml.getroot().iter('assembly'):
        if node.attrib['prefered'] == 'True':
            return node.attrib['id']

    raise AttributeError(f'Preferred id not found for the file {path}.')


def __get_cs_structure(config, wd, pdb):
    """Use instance of coordinate server to obtain biologicall assembly
    of the given protein.

    Args:
        config (argparse.Namespace): Application configuration.
        wd (str): working directory
        pdb (str): 4-letter PDB id.

    Returns:
        str: Path to the stored assembly structure
    """
    logger = logging.getLogger(__name__)
    input_str_path = os.path.join(config.input_dir, pdb[1:3], pdb, "clean_mmcif", f'{pdb}_updated.cif.gz')
    assembly_str_path = os.path.join(wd, f'{pdb}_assembly.cif')
    xml_path = os.path.join(config.input_dir,
                            pdb[1:3], pdb,
                            'assembly_generation',
                            f'{pdb}-assembly.xml')
    assembly_id = __get_assembly_id(xml_path)

    cs_config_data = [
        {
            "inputFilename": input_str_path,
            "outputFilename": assembly_str_path,
            "query": "assembly",
            "params": {
                "id": assembly_id
            }
        }
    ]

    cs_structure_generation_success = __run_cs(config, wd, cs_config_data)

    # If the structure generation fails it is likely because of the fact
    # that the structure is NMR and assembly cannot be generated.
    # Using /full instead.
    if not cs_structure_generation_success:
        cs_config_data[0]['params'] = {}
        cs_config_data[0]['query'] = 'full'

        cs_structure_generation_success = __run_cs(config, wd, cs_config_data)

        if not cs_structure_generation_success:
            raise CCDUtilsError('Structure generation using CoordinateServer failed.')
        else:
            logger.debug(f'{pdb} is likely to be NMR structure, generating /full structue.')
    else:
        logger.debug(f'Generated {pdb} with the assembly id {assembly_id}.')

    return assembly_str_path


def __run_cs(config, wd, cs_config_data):
    """Runs the coordinate server given the configuration file.

    Args:
        config (argparse.Namespace): Application configuration.
        wd (str): working directory.
        cs_config_data (dict of str): CoordinateServer config.

    Raises:
        CCDUtilsError: If the CoordinateServer run failed badly.

    Returns:
        bool: Whether or not the generation of the structure was succesfull.
    """
    cs_config_path = os.path.join(wd, 'cs_config.json')
    with open(cs_config_path, 'w') as f:
        json.dump(cs_config_data, f)

    try:
        out = subprocess.check_output([config.node, config.coordinate_server, cs_config_path],
                                      stderr=subprocess.STDOUT)
        out = out.decode('utf-8')

        return 'Failed' not in out

    except Exception:
        raise CCDUtilsError('Error while generating assembly file.')


def __add_hydrogens(config, wd, structure, protonated_cif, protonated_pdb):
    """Use ChimeraX to add hydrogens to the potein structure.
    The command being used is `addh hbonds true`.

    Args:
        config (argparse.Namespace): Application configuration.
        wd (str): Working directory
        structure (str): Path to the original struture
        protonated_cif (str): Path to the CIF protonated structure.
        protonated_pdb (str): Path to the PDB protonated structure.

    Raises:
        CCDUtilsError: If the ChimeraX fails or the protein structure
        does not exist.
    """
    logger = logging.getLogger(__name__)
    logger.debug('Protonation started.')
    cmd_file = os.path.join(wd, 'chimera_config.cxc')
    with open(cmd_file, 'w') as f:

        f.write("\n".join([
            f'open {structure}',
            'addh hbond true',
            f'save {protonated_cif} format mmcif',
            f'save {protonated_pdb} format pdb',
            'exit'
        ]))
    try:
        subprocess.call([config.chimerax, '--nogui', '--cmd', f'open {cmd_file}', '--silent'])
        logger.debug('Protonation finished.')
    except Exception:
        raise CCDUtilsError('Error while protonating file.')

    if not os.path.isfile(protonated_cif):
        raise CCDUtilsError('CIF protonated file was not created.')

    if not os.path.isfile(protonated_pdb):
        raise CCDUtilsError('PDB protonated file was not created.')


def __create_ligand_layout(depictor, ligand_path):
    """Retrieve 2D coordinates for a given ligand

    Args:
        depictor (DepictionManager): Helper class to aid nice 2d layouts.
        ligand_path (str): Path to the component.

    Returns:
        [dict of str]: json-like representation of the component layout.
    """

    component = ccd_reader.read_pdb_cif_file(ligand_path).component
    component.sanitize()
    component.compute_2d(depictor)

    return ccd_writer.to_json_dict(component, remove_hs=True, conf_type=ConformerType.Depiction)


def __gzip_folder(path):
    """Gzip all content of the calculation folder into gzip with the
    exception of config files.

    Args:
        path (str): Path to the working directory
    """
    for input_file in os.listdir(path):
        path_f = os.path.join(path, input_file)
        if pathlib.Path(path_f).suffix in to_gzip:
            with open(path_f, 'rb') as f_in:
                with gzip.open(f'{path_f}.gz', 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
            os.remove(path_f)


def main():
    """Main method to execute protein-ligand interactions pipeline
    """
    parser = create_parser()
    args = parser.parse_args()

    check_args(args)
    interactions_pipeline(args)
