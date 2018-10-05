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
import json
import logging
import os
import subprocess
import sys
import xml.etree.ElementTree as ET
from collections import namedtuple

import rdkit

import pdbeccdutils
from pdbeccdutils.core import (CCDUtilsError, ConformerType, ccd_reader,
                               structure_writer)
from pdbeccdutils.computations import DepictionManager
from pdbeccdutils.computations.interactions import ProtLigInteractions

# region logging


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
    logger.debug(f'PDBe protein-ligand interactions pipeline using:')
    logger.debug(f'pdbeccdutils core v. {pdbeccdutils.__version__} RDKit v. {rdkit.__version__}')

    return logger


def _log_settings(args):
    """Compose initial message about application settings.

    Args:
        args (ArgumentParser): Application arguments.

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
    settings += '{:29s}{:25s}{}\n'.format('', 'DEBUG', ('ON' if args.debug else 'OFF'))
    settings += '{:29s}{:25s}{}'.format('', 'running on ' + str(len(args.pdbs)) + ' structures.', '')

    return settings
# endregion


def create_parser():
    """
    Sets up a parser to get command line options.

    Returns:
         ArgumentParser parser
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

    }

    parser.add_argument('--config', required=True,
                        help=f'Configuration file in the following JSON format:\n {json.dumps(temp_config, sort_keys=True, indent=4)}')

    selection_group = parser.add_mutually_exclusive_group(required=True)
    selection_group.add_argument('-p', '--pdbs', type=str, nargs='+', help='List of pdb ids to be processed.')
    selection_group.add_argument('-sf', '--selection-file', type=str, help='Selections as above, but listed in a file.')
    parser.add_argument('--debug', action='store_true', help='Turn on debug message logging output')

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
        args (ArgumentParser): an argparse namespace containing the required arguments
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


def interactions_pipeline(args):
    """Given arguments compute protein-ligand interactions.

    Args:
        args (ArgumentParser): Command line arguments
    """
    logger = _set_up_logger(args)
    settings = _log_settings(args)
    logger.debug(settings)
    depictor = DepictionManager(pubchem_templates_path=args.config.pubchem_templates)

    for pdb in args.pdbs:
        logging.debug(f'Processing... {pdb}.')
        try:
            _process_single_structure(args.config, depictor, pdb)
        except Exception as e:
            logging.error(f'FAILED {pdb} | {str(e)}', file=sys.stderr)


def _process_single_structure(config, depictor, pdb):
    """Process a single pdb file and get composition of bound molecules
    along with the protein-ligand contacts and nice 2D depictions of
    ligands present in the bound molecules.

    TODO
    This is temporary implementation before the process is split into two

    Args:
        config (Namedtuple): Application configuration from the json file.
        depictor (DepictionManager): Helper class to get nice 2D images.
        pdb (str): Pdb id.
    """
    result_bag = {'depictions': {}}
    wd = os.path.join(config.output_dir, pdb[1:3], pdb)
    os.makedirs(wd, exist_ok=True)

    assembly_path = __get_cs_structure(config, wd, pdb)
    protonated_cif_path = os.path.join(wd, f'{pdb}_h.cif')
    protonated_pdb_path = os.path.join(wd, f'{pdb}_h.pdb')

    interactions = ProtLigInteractions(assembly_path, ['HOH', 'SO4'])  # encode res. names as a parameter!!

    if len(interactions.bound_molecules) == 0:
        logging.debug('No bound molecules found. Skipping entry.')
        return
    else:
        logging.debug(f'{len(interactions.bound_molecules)} bound molecule(s) found.')

    __add_hydrogens(config, wd, assembly_path, protonated_cif_path, protonated_pdb_path)

    i = 0
    interactions.path = protonated_cif_path
    for bm in interactions.bound_molecules:
        i += 1
        logging.debug(f'Bound molecule composition: {str(bm)}')

        for ligand in map(lambda l: l.name, bm.residues):
            if ligand not in result_bag['depictions']:
                ligand_path = os.path.join(config.components, ligand[0], ligand, f'{ligand}.cif')
                ligand_layout = __create_ligand_layout(depictor, ligand_path)
                result_bag['depictions'][ligand] = ligand_layout[ligand]

        print('getting contacts')
        print(bm.to_arpeggio_format())
        print(type(bm.to_arpeggio_format()))
        contacts = interactions.get_interaction(bm.to_arpeggio_format())
        result_bag[f'bm{i}'] = {'contacts': contacts}
        result_bag[f'bm{i}']['composition'] = bm.to_dict()

    with open(os.path.join(wd, 'contacts.json'), 'w') as f:
        json.dump(result_bag, f, sort_keys=True, indent=4)


def __get_assembly_id(path):
    """

    Args:
        path (str): Path to the assembly configuration in XML format.

    Raises:
        AttributeError: If the configuration file does not contain
            preferred id. This should not happen.

    Returns:
        str: preferred assembly id for a given PDB entry.
    """
    xml = ET.parse(path)

    for node in xml.getroot().iter('assembly'):
        if node.attrib['prefered'] == 'True':
            return node.attrib['id']

    raise AttributeError(f'Preferred id not found for the file {path}.')


def __get_cs_structure(config, wd, pdb):
    """Use instance of coordinate server to obtain biologicall assembly
    of the given protein.

    Args:
        config (ArgumentParser): Application configuration.
        wd (str): working directory
        pdb (str): 4-letter PDB id.

    Returns:
        str: Path to the stored assembly structure
    """
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
            logging.debug(f'{pdb} is likely to be NMR structure, generating /full structue.')
    else:
        logging.debug(f'Generated {pdb} with the assembly id {assembly_id}.')

    return assembly_str_path


def __run_cs(config, wd, cs_config_data):
    """Runs the coordinate server given the configuration file.

    Args:
        config (ArgumentParser): Application configuration.
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
        config (ArgumentParser): Application configuration.
        wd (str): Working directory
        structure (str): Path to the original struture
        protonated_cif (str): Path to the CIF protonated structure.
        protonated_pdb (str): Path to the PDB protonated structure.

    Raises:
        CCDUtilsError: If the ChimeraX fails or the protein structure
        does not exist.
    """
    logging.debug('Protonation started.')
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
        logging.debug('Protonation finished.')
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

    return structure_writer.to_json_dict(component, remove_hs=True, conf_type=ConformerType.Depiction)


def main():
    """Main method to execute protein-ligand interactions pipeline
    """
    parser = create_parser()
    args = parser.parse_args()

    check_args(args)
    interactions_pipeline(args)
