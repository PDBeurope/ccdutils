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
import urllib.request
import xml.etree.ElementTree as ET
from collections import namedtuple

import rdkit
from mmCif.mmcifIO import MMCIF2Dict

import pdbeccdutils
import arpeggio
from arpeggio.core import InteractionComplex
from pdbeccdutils.core import CCDUtilsError, ccd_reader, structure_writer, ConformerType
from pdbeccdutils.utils import DepictionManager

IGNORED_LIGANDS = ('HOH', 'SO4')

# region bound molecules parsing


class Graph:
    """Data structure to house interactions. Can be used to infere
    connected components.
    """

    def __init__(self, nodes=set(), edges=set()):
        self.nodes = nodes
        self.edges = edges

    def __str__(self):
        return '-'.join(map(lambda l: str(l), self.nodes))

    def to_arpeggio_format(self):
        return list(map(lambda l: l.to_arpeggio(), self.nodes))

    def add_edge(self, e):
        self.edges.add(e)

    def add_node(self, n):
        self.nodes.add(n)

    def to_dict(self):
        nodes = sorted(list(self.nodes), key=lambda l: (int(l.res_id), l.chain))

        results_bag = {'ligands': [], 'links': []}
        results_bag['ligands'] = list(map(lambda l: l.to_dict(), nodes))
        results_bag['links'] = list(map(lambda l: (nodes.index(l[0]), nodes.index(l[1])), self.edges))

        return results_bag

    def pop_component(self):
        """Get connected component.

        Returns:
            Graph: connected component representing a bound molecule.
        """
        if len(self.nodes) == 0:
            return None
        visited_edges = set()
        visited_nodes = set()
        stack = set([self.nodes.pop()])

        while len(stack) > 0:
            processed_node = stack.pop()
            visited_nodes.add(processed_node)
            edges = list(filter(lambda l: processed_node in l.nodes, self.edges))

            for e in edges:
                visited_edges.add(e)
                self.edges.remove(e)
                stack.add(e.get_other(processed_node))

        return Graph(nodes=visited_nodes, edges=visited_edges)


class Node:
    """Represents a single ligand.
    """

    def __init__(self, name, chain, res_id):
        self.name = name
        self.chain = chain
        self.res_id = res_id

    def __eq__(self, other):
        if self.name != other.name:
            return False
        if self.chain != other.chain:
            return False
        if self.res_id != other.res_id:
            return False

        return True

    def to_dict(self):
        return {
            'label_comp_id': self.name,
            'auth_asym_id': self.chain,
            'auth_seq_id': self.res_id
        }

    def __hash__(self):
        return hash(self.chain + self.res_id + self.name)

    def __str__(self):
        return f'/{self.name}/{self.res_id}/{self.chain}/'

    def to_arpeggio(self):
        return f'/{self.chain}/{self.res_id}/'


class Edge:
    """Class representing a bond parsed fromt he _struct_conn namespace
    """

    def __init__(self, a, b):
        self.nodes = (a, b)

    def __eq__(self, other):
        return self[0] in other.nodes and self[1] in other.nodes

    def __hash__(self):
        return hash(self[0]) * hash(self[1])

    def get_other(self, a):
        """Get the other atom from the bond.

        Args:
            a Node: A node whose partner we are looking for

        Returns:
            Node: The other atom which is a part of the bond
        """
        if a in self.nodes:
            return self[1] if self[0] == a else self[0]
        else:
            return None

    def __str__(self):
        return f'{self[0]} - {self[1]}'

    def __getitem__(self, i):
        """Indexer method to access n-th atom in the bond.

        Args:
            i (int): atom index

        Raises:
            AttributeError: if index is other than 0,1

        Returns:
            [type]: [description]
        """
        if i < 0 or i > 1:
            raise AttributeError('Bond has just two atoms')
        return self.nodes[i]


def __infer_bound_molecules(config, structure):
    """Identify bound molecules in the input protein structure.

    Args:
        config (ArgumentParser): Application configuration.
        structure (str): Path to the structure.

    Returns:
        [list of Graph]: Bound molecules found in the pdb entry.
    """
    bms = []
    temp = _parse_bound_molecules(structure)

    while True:
        component = temp.pop_component()
        if component is None:
            return bms
        else:
            bms.append(component)


def _parse_bound_molecules(path):
    """Parse information from the information about HETATMS from the
    `_pdbx_nonpoly_scheme` (or `_atom_sites` if the former one is not
    present) and connectivity among them from `_struct_conn`.

    Args:
        path (str): Path to the mmCIF structure

    Returns:
        Graph: All the bound molecules in a given entry.
    """
    def __parse_ligands_from_nonpoly_schema(schema):
        g = Graph()
        for i in range(len(schema['asym_id'])):
            n = Node(
                schema['auth_mon_id'][i],  # aka label_comp_id
                schema['pdb_strand_id'][i],  # aka auth_asym_id
                schema['pdb_seq_num'][i])  # aka auth_seq_id

            if n.name not in IGNORED_LIGANDS:
                g.add_node(n)

        return g

    def __parse_ligands_from_atom_sites(atom_sites):
        g = Graph()
        for i in range(len(atom_sites['id'])):
            if atom_sites['group_PDB'][i] == 'HETATM':
                n = Node(
                    atom_sites['label_comp_id'][i],  # aka label_comp_id
                    atom_sites['auth_asym_id'][i],  # aka auth_asym_id
                    atom_sites['auth_seq_id'][i])  # aka auth_seq_id

                if n.name not in IGNORED_LIGANDS:
                    g.add_node(n)

        return g

    def __add_connections(g, struct_conn):
        for i in range(len(struct_conn['id'])):
            a = filter(lambda l:
                       l.name == struct_conn['ptnr1_label_comp_id'][i] and
                       l.chain == struct_conn['ptnr1_auth_asym_id'][i] and
                       l.res_id == struct_conn['ptnr1_auth_seq_id'][i], g.nodes)
            a = next(a, None)
            b = filter(lambda l:
                       l.name == struct_conn['ptnr2_label_comp_id'][i] and
                       l.chain == struct_conn['ptnr2_auth_asym_id'][i] and
                       l.res_id == struct_conn['ptnr2_auth_seq_id'][i], g.nodes)
            b = next(b, None)

            if a is not None and b is not None:
                g.add_edge(Edge(a, b))

    parsed_str = list(MMCIF2Dict().parse(path).values())[0]

    graph = __parse_ligands_from_atom_sites(parsed_str['_atom_site'])
    # graph = __parse_ligands_from_nonpoly_schema(parsed_str['_pdbx_nonpoly_scheme'])
    if '_struct_conn' in parsed_str:
        __add_connections(graph, parsed_str['_struct_conn'])

    return graph


# endregion

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
    logger.debug(f'pdbeccdutils core v. {pdbeccdutils.__version__} RDKit v. {rdkit.__version__}, Arpeggio v. {arpeggio.__version__}')

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

    Args:
        config (Namedtuple): Application configuration from the json file.
        depictor (DepictionManager): Helper class to get nice 2D images.
        pdb (str): Pdb id.
    """
    result_bag = {'depictions': {}}
    wd = os.path.join(config.output_dir, pdb[1:3], pdb)
    os.makedirs(wd, exist_ok=True)

    source_path = os.path.join(config.input_dir, pdb[1:3], pdb, "clean_mmcif", f'{pdb}_updated.cif.gz')
    assembly_path = __get_cs_structure(config, wd, pdb)
    protonated_cif_path = os.path.join(wd, f'{pdb}_h.cif')
    protonated_pdb_path = os.path.join(wd, f'{pdb}_h.pdb')

    bound_molecules = __infer_bound_molecules(config, assembly_path)

    if len(bound_molecules) == 0:
        logging.debug('No bound molecules found. Skipping entry.')
        return
    else:
        logging.debug(f'{len(bound_molecules)} bound molecules found.')

    __add_hydrogens(config, wd, assembly_path, protonated_cif_path, protonated_pdb_path)

    i = 0
    for bm in bound_molecules:
        i += 1
        logging.debug(f'Bound molecule composition: {str(bm)}')

        for ligand in map(lambda l: l.name, bm.nodes):
            if ligand not in result_bag['depictions']:
                ligand_path = os.path.join(config.components, ligand[0], ligand, f'{ligand}.cif')
                ligand_layout = __create_ligand_layout(depictor, ligand_path)
                result_bag['depictions'][ligand] = ligand_layout[ligand]

        contacts = __run_arpeggio(config, protonated_cif_path, bm)
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
    # that the structure is NMR and assembly cannot be generated
    # using /full instead.
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


def __run_arpeggio(config, structure_path, bm):
    """Run Arpeggio algorithm to retrieve contacts between bound
    molecule and protein.

    Args:
        config (ArgumentParser): Application configuration.
        structure_path (str): Path to the *.cif file to be processed.
        bm (Graph): Bound molecule

    Returns:
        [dict of str]: Results of the arpeggio run with the contacts
            between selection and the protein.
    """
    selection = bm.to_arpeggio_format()

    i_complex = InteractionComplex(structure_path)
    i_complex.structure_checks()
    i_complex.address_ambiguities()
    i_complex.run_arpeggio(selection, 5.0, 0.1, False)

    return i_complex.get_contacts()


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
