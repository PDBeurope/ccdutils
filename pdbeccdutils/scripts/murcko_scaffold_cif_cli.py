import argparse
import csv
import logging
import os
import sys
from functools import reduce

import rdkit

import pdbeccdutils
from pdbeccdutils.core import ccd_reader
from pdbeccdutils.core.models import ScaffoldingMethod


def _create_parser():
    """
    Sets up parse the command line options.

    Returns:
        argparse.Namespace: Parsed arguments.
    """
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    add_arg = parser.add_argument
    add_arg('components_cif', help='Input PDB-CCD components.cif file (must be specified)')

    add_arg('--output_dir', '-o', required=True,
            help='Create an output directory with output files to load in the database.')

    add_arg('--debug', action='store_true', help='Turn on debug message logging output')

    return parser


def _check_args(args):
    """Validate suplied arguments.

    Args:
        args (argparse.Namespace): Parsed arguments.
    """

    if not os.path.isfile(args.components_cif):
        print(f'{args.components_cif} does not exist', file=sys.stderr)
        sys.exit(os.EX_NOINPUT)

    if not os.path.isdir(args.output_dir):
        os.makedirs(args.output_dir, exist_ok=True)

    log = logging.getLogger(__name__)

    log.info('Settings:')
    for k, v in vars(args).items():
        log.info(f'{"":5s}{k:25s}{v}')


def _set_up_logger(args):
    """Set up application level logging.

    Args:
        args (argparse.Namespace): Parsed arguments.

    Returns:
        logging.Logger: Application log.
    """

    logger = logging.getLogger(__name__)

    level = logging.DEBUG if args.debug else logging.WARNING
    fmt = '[%(asctime)-15s]  %(message)s'
    logging.basicConfig(level=level, format=fmt, datefmt='%a, %d %b %Y %H:%M:%S')

    logger.debug(f'PDBe Murcko scaffold pipeline using:')
    logger.debug(f'pdbeccdutils core v. {pdbeccdutils.__version__} RDKit v. {rdkit.__version__}')


def read_component_cif(args):
    """Reads component cif file and creates a generator.

    Args:
        args (argparse.Namespace): Parsed arguments.
    """
    dictionary = ccd_reader.read_pdb_components_file(args.components_cif)
    for v in dictionary.values():
        yield v.component


def write_csv(scaffolds, filename):
    """Write scaffolds in the CSV format.

    Args:
        scaffolds (dict of str: SubstructureMapping): All the calculated scaffolds.
        filename (str): Path to the CSV file destination.
    """
    # we need header of the csv file
    with open(filename, "w") as f:
        fileWriter = csv.writer(f, delimiter=',', quoting=csv.QUOTE_MINIMAL)
        fileWriter.writerow(['het_id', 'murcko_smiles', 'murcko_atom_mapping', 'murcko_generic_smiles', 'murcko_generic_atom_mapping'])
        for k, v in scaffolds.items():
            murcko = [i for i in v if i.name == 'MurckoScaffold']
            generic = [i for i in v if i.name == 'MurckoGeneric']            
            
            murcko_smiles = murcko[0].smiles if murcko else ''
            murcko_mapping = reduce(lambda l, m: f"{str(l)};{str(m)}", murcko[0].mappings)
            
            generic_smiles = generic[0].smiles if generic else ''
            generic_mapping = reduce(lambda l,m: f"{str(l)};{str(m)}", generic[0].mappings)

            fileWriter.writerow([k, murcko_smiles, murcko_mapping, generic_smiles, generic_mapping])


def calculate_scaffold(args):
    """Generate RDKit murcko smile and murcko generic smile for
    each component and write out csv file for FTP area.

    Args:
        args (argparse.Namespace): Parsed arguments.
    """
    logger = logging.getLogger(__name__)
    murcko_scaffolds = {}
    component_generator = read_component_cif(args)

    for component in component_generator:
        try:
            component.get_scaffolds(ScaffoldingMethod.MurckoScaffold)
            component.get_scaffolds(ScaffoldingMethod.MurckoGeneric)

            scaffolds = component.scaffolds

        except Exception as e:
            logger.error(f'{component.id} | FAILED with {str(e)}.')
            continue

        murcko_scaffolds[component.id] = scaffolds

        murcko_scaffolds_count = sum(i.name == ScaffoldingMethod.MurckoScaffold.name for i in scaffolds)
        murcko_generic_scaffolds_count = sum(i.name == ScaffoldingMethod.MurckoGeneric.name for i in scaffolds)

        logger.info((f'{component.id} | '
                     f'Murcko scaffolds: {murcko_scaffolds_count}, '
                     f'Murcko generic scaffolds: {murcko_generic_scaffolds_count}'))

    write_csv(murcko_scaffolds, os.path.join(args.output_dir, 'scaffold_murcko_ccd.csv'))


def main():
    """Runs the Murcko scaffolds pipeline
    """
    parser = _create_parser()
    args = parser.parse_args()

    _set_up_logger(args)
    _check_args(args)

    calculate_scaffold(args)
