import argparse
import csv
import logging
import os
import sys
from typing import List, NamedTuple

import rdkit
from rdkit import Chem

import pdbeccdutils
from pdbeccdutils.core import ccd_reader
from pdbeccdutils.core.models import ScaffoldingMethod

MurckoResults = NamedTuple('MurckoResults',
                           [('het_code', str),
                            ('smiles', List[str]),
                            ('atom_mapping', List[str]),
                            ('generic_smiles', str)])


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
    print
    if not os.path.isfile(args.components_cif):
        print(f'{args.components_cif} does not exist', file=sys.stderr)
        sys.exit(os.EX_NOINPUT)

    if not os.path.isdir(args.output_dir):
        os.makedirs(args.output_dir, exist_ok=True)


def _set_up_logger(args):
    """Set up application level logging.

    Args:
        args (argparse.Namespace): Parsed arguments.

    Returns:
        logging.Logger: Application log.
    """

    logger = logging.getLogger(__name__)

    level = logging.DEBUG if args.debug else logging.WARNING
    format = '[%(asctime)-15s]  %(message)s'
    logging.basicConfig(level=level, format=format, datefmt='%a, %d %b %Y %H:%M:%S')

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


def write_csv(data, filename):
    """Write scaffolds in the CSV format.

    Args:
        data (MurckoResults): All the calculated scaffolds
        filename (str): Path to the CSV file destination.
    """
    # we need header of the csv file
    with open(filename, "w") as f:
        fileWriter = csv.writer(f, delimiter=',', quoting=csv.QUOTE_MINIMAL)
        fileWriter.writerow(['het_id', 'murcko_smiles', 'murcko_atom_mapping', 'murcko_generic_smiles'])
        for row in data:
            fileWriter.writerow(row)


def calculate_scaffold(args):
    """Generate RDKit murcko smile and murcko generic smile for
    each component and write out csv file for FTP area.

    Args:
        args (argparse.Namespace): Parsed arguments.
    """
    logger = logging.getLogger(__name__)
    murcko_scaffolds = []
    component_generator = read_component_cif(args)

    for component in component_generator:
        try:
            list_of_scaffolds = component.get_scaffolds()
            temp = component.locate_fragment(list_of_scaffolds[0])  # get the atoms of the scaffolds
            murcko_atom_names = list(map(lambda l: l.GetProp('name'), temp[0]))  # get the names of the atoms
            scaffold_generic = component.get_scaffolds(ScaffoldingMethod.MurckoGeneric)  # get generic scaffold smiles

        except Exception as e:
            logger.error(f'{component.id} | FAILED with {str(e)}.')
            pass

        murcko_data = MurckoResults(het_code=component.id,
                                    smiles=[Chem.MolToSmiles(x) for x in list_of_scaffolds if x.GetNumAtoms() > 0],
                                    atom_mapping=murcko_atom_names,
                                    generic_smiles=[Chem.MolToSmiles(x) for x in scaffold_generic if x.GetNumAtoms() > 0])
        murcko_scaffolds.append(murcko_data)

        logger.info((f'{component.id} | '
                     f'Murcko: {len(murcko_data.smiles)}, '
                     f'Murcko Generic: {len(murcko_data.generic_smiles)}'))

    write_csv(murcko_scaffolds, os.path.join(args.output_dir, 'scaffold_murcko_ccd.csv'))


def main():
    """Runs the Murcko scaffolds pipeline
    """
    parser = _create_parser()
    args = parser.parse_args()
    print('Hi!')
    print(args)
    _set_up_logger(args)
    _check_args(args)

    calculate_scaffold(args)
