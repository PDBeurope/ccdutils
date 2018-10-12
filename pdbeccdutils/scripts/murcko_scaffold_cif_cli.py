from rdkit import Chem
import sys
from pdbeccdutils.core import ccd_reader
from pdbeccdutils.core.models import ScaffoldingMethod
import argparse
import os
from collections import namedtuple
import csv

MurckoScaffold_tuple = namedtuple('Scaffolds', 'het_code murcko_smiles murcko_atom_mapping murcko_generic_smiles')


def create_parser():
    """
    Sets up parse the command line options.

    Returns:
         ArgumentParser parser
    """
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    add_arg = parser.add_argument
    add_arg('components_cif', help='Input PDB-CCD components.cif file (must be specified)')

    add_arg('--output_dir', '-o', required=True,
            help='Create an output directory with output files to laod in the database')

    add_arg('--debug', action='store_true', help='Turn on debug message logging output')

    return parser


def check_args(args):
    """Validate suplied arguments.

    Args:
        args (ArgumentParser): an argparse namespace containing the required arguments
    """
    if not os.path.isfile(args.component_cif):
        print(f'{args.component_cif} does not exist', file=sys.stderr)
        sys.exit(os.EX_NOINPUT)

    if not os.path.isdir(args.output_dir):
        os.makedirs(args.output_dir)


def read_component_cif(args):
    """
    read component cif file and create a generator
    Args:
        args: an argparse namespace containing the required arguments

    Returns:  sanitize component as iterator

    """

    dictionary = ccd_reader.read_pdb_components_file(args.components_cif)
    for k, v in dictionary.items():
        rdkit_inv = v.component
        rdkit_inv.sanitize()
        yield rdkit_inv


def write_csv(data, filename):
    """
    create csv file from named tuple
    Args:
        data: namedtuple
        filename: csv filename

    Returns: csv file

    """
    # we need header of the csv file
    with open(filename, "w") as f:
        fileWriter = csv.writer(f, delimiter=',', quoting=csv.QUOTE_MINIMAL)
        fileWriter.writerow(['het_id', 'murcko_smiles', 'murcko_atom_mapping', 'murcko_generic_smiles'])
        for row in list(data):
            fileWriter.writerow(row)

    return filename


def calculate_scaffold(args):
    """
    generate rdkit murcko smile  and murcko generic smile for each component and write out csv file for ftp
    Args:
        args:

    Returns: csv file for each method, e.g.
    'TDP,['c1ncc(C[n+]2cscc2)cn1'],"[""N1'"", ""C2'"", ""N3'"", ""C4'"", ""C5'"", ""C6'"", 'C35', 'N3', 'C2', 'S1', 'C5'
        , 'C4']",['CC1CCC(CC2CCC(CCCC(C)(C)CC(C)(C)C)C2C)C(C)C1']'


    """

    murcko_scaffolds = []

    component_generator = read_component_cif(args)

    for component in component_generator:

        # this part calculate murcko scaffold and create the data structure

        try:

            list_of_scaffolds = component.get_scaffolds()  # get scaffolds
            temp = component.locate_fragment(list_of_scaffolds[0])  # get the atoms of the scaffolds
            murcko_atom_names = list(map(lambda l: l.GetProp('name'), temp[0]))  # get the names of the atoms
            scaffold_generic = component.get_scaffolds(ScaffoldingMethod.Murcko_generic) # get generic scaffold smiles

        except Exception:
            pass

        finally:
            murcko_data = MurckoScaffold_tuple(het_code=component.id,
                                            murcko_smiles=[Chem.MolToSmiles(y) for y in list_of_scaffolds],
                                            murcko_atom_mapping=murcko_atom_names,
                                            murcko_generic_smiles= [Chem.MolToSmiles(z) for z in scaffold_generic])

        murcko_scaffolds.append(murcko_data)


    write_csv(murcko_scaffolds, os.path.join(args.output_dir, 'scaffold_murcko_ccd.csv'))


def main():
    """Runs the PDBeChem pipeline
    """
    parser = create_parser()
    args = parser.parse_args()

    check_args(args)

    calculate_scaffold(args)
