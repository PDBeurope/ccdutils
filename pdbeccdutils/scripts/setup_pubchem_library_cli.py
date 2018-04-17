import argparse
import os
import sys

from pdbeccdutils.core import structure_reader
from pdbeccdutils.utils import PubChemDownloader


def main():
    parser = argparse.ArgumentParser(description='PDBe downloader of pubchem depictions')
    parser.add_argument('-components_dir', type=str, default='',
                        help='Path to the directory with CCD files', required=False)
    parser.add_argument('-ccd', type=str, default='', help='Path to the CCD file', required=False)
    parser.add_argument('-pubchem_templates', type=str, help='Path to the pubchem templates.',
                        required=True)

    config = parser.parse_args()

    pubchem = PubChemDownloader(config.pubchem_templates)
    if os.path.isdir(config.components_dir):
        pubchem.update_ccd_dir(config.components_dir)
        sys.exit(os.EX_OK)
    elif os.path.isfile(config.ccd):
        pubchem.update_ccd_file(config.ccd)
        sys.exit(os.EX_OK)
    else:
        print('Either components_dir or components.cif file needs to be set as a source of data.')
        sys.exit(os.EX_USAGE)
