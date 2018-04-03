import argparse
from pdbeccdutils.core import structure_reader
from pdbeccdutils.utils import PubChemDownloader


def main():
    parser = argparse.ArgumentParser(description='PDBe downloader of pubchem depictions')
    parser.add_argument('-components', type=str, help='Path to the component library', required=True)
    parser.add_argument('-pubchem_templates', type=str, help='Path to the pubchem templates.',
                        required=True)

    config = parser.parse_args()
    PubChemDownloader(config.components, config.pubchem_templates).run()
