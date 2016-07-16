import os
import sys
import pprint
import logging


def load_smiles_from_smi_text_file():
    """
    loads the smiles to fragment name dictionary from the file smi.text in the data directory

    Returns:
         {} smiles to fragment name dict {str,str}

    Note:
        smi.txt has format
        cyclopropane:C1CC1
        phenyl:c1ccccc1
    """
    this_script_dir = os.path.dirname(os.path.abspath(__file__))
    data_dir = os.path.join(this_script_dir, 'data')
    fragment_file_name = os.path.join(data_dir, 'smi.text')
    try:
        fragment_file = open(fragment_file_name, 'r')
    except IOError as err:
        print "Error cannot open fragment file {0} error is '{1}'".format(fragment_file_name, err.strerror)
        sys.exit(1)
    lines =  fragment_file.read().splitlines()
    smiles_to_fragment_name = dict(((smile, name) for name, smile in (line.split(':') for line in lines)))
    return smiles_to_fragment_name


def main():
    print "ccd_find_fragments"

    logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.DEBUG)

    smiles_to_fragment_name = load_smiles_from_smi_text_file()
    logging.debug('Have loaded {} smiles, fragment names: '.format(len(smiles_to_fragment_name)))
    logging.debug(pprint.pformat(smiles_to_fragment_name))


if __name__ == "__main__":
    main()