import os
import sys
import pprint
import logging
from rdkit import Chem

def this_script_dir():
    """
    gives the name of directory were script being run resides

    Returns:
        str the name of the directory is

    """
    return os.path.dirname(os.path.abspath(__file__))

def load_smiles_to_fragment_name_from_file( fragment_file_name):
    """
    loads the smiles to fragment name dictionary from the given file.

    Args:
        fragment_file_name: the fragment file name (normally smi.text in the data directory)

    Returns:
         {} smiles to fragment name dict {str,str}

    Note:
        smi.txt has format
        cyclopropane:C1CC1
        phenyl:c1ccccc1
    """
    try:
        fragment_file = open(fragment_file_name, 'r')
    except IOError as err:
        print "Error cannot open fragment file {0} error is '{1}'".format(fragment_file_name, err.strerror)
        sys.exit(1)
    smiles_to_fragment_name = {}
    lines =  fragment_file.read().splitlines()
    for line in lines:
        name, smile = line.split(':')
        smile = smile.replace('\t', '')  # take out tabs
        smiles_to_fragment_name[smile] = name
    number_of_entries = len(smiles_to_fragment_name)
    logging.debug('method load_smiles_from_smi_text_file:')
    logging.info('Have loaded smiles_to_fragment_name dictionary with {} entries from file {}'.
                 format(number_of_entries,fragment_file_name))
    logging.debug('\tdump entries:\n' + pprint.pformat(smiles_to_fragment_name))
    return smiles_to_fragment_name


def create_smiles_to_rdkit_mol( smiles_list):
    """
    creates a dictionary with an rdkit molecule for each SMILES string in the input list.

    Args:
        smiles_list: list of SMILES strings

    Returns:
        {} dictionary SMILES to rdkit molecule {string, Mol}

    """
    smiles_to_rdkit_mol = {}
    for smile in smiles_list:
        rdkit_mol = Chem.MolFromSmiles(smile)
        # Sameer had commented out next sanitize operation?
        Chem.SanitizeMol(rdkit_mol, sanitizeOps=Chem.SanitizeFlags.SANITIZE_ALL ^
                                                Chem.SanitizeFlags.SANITIZE_KEKULIZE ^
                                                Chem.SanitizeFlags.SANITIZE_SETAROMATICITY)
        smiles_to_rdkit_mol[smile] = rdkit_mol
    return smiles_to_rdkit_mol



def main():
    logging_level = logging.INFO
    logging_level = logging.DEBUG
    logging.basicConfig(format='%(levelname)s:%(message)s', level=logging_level)
    logging.info('ccd_find_fragments: start')

    data_dir = os.path.join(this_script_dir(), 'data')
    fragment_file_name = os.path.join(data_dir, 'smi.text')
    smiles_to_fragment_name = load_smiles_to_fragment_name_from_file(fragment_file_name)
    smiles_to_rdkit_mol = create_smiles_to_rdkit_mol(smiles_to_fragment_name.keys())

    logging.debug('check rdkit_mol conversion worked by converting back to SMILES')
    for smile in smiles_to_fragment_name:
        logging.debug('\tfragment {} original smiles:{} rdkit smiles:{}'.
                      format(smiles_to_fragment_name[smile],
                             smile,
                             Chem.MolToSmiles(smiles_to_rdkit_mol[smile])))

    #f_file = "d"

if __name__ == "__main__":
    main()