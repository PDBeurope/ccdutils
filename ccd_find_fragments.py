import os
import sys
import pprint
import logging
from utilities import this_script_dir
from rdkit import Chem


def load_smiles_to_fragment_name_from_file( fragment_file_name):
    """
    loads the smiles to fragment name dictionary from the given file.

    Args:
        fragment_file_name (str): the fragment file name (normally smi.text in the data directory)

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

def cif_file_to_rwmol( cif_file):
    """

    Args:
        cif_file (str): filename of input PDB Chemical Component Definition CIF file containing the definition

    Returns:
        An RDKIT RWMol Mol object for the component
    """
    logging.debug('TODO write cif_file_to_rwmol')
    return None  # TODO write method!

def main():
    logging_level = logging.INFO
    logging_level = logging.DEBUG
    logging.basicConfig(format='%(levelname)s:%(message)s', level=logging_level)
    logging.info('ccd_find_fragments: start')

    data_dir = os.path.join(this_script_dir(), 'data')
    fragment_file_name = os.path.join(data_dir, 'smi.text')
    smiles_to_fragment_name = load_smiles_to_fragment_name_from_file(fragment_file_name)
    smiles_to_rdkit_mol = create_smiles_to_rdkit_mol(smiles_to_fragment_name.keys())

    if (logging_level == logging.DEBUG):
        logging.debug('check rdkit_mol conversion worked by converting back to SMILES')
        for smile in smiles_to_fragment_name:
            # try out inchi and inchikey
            test_inchi = Chem.inchi.MolToInchi(smiles_to_rdkit_mol[smile])
            test_inchikey = Chem.inchi.InchiToInchiKey(test_inchi)
            logging.debug('\tfragment {} original smiles:{} rdkit smiles:{} inchi:{} inchikey:{}'.
                          format(smiles_to_fragment_name[smile],
                                 smile,
                                 Chem.MolToSmiles(smiles_to_rdkit_mol[smile]),
                                 test_inchi, test_inchikey))

    cif_file = os.path.join(data_dir, '007.cif')  # temporary for testing.
    logging.debug('temporary for testing using cif file {}'.format(cif_file))
    my_rwmol = cif_file_to_rwmol(cif_file)
    logging.debug('temporary for testing my_rwmol returned as {}'.format(my_rwmol))

if __name__ == "__main__":
    main()