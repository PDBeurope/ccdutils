import os
import sys
import pprint


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
    fragment_file_name = os.path.join(data_dir, 'smi.text2')
    try:
        open(fragment_file_name, 'r')
    except IOError as err:
        print "Error cannot open fragment file {0} error is '{1}'".format(fragment_file_name, err.strerror)
        sys.exit(1)

    lines = [line.strip("\n").strip("\t").strip("\r").replace("\n", "").replace("\t", "").replace("\r", "") for line in
             open(fragment_file_name)]
    smiles = dict(((smile, name) for name, smile in (line.split(':') for line in lines)))
    return smiles


def main():
    print "ccd_find_fragments"

    smiles_2_fragment_name = load_smiles_from_smi_text_file()
    pprint.pprint(smiles_2_fragment_name)


if __name__ == "__main__":
    main()