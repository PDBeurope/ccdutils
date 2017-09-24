# software from PDBe: Protein Data Bank in Europe; http://pdbe.org
#
# Copyright 2017 EMBL - European Bioinformatics Institute
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
#
"""
process a PDB-CCD mmcif file.
"""
import argparse
import logging
from argparse import RawTextHelpFormatter

from pdbeccdutils.fragment_library import FragmentLibrary
from pdbeccdutils.pdb_chemical_components_rdkit import PdbChemicalComponentsRDKit


def __parse_command_line_args():
    """
    Sets up and parses the command line options.

    Returns:
         the arguments name space
    """
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=RawTextHelpFormatter)
    parser.add_argument('CIF', help='input PDB-CCD mmcif file (must be specified)')
    parser.add_argument('--library',
                        help='use this fragment library in place of the one supplied with the code.')
    parser.add_argument('--svg', help='write an 2D image of the compound to svg file')
    atom_labels_parser = parser.add_mutually_exclusive_group(required=False)
    atom_labels_parser.add_argument('--atom_labels', dest='atom_labels', action='store_true',
                                    help='turn on atom labels for 2D image.')
    atom_labels_parser.add_argument('--no_atom_labels', dest='atom_labels', action='store_false',
                                    help='turn off atom labels for 2D image (the default).')
    parser.set_defaults(atom_labels=False)
    parser.add_argument('--debug', action='store_true', help='turn on debug message logging output')
    return parser.parse_args()


def main():
    logger = logging.getLogger(' ')
    args = __parse_command_line_args()
    if args.debug:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)
    cif_file = args.CIF
    svg_file = args.svg
    atom_labels = args.atom_labels
    library = args.library
    logger.debug('CIF is {}'.format(cif_file))

    try:
        pdb_cc_rdkit = PdbChemicalComponentsRDKit(file_name=cif_file)
    except IOError as e_detail:
        raise SystemExit(e_detail)
    chem_comp_id = pdb_cc_rdkit.chem_comp_id
    logger.debug('chem_comp_id={}'.format(chem_comp_id))

    frag_lib = FragmentLibrary(override_fragment_library_file_path=library)
    fragments = frag_lib.fragments_for_pdb_chemical_components_rdkit(pdb_cc_rdkit)
    logger.info(' fragments:')
    for (name, list_of_atom_list) in fragments.items():
        logger.info('    {} occurs {} times:'.format(name, len(list_of_atom_list)))
        for atoms in list_of_atom_list:
            logger.info('         {}'.format(' '.join(atoms)))

    if svg_file is not None:
        pdb_cc_rdkit.image_file_or_string(file_name=svg_file, atom_labels=atom_labels)


if __name__ == "__main__":
    main()
