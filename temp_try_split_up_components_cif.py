#!/usr/bin/env python
# try reading in a cif file with 5 PDB CCD using the PDBeCIF parser - can we split it up?
# see issue #9
import os
import mmCif.mmcifIO as mmcifIO
from mmCif import CIFWrapper
from utilities import file_name_in_tsts_out
from pdb_chemical_components import PdbChemicalComponents
cif_parser = mmcifIO.CifFileReader(input='data', preserve_order=True)
# test full components.cif by commenting/uncommenting next 2 lines:
test_file = 'tests/components_cif/components.cif.first_five_comps'
#test_file = '/homes/osmart/Downloads/components.cif'
cif_dictionary = cif_parser.read(test_file, output='cif_dictionary')
print('debug cif_obj is type {}'.format(type(cif_dictionary)))
for data_block_id, data_block in cif_dictionary.items():
    print('data_block_id = {}'.format(data_block_id))
    if data_block_id == '000':
        print('data_block = {}'.format(data_block))
    chem_comp_id = data_block['_chem_comp']['id']
    chem_comp_name = data_block['_chem_comp']['name']
    print('chem_comp_id={} chem_comp_name={}'.format(chem_comp_id, chem_comp_name))
    just_this_data_block = {}
    just_this_data_block[chem_comp_id] = data_block
    if len(cif_dictionary) < 10:
        output_cif_filename = file_name_in_tsts_out(chem_comp_id + '_split.cif')
    else:
        output_cif_filename = os.path.join('/var/tmp/test', chem_comp_id + '_split.cif')
    cfd = mmcifIO.CifFileWriter(output_cif_filename, preserve_order=True)
    cfd.write(just_this_data_block)
    print('have written file {}'.format(output_cif_filename))
    # create a pdbecif cij_obj from the dictionary
    mmcif_dict = just_this_data_block
    token_ordering = True
    # next line taken from CifFileReader method read
    cif_obj = dict(((block_id, CIFWrapper(block_data, data_id=block_id, preserve_token_order=token_ordering))
                    for block_id, block_data in list(mmcif_dict.items())))
    # can we parse it using PdbChemicalComponents?
    pdb_ccd = PdbChemicalComponents()
    pdb_ccd.read_ccd_from_pdbecif_cif_obj(cif_obj)
    print('pdb_ccd.chem_comp_name={}'.format(pdb_ccd.chem_comp_name))
    print('pdb_ccd.inchikey={}'.format(pdb_ccd.inchikey))
