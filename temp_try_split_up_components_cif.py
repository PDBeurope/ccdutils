#!/usr/bin/env python
# try reading in a cif file with 5 PDB CCD using the PDBeCIF parser - can we split it up?
# see issue #9
import os
import mmCif.mmcifIO as mmcifIO
from utilities import file_name_in_tsts_out
cif_parser = mmcifIO.CifFileReader(input='data', preserve_order=True)
# test full components.cif by commenting/uncommenting next 2 lines:
test_file = 'tests/components_cif/components.cif.first_five_comps'
test_file = '/homes/osmart/Downloads/components.cif'
cif_dictionary = cif_parser.read(test_file, output='cif_dictionary')
print('debug cif_obj is type {}'.format(type(cif_dictionary)))
for data_block_id, data_block in cif_dictionary.items():
    print 'data_block_id = ', data_block_id
    if data_block_id == '000':
        print 'data_block = ', data_block
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