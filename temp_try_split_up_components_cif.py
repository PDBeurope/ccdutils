#!/usr/bin/env python
# try reading in a cif file with 5 PDB CCD using the PDBeCIF parser - can we split it up?
# see issue #9
import os
import mmCif.mmcifIO as mmcifIO
from utilities import file_name_in_tsts_out
from pdb_chemical_components_rdkit import PdbChemicalComponentsRDKit
cif_parser = mmcifIO.CifFileReader(input='data', preserve_order=True)
# test full components.cif by commenting/uncommenting next 2 lines:
test_file = 'tests/components_cif/components.cif.first_five_comps'
#test_file = '/homes/osmart/Downloads/components.cif'
cif_dictionary = cif_parser.read(test_file, output='cif_dictionary')
for data_block_id, data_block in cif_dictionary.items():
    print('data_block_id = {}'.format(data_block_id))
    # if data_block_id == '000':
    #     print('data_block = {}'.format(data_block))
    # chem_comp_id = data_block['_chem_comp']['id']
    # chem_comp_name = data_block['_chem_comp']['name']
    # print('chem_comp_id={} chem_comp_name={}'.format(chem_comp_id, chem_comp_name))
    individual_cif_dictionary = {}
    individual_cif_dictionary[data_block_id] = data_block
    if len(cif_dictionary) < 10:
        output_cif_filename = file_name_in_tsts_out(data_block_id + '_split.cif')
    else:
        output_cif_filename = os.path.join('/var/tmp/test', data_block_id  + '_split.cif')
    cfd = mmcifIO.CifFileWriter(output_cif_filename, preserve_order=True)
    cfd.write(individual_cif_dictionary)
    print('have written file {}'.format(output_cif_filename))
    pdb_ccd = PdbChemicalComponentsRDKit(cif_dictionary=individual_cif_dictionary)
    print('pdb_ccd.chem_comp_name={}'.format(pdb_ccd.chem_comp_name))
    print('pdb_ccd.inchikey={}'.format(pdb_ccd.inchikey))
    output_sdf = output_cif_filename + '.sdf'
    pdb_ccd.sdf_file_or_string(file_name=output_sdf)
    print('have written file {}'.format(output_sdf))
    print('\n')
