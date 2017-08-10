#!/usr/bin/env python
# try reading in a cif file with 5 PDB CCD using the PDBeCIF parser - can we split it up?
# see issue #9
import mmCif.mmcifIO as mmcifIO
cif_parser = mmcifIO.CifFileReader(input='data', preserve_order=True)
cif_obj = cif_parser.read('tests/components_cif/components.cif.first_five_comps', output='cif_wrapper')
for data_block in list(cif_obj.values()):
    chem_comp = data_block._chem_comp
    id = chem_comp['id'][0]
    print id