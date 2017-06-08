#!/usr/bin/env python
from rdkit import Chem
from rdkit.Chem import Descriptors
try:
    import xml.etree.cElementTree as ETree
except ImportError:
    import xml.etree.ElementTree as ETree
tree = ETree.ElementTree(file='temp_sample_chem_comp.xml')
root = tree.getroot()
for chem_comp in root:
    comp_id = chem_comp.find('id').text
    smiles = chem_comp.find('stereoSmiles').text
    nrb = Descriptors.NumRotatableBonds(Chem.MolFromSmiles(smiles))
    print 'debug id: {} smiles: {} nrb: {}'.format(comp_id, smiles, nrb)
