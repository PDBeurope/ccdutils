#!/usr/bin/env python
# temporary jiffy to calculate number of rotatable bonds for each chemical component -
# use SMILES parsed from  http://ftp.ebi.ac.uk/pub/databases/msd/pdbechem/chem.xml
# write to csv with -1
from rdkit import Chem
from rdkit.Chem import Descriptors
try:
    import xml.etree.cElementTree as ETree
except ImportError:
    import xml.etree.ElementTree as ETree
tree = ETree.ElementTree(file='/nfs/ftp/pub/databases/msd/pdbechem/chem.xml')
root = tree.getroot()
for chem_comp in root:
    comp_id = chem_comp.find('id').text
    smiles = chem_comp.find('stereoSmiles')
    nrb = -1  # means there is an error
    if smiles is not None:
        smiles = smiles.text
        rdkit_mol = Chem.MolFromSmiles(smiles)
        if rdkit_mol is not None:
            nrb = Descriptors.NumRotatableBonds(Chem.MolFromSmiles(smiles))
    print '{},{}'.format(comp_id, nrb)
