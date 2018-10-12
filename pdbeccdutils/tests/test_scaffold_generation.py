import pytest

from rdkit import Chem

from pdbeccdutils.tests.tst_utilities import cif_filename
from pdbeccdutils.core import ccd_reader


"""
Test of the Scaffolding generation
"""

test_inputs = [
    ('IBP', 'c1ccccc1'),
    ('NAG', 'C1CCOCC1'),
    ('VIA', 'O=c1[nH]c(-c2cc(S(=O)(=O)N3CCNCC3)ccc2)nc2c1[nH]nc2'),
    ('007', 'c1ccc(C2CCCC2)cc1'),
    ('BCD', 'C1CC2OCC1OC1CCC(CO1)OC1CCC(CO1)OC1CCC(CO1)OC1CCC(CO1)OC1CCC(CO1)OC1CCC(CO1)O2'),
    ('DMS', ''),
    ('EOH', '')
]



class TestIn:
    @staticmethod
    @pytest.mark.parametrize('id,smiles', test_inputs)
    def test_scaffold_present(id, smiles):

        component = ccd_reader.read_pdb_cif_file(cif_filename(id)).component
        component.sanitize()

        result = component.get_scaffolds()

        assert Chem.MolToSmiles(result[0])== smiles
