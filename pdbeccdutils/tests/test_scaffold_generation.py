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
    ('VIA', 'O=c1[nH]c(-c2cccc(S(=O)(=O)N3CCNCC3)c2)nc2cn[nH]c12'),
    ('007', 'c1ccc(C2CCCC2)cc1'),
    ('BCD', 'C1C[C@H]2OC[C@H]1O[C@H]1CC[C@@H](CO1)O[C@H]1CC[C@@H](CO1)O[C@H]1CC[C@@H](CO1)O[C@H]1CC[C@@H](CO1)O[C@H]1CC[C@@H](CO1)O[C@H]1CC[C@@H](CO1)O2'),
    ('DMS', ''),
    ('EOH', '')
]


class TestScaffold:
    @staticmethod
    @pytest.mark.parametrize('id,smiles', test_inputs)
    def test_scaffold_present(id, smiles):
        component = ccd_reader.read_pdb_cif_file(cif_filename(id)).component        
        result = component.get_scaffolds()

        assert Chem.MolToSmiles(result[0]) == smiles
