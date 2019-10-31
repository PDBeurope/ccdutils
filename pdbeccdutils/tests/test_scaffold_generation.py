import pytest
from rdkit import Chem

from pdbeccdutils.core import ccd_reader
from pdbeccdutils.tests.tst_utilities import cif_filename

test_inputs = [
    ('IBP', 'c1ccccc1'),
    ('NAG', 'C1CCOCC1'),
    ('VIA', 'O=c1[nH]c(-c2cccc(S(=O)(=O)N3CCNCC3)c2)nc2cn[nH]c12'),
    ('007', 'c1ccc(C2CCCC2)cc1'),
    ('BCD', 'C1C[C@@H]2OC[C@H]1O[C@@H]1CC[C@@H](CO1)O[C@@H]1CC[C@@H](CO1)O[C@@H]1CC[C@@H](CO1)O[C@@H]1CC[C@@H](CO1)O[C@@H]1CC[C@@H](CO1)O[C@@H]1CC[C@@H](CO1)O2')
]


class TestScaffold:
    """
    Test of the Scaffolding generation
    """
    @staticmethod
    @pytest.mark.parametrize('ccd_id,smiles', test_inputs)
    def test_scaffold_present(ccd_id, smiles):
        component = ccd_reader.read_pdb_cif_file(cif_filename(ccd_id)).component
        result = component.get_scaffolds()

        mol = result[0]

        assert Chem.MolToSmiles(mol) == smiles
        assert component.scaffolds
        assert component.scaffolds[0].smiles == smiles

    @staticmethod
    @pytest.mark.parametrize('ccd_id', ['EOH', 'DMS'])
    def test_scaffold_absent(ccd_id):
        component = ccd_reader.read_pdb_cif_file(cif_filename(ccd_id)).component
        result = component.get_scaffolds()

        assert result[0].GetNumAtoms() == 0
        assert not component.scaffolds
