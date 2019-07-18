"""Test of PARITY method
"""

import pytest

from pdbeccdutils.computations.parity_method import compare_molecules
from pdbeccdutils.core import ccd_reader
from pdbeccdutils.tests.tst_utilities import cif_filename


class TestParity:

    @staticmethod
    @pytest.mark.parametrize('id_1,id_2,score', [
        ("MAN", 'NAG', 0.74),
        ("NAD", 'NAG', 0.09),
        ("000", '000', 1.00),
        ("SAC", 'GOL', 0.50),
        ("MAN", 'GLC', 1.00),
    ])
    def test_parity_method(id_1, id_2, score):
        c1 = ccd_reader.read_pdb_cif_file(cif_filename(id_1)).component
        c2 = ccd_reader.read_pdb_cif_file(cif_filename(id_2)).component

        parity = compare_molecules(c1.mol, c2.mol)

        assert round(parity.similarity_score, 2) == score
