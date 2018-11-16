import pytest

from rdkit import Chem

from pdbeccdutils.tests.tst_utilities import cif_filename
from pdbeccdutils.core import ccd_reader


"""
test of the property generation of hetcode NAG
"""

test_inputs = [('NAG', 'logp_v', -3.078)]
               # ('heavy_atom_count', '26'),
               # ('numH_acceptors', '4'),
               # ('numH_donors', '2'),
               # ('num_rotable_bonds', '7'),
               # ('rings_count', '3'),
               # ('TPSA', '66.75999999999999'),
               # ('molwt', '352.4299999999996')]

class TestIn:
    @staticmethod
    @pytest.mark.parametrize('id, logp_v, value', test_inputs)
    def test_property_calculation(id, logp_v, value):

        component = ccd_reader.read_pdb_cif_file(cif_filename(id)).component

        result = component.properties.logP
        assert result == value


