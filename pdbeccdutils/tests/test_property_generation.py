import pytest

from pdbeccdutils.core import ccd_reader
from pdbeccdutils.tests.tst_utilities import cif_filename

test_inputs = {
    'NAG': {
        'logp': -3.078,
        'heavy_atom_count': 26
    }
}

# ('numH_acceptors', '4'),
# ('numH_donors', '2'),
# ('num_rotable_bonds', '7'),
# ('rings_count', '3'),
# ('TPSA', '66.75999999999999'),
# ('molwt', '352.4299999999996')]


class TestPropertyCalculation:

    @staticmethod
    @pytest.mark.parametrize('id', test_inputs)
    def test_properties(id):
        component = ccd_reader.read_pdb_cif_file(cif_filename(id)).component

        assert test_inputs[id]['logp'] == round(component.properties.logP, 3)
        # TODO other tests in a similar manner as before
