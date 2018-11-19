import pytest

from pdbeccdutils.core import ccd_reader
from pdbeccdutils.tests.tst_utilities import cif_filename

test_inputs = {
    'NAG': {
        'logp': -3.078,
        'heavy_atom_count': 15,
        'numH_acceptors': 6,
        'numH_donors': 5,
        'num_rotable_bonds': 7,
        'rings_count': 1,
        'TPSA': 119.250,
        'molwt': 221.209
    },
    'UNL': {
        'logp': 0.0,
        'heavy_atom_count': 0,
        'numH_acceptors': 0,
        'numH_donors': 0,
        'num_rotable_bonds': 0,
        'rings_count': 0,
        'TPSA': 0.0,
        'molwt': 0.0
    }
}



class TestPropertyCalculation:

    @staticmethod
    @pytest.mark.parametrize('id', test_inputs)
    def test_properties(id):
        component = ccd_reader.read_pdb_cif_file(cif_filename(id)).component

        assert test_inputs[id]['logp'] == round(component.properties.logP, 3)
        assert test_inputs[id]['heavy_atom_count'] == component.properties.heavy_atom_count
        assert test_inputs[id]['numH_acceptors'] == component.properties.numH_acceptors
        assert test_inputs[id]['numH_donors'] == component.properties.numH_donors
        assert test_inputs[id]['num_rotable_bonds'] == component.properties.num_rotable_bonds
        assert test_inputs[id]['rings_count'] == component.properties.ring_count
        assert test_inputs[id]['TPSA'] == round(component.properties.TPSA, 3)
        assert test_inputs[id]['molwt'] == round(component.properties.molwt, 3)