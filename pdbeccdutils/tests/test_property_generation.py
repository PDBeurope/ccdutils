import pytest

from pdbeccdutils.core import ccd_reader
from pdbeccdutils.tests.tst_utilities import cif_filename, unl_model

test_inputs = {
    'ATP': {
        'logp': -2.438,
        'heavy_atom_count': 31,
        'numH_acceptors': 18,
        'numH_donors': 7,
        'num_rotable_bonds': 15,
        'rings_count': 3,
        'TPSA': 279.130,
        'molwt': 506.996
    },
    'NAG': {
        'logp': -3.078,
        'heavy_atom_count': 15,
        'numH_acceptors': 6,
        'numH_donors': 5,
        'num_rotable_bonds': 7,
        'rings_count': 1,
        'TPSA': 119.250,
        'molwt': 221.09
    }
}


class TestPropertyCalculation:

    @staticmethod
    @pytest.mark.parametrize('key', test_inputs)
    def test_valid_properties(key):
        physchem_props = ccd_reader.read_pdb_cif_file(cif_filename(key)).component.physchem_properties

        assert test_inputs[key]['logp'] == round(physchem_props['CrippenClogP'], 3)
        assert test_inputs[key]['heavy_atom_count'] == physchem_props['NumHeavyAtoms']
        assert test_inputs[key]['numH_acceptors'] == physchem_props['NumHBA']
        assert test_inputs[key]['numH_donors'] == physchem_props['NumHBD']
        assert test_inputs[key]['num_rotable_bonds'] == physchem_props['NumRotatableBonds']
        assert test_inputs[key]['rings_count'] == physchem_props['NumRings']
        assert test_inputs[key]['TPSA'] == round(physchem_props['tpsa'], 3)
        assert test_inputs[key]['molwt'] == round(physchem_props['exactmw'], 3)

    @staticmethod
    @pytest.mark.parametrize('key', ['10R', '08T'])
    def test_invalid_properties(key):
        physchem_props = ccd_reader.read_pdb_cif_file(cif_filename(key)).component.physchem_properties

        assert physchem_props == {}
