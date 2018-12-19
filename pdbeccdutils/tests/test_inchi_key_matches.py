import os
import pytest
from pdbeccdutils.tests.tst_utilities import supply_list_of_sample_cifs
from pdbeccdutils.core import ccd_reader

sample_ccd_cifs = supply_list_of_sample_cifs()
sample_ccd_with_inchi_problems = ['7OM', 'ASX', 'CDL', '0OD']


class TestIn:
    @staticmethod
    @pytest.mark.parametrize('test_ccd_cif', sample_ccd_cifs)
    def test_inchikeys_from_rdkit_and_ccd_match(test_ccd_cif):
        assert os.path.isfile(test_ccd_cif)
        reader = ccd_reader.read_pdb_cif_file(test_ccd_cif)
        assert reader.errors == []
        # do not:
        # assert reader.warnings == []
        # as there are warnings from 10R, ASX, NA, SY9 and UNL about
        # "missing" cif Namespaces.
        component = reader.component
        if component.id in sample_ccd_with_inchi_problems:
            pass
        else:
            assert component.inchikey == component.inchikey_from_rdkit
