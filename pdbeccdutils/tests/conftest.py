"""Common fixtures shared among all the tests
"""
import pytest
from pdbeccdutils.core import ccd_reader
from pdbeccdutils.tests.tst_utilities import supply_list_of_sample_cifs

sample_ccd_cifs = supply_list_of_sample_cifs()
problematic_ids = ['UNL', 'NA', 'SY9', '10R', 'ASX']


@pytest.fixture(scope='session', params=sample_ccd_cifs)
def component(request):
    reader = ccd_reader.read_pdb_cif_file(request.param)
    c = reader.component

    if c.id not in problematic_ids:
        assert reader.warnings == []

    return c
