"""Common fixtures shared among all the tests
"""
import pytest
from pdbeccdutils.core import ccd_reader, bm_reader
from pdbeccdutils.core.fragment_library import FragmentLibrary
from pdbeccdutils.tests.tst_utilities import supply_list_of_sample_cifs
from pdbeccdutils.tests.tst_utilities import supply_list_of_cifs_with_boundmolecules

sample_ccd_cifs = supply_list_of_sample_cifs()
sample_cifs_with_boundmolecules = supply_list_of_cifs_with_boundmolecules()
problematic_ids = ["UNL", "NA", "SY9", "10R", "ASX", "0KA"]


@pytest.fixture(scope="session", params=sample_ccd_cifs)
def component(request):
    reader = ccd_reader.read_pdb_cif_file(request.param)
    c = reader.component

    if c.id not in problematic_ids:
        assert reader.warnings == []

    return c


@pytest.fixture(scope="session", params=sample_cifs_with_boundmolecules)
def boundMolecules(request):
    reader_list = bm_reader.read_pdb_updated_cif_file(request.param)
    bm_components = []
    for reader_result in reader_list:
        c = reader_result.component
        assert reader_result.warnings == []
        bm_components.append(c)

    return bm_components


@pytest.fixture(scope="session")
def library():
    return FragmentLibrary()
