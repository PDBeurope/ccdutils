"""Common fixtures shared among all the tests
"""
import pytest
from pdbeccdutils.core import ccd_reader, clc_reader
from pdbeccdutils.core.fragment_library import FragmentLibrary
from pdbeccdutils.tests.tst_utilities import supply_list_of_sample_cifs
from pdbeccdutils.tests.tst_utilities import updated_mmcif_filename

sample_ccd_cifs = supply_list_of_sample_cifs()
problematic_ids = ["UNL", "NA", "SY9", "10R", "ASX", "0KA"]


@pytest.fixture(scope="session", params=sample_ccd_cifs)
def component(request):
    reader = ccd_reader.read_pdb_cif_file(request.param)
    c = reader.component

    if c.id not in problematic_ids:
        assert reader.warnings == []

    return c


@pytest.fixture(scope="session")
def library():
    return FragmentLibrary()


@pytest.fixture(scope="session")
def component_globotriose():
    """
    load 1c4q_processed.cif.gz and returns the
    component object of globotriose

    Returns:
        pdbeccdutils.core.component.Component: component object for
        globotriose
    """
    cif_file = updated_mmcif_filename("1c4q")
    reader_result = clc_reader.read_pdb_cif_file(cif_file)[0]
    assert reader_result.warnings == []
    assert reader_result.errors == []
    component = reader_result.component

    return component
