"""
load 6lq4_updated.cif from file and test important cif item.
"""
import pytest
from pdbeccdutils.core import clc_reader
from pdbeccdutils.core.models import ReleaseStatus
from pdbeccdutils.tests.tst_utilities import updated_mmcif_filename


def test_clc_reader_wo_clc():
    path_to_cif = updated_mmcif_filename("1tqh")
    clc_reader_result = clc_reader.read_pdb_cif_file(path_to_cif)
    assert clc_reader_result == []


def test_clc_reader_with_clc():
    path_to_cif = updated_mmcif_filename("1c4q")
    clc_reader_result = clc_reader.read_pdb_cif_file(path_to_cif)
    assert len(clc_reader_result) > 0


clc_items = [
    ("formula", "C18H32O16"),
    ("pdbx_release_status", ReleaseStatus.REL),
    ("inchikey", "FYGDTMLNYKFZSV-SKWQFERISA-N"),
    (
        "inchi",
        "InChI=1S/C18H32O16/c19-1-4-7(22)8(23)12(27)17(31-4)34-15-6(3-21)32-18(13(28)10(15)25)33-14-5(2-20)30-16(29)11(26)9(14)24/h4-29H,1-3H2/t4-,5-,6-,7+,8+,9-,10-,11-,12-,13-,14-,15+,16-,17-,18+/m1/s1",
    ),
    ("number_atoms", 66),
]


@pytest.mark.parametrize("attribute, expected", clc_items)
def clc_items(attribute, expected, component_globotriose):
    assert hasattr(component_globotriose, attribute)
    assert getattr(component_globotriose, attribute) == expected
