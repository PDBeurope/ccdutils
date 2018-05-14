"""
load EOH.cif from file and test important cif item.
"""
import pytest
from pdbeccdutils.tests.tst_utilities import cif_filename
from pdbeccdutils.core import ccd_reader
from pdbeccdutils.core import ReleaseStatus


class TestLoadEOH:
    chem_comp_items = [('id', 'EOH'),
                       ('name', 'ETHANOL' ),
                       ('formula', 'C2 H6 O'),
                       ('pdbx_release_status', ReleaseStatus.REL),
                       ('released', True)
                       ]
    @staticmethod
    @pytest.mark.parametrize('attribute, expected', chem_comp_items)
    def test_chem_comp_items(attribute, expected):
        """ for EOH.cif tests items from _chem_comp:
            _chem_comp.id                                    EOH
            _chem_comp.name                                  ETHANOL
            _chem_comp.formula                               "C2 H6 O"
            _chem_comp.pdbx_release_status                   REL
            and that released is True
        """
        cif_file = cif_filename('EOH')
        reader = ccd_reader.read_pdb_cif_file(cif_file)
        assert reader.warnings == []
        assert reader.errors == []
        component = reader.component
        assert hasattr(component, attribute)
        if type(expected) is type(True):
            assert getattr(component, attribute) is expected
        else:
            assert getattr(component, attribute) == expected
