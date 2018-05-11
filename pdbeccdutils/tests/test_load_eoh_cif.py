"""
load EOH.cif from file and test important cif item.
"""
from pdbeccdutils.tests.tst_utilities import cif_filename
from pdbeccdutils.core import structure_reader as sr

class TestLoadEOH:
    @staticmethod
    def test_load_eoh():
        """ test eoh?? """
        cif_file = cif_filename('EOH')
        reader = sr.read_pdb_cif_file(cif_file)
        assert reader.warnings == []
        assert reader.errors == []
        component = reader.component
        assert component.id == 'EOH'
        assert component.name == 'ETHANOL'
        assert component.formula == 'C2 H6 O'
        # really want to test _chem_comp.pdbx_release_status
        assert component.released == True
