import os
import pytest
from pdbeccdutils.tests.tst_utilities import supply_list_of_sample_cifs, file_name_in_tsts_out
from pdbeccdutils.core import ccd_reader
from pdbeccdutils.core import structure_writer


sample_ccd_cifs = supply_list_of_sample_cifs()

class TestSDF:
    @staticmethod
    @pytest.mark.parametrize('test_ccd_cif', sample_ccd_cifs)
    def test_write_sdf(test_ccd_cif):
        assert os.path.isfile(test_ccd_cif)
        reader = ccd_reader.read_pdb_cif_file(test_ccd_cif)
        assert reader.errors == []
        component = reader.component
        sdf_file = file_name_in_tsts_out(component.id + '.ideal_withH.sdf')
        if os.path.isfile(sdf_file):
            os.remove(sdf_file)
        if component.id in ('10R', 'UNL'): # known problem codes.
            with pytest.raises(ValueError):
                structure_writer.write_molecule(path=sdf_file, component=component,
                                                remove_hs=False)
        else:
            structure_writer.write_molecule(path=sdf_file, component=component,
                                        remove_hs=False)
            assert os.path.isfile(sdf_file)
            assert os.path.getsize(sdf_file) > 0

