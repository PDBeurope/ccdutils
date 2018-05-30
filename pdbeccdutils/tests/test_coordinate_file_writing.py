import os
import pytest
from pdbeccdutils.tests.tst_utilities import supply_list_of_sample_cifs, file_name_in_tsts_out
from pdbeccdutils.core import ccd_reader
from pdbeccdutils.core import structure_writer
from pdbeccdutils.core import ConformerType


sample_ccd_cifs = supply_list_of_sample_cifs()


class TestSDF:
    @staticmethod
    @pytest.mark.parametrize('test_ccd_cif', sample_ccd_cifs)
    def test_write_sdf(test_ccd_cif):
        assert os.path.isfile(test_ccd_cif)
        reader = ccd_reader.read_pdb_cif_file(test_ccd_cif)
        assert reader.errors == []
        component = reader.component
        for ideal in True, False:
            for remove_hs in True, False:
                sdf_file = file_name_in_tsts_out(component.id)
                if ideal:
                    sdf_file += '_ideal'
                    conf_type = ConformerType.Ideal
                else:
                    sdf_file += '_model'
                    conf_type = ConformerType.Model
                if remove_hs:
                    sdf_file += '_no_h'
                else:
                    sdf_file += '_with_h'
                sdf_file += '.sdf'
                if os.path.isfile(sdf_file):
                    os.remove(sdf_file)
                if component.id in ('10R', 'UNL'):  # known problem codes.
                    with pytest.raises(Exception):
                        structure_writer.write_molecule(path=sdf_file,
                                                        component=component,
                                                        conf_type=conf_type,
                                                        remove_hs=remove_hs)
                else:
                    structure_writer.write_molecule(path=sdf_file,
                                                    component=component,
                                                    conf_type=conf_type,
                                                    remove_hs=remove_hs)
                    assert os.path.isfile(sdf_file)
                    assert os.path.getsize(sdf_file) > 0
