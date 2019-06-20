import os

import pytest

from pdbeccdutils.core import ccd_reader, ccd_writer
from pdbeccdutils.core.models import ConformerType
from pdbeccdutils.tests.tst_utilities import supply_list_of_sample_cifs

sample_ccd_cifs = supply_list_of_sample_cifs()


class TestSDF:
    @staticmethod
    @pytest.mark.parametrize('test_ccd_cif', sample_ccd_cifs)
    def test_write_sdf(test_ccd_cif, tmpdir_factory):
        wd = tmpdir_factory.mktemp('sdf_test')
        reader = ccd_reader.read_pdb_cif_file(test_ccd_cif)
        component = reader.component

        assert reader.errors == []
        for ideal in True, False:
            for remove_hs in True, False:
                sdf_file = os.path.join(wd, component.id)
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

                ccd_writer.write_molecule(path=sdf_file,
                                          component=component,
                                          conf_type=conf_type,
                                          remove_hs=remove_hs)
                assert os.path.isfile(sdf_file)

                # UNL component does not have coordinates. Only the default ideal coords are created.
                if component.id != 'UNL':
                    assert os.path.getsize(sdf_file) > 0
