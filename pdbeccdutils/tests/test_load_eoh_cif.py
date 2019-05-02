"""
load EOH.cif from file and test important cif item.
"""
import pytest

from pdbeccdutils.core import ccd_reader, ccd_writer
from pdbeccdutils.core.models import ConformerType, ReleaseStatus
from pdbeccdutils.tests.tst_utilities import cif_filename


class TestLoadEOH:

    @staticmethod
    @pytest.fixture(scope='class')
    def component_eoh():
        """
        loads PDB-CCD for ethanol eoh.cif and returns the
        component object

        Returns:
            pdbeccdutils.core.component.Component: component object for
            eoh.cif
        """
        cif_file = cif_filename('EOH')
        reader = ccd_reader.read_pdb_cif_file(cif_file)
        assert reader.warnings == []
        assert reader.errors == []
        component = reader.component
        return component

    @staticmethod
    def test_chem_comp_dot_id(component_eoh):
        """test retrieval of _chem_comp.id EOH"""
        assert component_eoh.id == 'EOH'

    chem_comp_items = [('id', 'EOH'),
                       ('name', 'ETHANOL'),
                       ('formula', 'C2 H6 O'),
                       ('pdbx_release_status', ReleaseStatus.REL),
                       ]

    @staticmethod
    @pytest.mark.parametrize('attribute, expected', chem_comp_items)
    def test_chem_comp_items(attribute, expected, component_eoh):
        """ for EOH.cif tests items from _chem_comp:
            _chem_comp.id                                    EOH
            _chem_comp.name                                  ETHANOL
            _chem_comp.formula                               "C2 H6 O"
            _chem_comp.pdbx_release_status                   REL
        """
        assert hasattr(component_eoh, attribute)
        assert getattr(component_eoh, attribute) == expected

    @staticmethod
    def test_released_is_true(component_eoh):
        assert component_eoh.released is True

    @staticmethod
    def test_inchikey(component_eoh):
        """
        test InChIKey retrieval from _pdbx_chem_comp_descriptor line:
        EOH InChIKey InChI 1.03  LFQSCWFLJHTTHZ-UHFFFAOYSA-N
        """
        assert component_eoh.inchikey == 'LFQSCWFLJHTTHZ-UHFFFAOYSA-N'

    @staticmethod
    def test_eoh_has_nine_atoms(component_eoh):
        """ test number of atoms in the _chem_comp_atom table"""
        assert component_eoh.number_atoms == 9

    @staticmethod
    def test_eoh_atom_ids(component_eoh):
        assert component_eoh.atoms_ids == ('C1', 'C2', 'O', 'H11', 'H12',
                                           'H21', 'H22', 'H23', 'HO')

    @staticmethod
    def test_inchikeys_from_rdkit_and_ccd_match(component_eoh):
        assert component_eoh.inchikey == component_eoh.inchikey_from_rdkit

    @staticmethod
    def test_to_sdf_string_ideal_no_h(component_eoh):
        sdf_string = ccd_writer.to_sdf_str(component_eoh)
        # there should be two carbons and no hydrogen atoms:
        assert sdf_string.count(' C ') == 2
        assert sdf_string.count(' O ') == 1
        assert sdf_string.count(' H ') == 0
        # check x and y coordinates of the Oxygen
        assert '1.130' in sdf_string
        assert '0.315' in sdf_string

    @staticmethod
    def test_to_sdf_string_ideal_with_h(component_eoh):
        sdf_string = ccd_writer.to_sdf_str(component_eoh, remove_hs=False)
        # six hydrogen atoms:
        assert sdf_string.count(' H ') == 6
        # check z coordinate of the first hydrogen atom
        assert '0.890' in sdf_string

    @staticmethod
    def test_to_sdf_string_model_no_h(component_eoh):
        sdf_string = ccd_writer.to_sdf_str(component_eoh, conf_type=ConformerType.Model)
        # there should be two carbons and no hydrogen atoms:
        assert sdf_string.count(' C ') == 2
        assert sdf_string.count(' O ') == 1
        assert sdf_string.count(' H ') == 0
        # check x and y coordinates of the Oxygen
        assert '15.861' in sdf_string
        assert '8.256' in sdf_string
