"""
Unittest of the process_components_cif_cli.py command line script

method based on
http://dustinrcollins.com/testing-python-command-line-apps
adapted to use nose then converted to pytest
"""
import json
import os
import xml.etree.ElementTree as ET

import pytest

from pdbeccdutils.scripts.process_components_cif_cli import (PDBeChemManager,
                                                             check_args,
                                                             create_parser)
from pdbeccdutils.tests.tst_utilities import (cif_filename,
                                              test_cut_down_components_cif)


class TestCommandLineArgs:
    @staticmethod
    def test_with_empty_args():
        """
        User passes no args, should produce a usage statement and then
        raise SystemExit. Usage statement will appear
        """
        parser = create_parser()
        with pytest.raises(SystemExit):
            parser.parse_args()

    @staticmethod
    def test_input_file_that_cannot_exist_raises_system_exit():
        parser = create_parser()
        args = parser.parse_args(['-o foo', '/////impossible_to_open_file', '--debug'])
        with pytest.raises(SystemExit):
            check_args(args)


class TestRegressionTest:
    @staticmethod
    def test_with_problematic_cif_7om(tmpdir):
        """
        7OM caused problems because it lacks an inchikey
        and this used to cause a crash
        """
        cif_file = cif_filename('7OM')
        test_output_dir = tmpdir.mkdir('test_process_components_cif_70M')
        parser = create_parser()
        args = parser.parse_args([cif_file, '-o', str(test_output_dir)])
        m = PDBeChemManager()
        m.run_pipeline(args)

        assert os.path.isdir(test_output_dir), f'output directory {str(test_output_dir)} must be created'


class TestCutDownComponentsCif:
    """
    run process_components_cif_cli on test file:

    cut_down_components.cif

    that is a cutdown components cif with just the first 5 chemical
    components definitions.

    use pytest fixture to do the run once but then have separate tests
    for the creation of each directory file that is required in the
    ftp area.
    """
    CHEM_COMP_IDS = ['000', '001', '002', '003', '004', 'ZPN']

    @pytest.fixture(scope='class')
    def pipeline_wd(self, tmpdir_factory):
        wd = tmpdir_factory.mktemp('pdbechem_test')

        parser = create_parser()
        args = parser.parse_args(['-o', str(wd), test_cut_down_components_cif])

        check_args(args)
        m = PDBeChemManager()
        m.run_pipeline(args)

        return str(wd)

    @staticmethod
    def test_output_dir_created(pipeline_wd):
        assert os.path.isdir(pipeline_wd)

    @staticmethod
    @pytest.mark.parametrize('chem_comp_id', CHEM_COMP_IDS)
    def test_subdir_tree_created(pipeline_wd, chem_comp_id):
        path = os.path.join(pipeline_wd, chem_comp_id[0], chem_comp_id)
        assert os.path.isdir(path)

    @staticmethod
    @pytest.mark.parametrize('chem_comp_id', CHEM_COMP_IDS)
    def test_all_files_created(pipeline_wd, chem_comp_id):
        path = os.path.join(pipeline_wd, chem_comp_id[0], chem_comp_id)
        files = os.listdir(path)
        assert len(files) == 19

        for f in files:
            assert os.path.getsize(os.path.join(path, f)) > 0

    @staticmethod
    @pytest.mark.parametrize('id_,name', [
        ("000", 'OA'),
        ("001", 'F11'),
        ("002", 'N3'),
        ("003", 'C17'),
        ("004", 'OXT')
    ])
    def test_images_with_names_created(pipeline_wd, id_, name):
        """Test if the depictions with names contain certain atom labels
        as expected.
        """
        path = os.path.join(pipeline_wd, id_[0], id_, '{}_100_names.svg'.format(id_))
        pattern = '<tspan>{}</tspan>'.format(name)

        with open(path) as f:
            str_repr = f.read()

            assert pattern in str_repr

    @staticmethod
    @pytest.mark.parametrize('id_,name,alt_name', [
        ('001', 'H021', '1H02'),
        ('002', 'H121', '1H12'),
        ('003', 'H121', '1H12')
    ])
    def test_correct_atom_naming_in_files(pipeline_wd, id_, name, alt_name):
        """Test if alternate names are used for exported model/ideal
        pdb files.
        """
        alts = [os.path.join(pipeline_wd, id_[0], id_, '{}_ideal_alt.pdb'.format(id_)),
                os.path.join(pipeline_wd, id_[0], id_, '{}_model_alt.pdb'.format(id_))]
        regular = [
            os.path.join(pipeline_wd, id_[0], id_, '{}_ideal.pdb'.format(id_)),
            os.path.join(pipeline_wd, id_[0], id_, '{}_model.pdb'.format(id_))
        ]

        for i in alts:
            with open(i) as f:
                str_repr = f.read()
                assert alt_name in str_repr and name not in str_repr

        for i in regular:
            with open(i) as f:
                str_repr = f.read()
                assert name in str_repr and alt_name not in str_repr

    @staticmethod
    @pytest.mark.parametrize('chem_comp_id', CHEM_COMP_IDS)
    def test_cml_files(pipeline_wd, chem_comp_id):
        """Test if the CML file is parsable. Each of the tested compounds
        should contain element with type C.
        """
        path = os.path.join(pipeline_wd, chem_comp_id[0], chem_comp_id, '{}.cml'.format(chem_comp_id))

        xml_root = ET.parse(path).getroot()
        atoms = xml_root.find('molecule').find('atomArray').findall('atom')

        assert any(a.attrib['elementType'] == 'C' for a in atoms)

    @staticmethod
    @pytest.mark.parametrize('chem_comp_id', CHEM_COMP_IDS)
    def test_annotation_file(pipeline_wd, chem_comp_id):
        """Test if the annotation.jsone is parsable. Each of the tested compounds
        should contain some data.
        """

        path = os.path.join(pipeline_wd, chem_comp_id[0], chem_comp_id, f'{chem_comp_id}_annotation.json')
        assert os.path.isfile(path)

        with open(path, 'r') as fp:
            data = json.load(fp)

            assert data
            assert data['atoms']
            assert data['bonds']
            assert data['resolution']