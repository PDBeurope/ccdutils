"""
Unittest of the process_components_cif_cli.py command line script

method based on
http://dustinrcollins.com/testing-python-command-line-apps
adapted to use nose then converted to pytest
"""
import glob
import os
import pytest
import shutil

from pdbeccdutils.scripts.process_components_cif_cli import create_parser, process_components_cif
from pdbeccdutils.tests.tst_utilities import test_cut_down_components_cif, \
    file_name_in_tsts_out, cif_filename

FILES_SUBDIRS = ('mmcif', 'sdf', 'sdf_nh', 'sdf_r', 'sdf_r_nh',
                 'pdb', 'pdb_r', 'cml', 'xyz', 'xyz_r')
IMAGES_SUBDIRS = 'svg_with_atom_labels', 'svg_without_atom_labels'


@pytest.mark.skip(reason='Not yet implemented in this branch')  # TODO implement osmart 14 May 2018
class TestCommandLineArgs(object):
    @staticmethod
    def test_with_empty_args():
        """
        User passes no args, should produce a usage statement and then
        raise SystemExit. Usage statement will appear
        """
        parser = create_parser()
        with pytest.raises(SystemExit):
            parser.parse_args([])

    @staticmethod
    def test_input_file_that_cannot_exist_raises_system_exit():
        parser = create_parser()
        args = parser.parse_args(['/////impossible_to_open_file', '--debug'])
        with pytest.raises(SystemExit):
            process_components_cif(args)


@pytest.mark.skip(reason='Not yet implemented in this branch')  # TODO implement osmart 14 May 2018
class TestRegressionTest(object):
    @staticmethod
    def test_with_problematic_cif_7om():
        """ 7OM caused problems because it lacks an inchikey
        and this caused a crash """
        cif_file = cif_filename('7OM')
        test_output_dir = file_name_in_tsts_out('test_process_components_cif_70M')
        if os.path.isdir(test_output_dir):
            shutil.rmtree(test_output_dir)
        parser = create_parser()
        args = parser.parse_args([cif_file, '-o', test_output_dir])
        process_components_cif(args)
        assert os.path.isdir(test_output_dir), \
            'output directory  {} must be created'.format(test_output_dir)


@pytest.mark.skip(reason='Not yet implemented in this branch')  # TODO implement osmart 14 May 2018
class TestCutDownComponentsCif(object):
    """
    run process_components_cif_cli on test file:

    cut_down_components.cif

    that is a cutdown components cif with just the first 5 chemical
    components definitions.

    use pytest fixture to do the run once but then have separate tests
    for the creation of each directory file that is required in the
    ftp area.
    """
    CHEM_COMP_IDS = ('000', '001', '002', '003', '004')

    @pytest.fixture(scope='class')
    def output_dir(self):
        parser = create_parser()
        test_components_cif = test_cut_down_components_cif
        output_dir = file_name_in_tsts_out('test_process_components_cif')
        if os.path.isdir(output_dir):
            shutil.rmtree(output_dir)
        args = parser.parse_args([test_components_cif, '-o', output_dir])
        process_components_cif(args)
        return output_dir

    @staticmethod
    def test_output_dir_created(output_dir):
        assert os.path.isdir(output_dir)

    @staticmethod
    def test_subdir_files_created(output_dir):
        files_path = os.path.join(output_dir, 'files')
        assert os.path.isdir(files_path)

    @staticmethod
    @pytest.mark.parametrize('subdir', FILES_SUBDIRS)
    def test_subdirs_in_files_are_created(output_dir, subdir):
        subdir_path = os.path.join(output_dir, 'files', subdir)
        assert os.path.isdir(subdir_path)

    @staticmethod
    @pytest.mark.parametrize('subdir', FILES_SUBDIRS)
    @pytest.mark.parametrize('chem_comp_id', CHEM_COMP_IDS)
    def test_there_is_a_file_for_chem_comp_id(output_dir, subdir, chem_comp_id):
        subdir_path = os.path.join(output_dir, 'files', subdir)
        assert os.path.isdir(subdir_path)
        # simple check that there is a single file starting with the chem_comp_id
        files_for_chem_comp_id = glob.glob1(subdir_path, chem_comp_id + '*')
        assert len(files_for_chem_comp_id) == 1

    @staticmethod
    def test_subdir_images_created(output_dir):
        files_path = os.path.join(output_dir, 'images')
        assert os.path.isdir(files_path)

    @staticmethod
    @pytest.mark.parametrize('subdir', IMAGES_SUBDIRS)
    def test_subdirs_in_images_are_created(output_dir, subdir):
        subdir_path = os.path.join(output_dir, 'images', subdir)
        assert os.path.isdir(subdir_path)

    @staticmethod
    @pytest.mark.parametrize('subdir', IMAGES_SUBDIRS)
    @pytest.mark.parametrize('chem_comp_id', CHEM_COMP_IDS)
    def test_there_is_a_svg_for_chem_comp_id(output_dir, subdir, chem_comp_id):
        subdir_path = os.path.join(output_dir, 'images', subdir, chem_comp_id[:1])
        assert os.path.isdir(subdir_path)
        this_svg = os.path.join(subdir_path, chem_comp_id + '.svg')
        assert os.path.isfile(this_svg) and os.path.getsize(this_svg) > 0

    def test_file_chem_comp_dot_list(self, output_dir):
        chem_comp_dot_list_file = os.path.join(output_dir, 'chem_comp.list')
        try:
            with open(chem_comp_dot_list_file, 'r') as chem_comp_file:
                lines = chem_comp_file.read().splitlines()
                assert lines == list(self.CHEM_COMP_IDS), \
                    'chem_comp.list file should contain list of ccd''s one per line'
        except IOError as message:
            assert False, 'problem opening chem_comp.list file "{}"'.format(message)

    @staticmethod
    def test_file_chem_dot_xml(output_dir):
        chem_dot_xml_file_name = os.path.join(output_dir, 'chem.xml')
        try:
            with open(chem_dot_xml_file_name, 'r') as chem_dot_xml_file:
                lines = chem_dot_xml_file.read().splitlines()
                strip_lines = [item.strip() for item in lines]
                for test_str in ('</chemCompList>', '<id>000</id>', '<id>004</id>',
                                 '<name>(2S)-amino(phenyl)ethanoic acid</name>',
                                 '<fragment id="2" name="phenyl">'):
                    assert test_str in strip_lines, 'chem.xml should contain {}'.format(test_str)
        except IOError as message:
            assert False, 'problem opening chem.xml file "{}"'.format(message)

    @staticmethod
    def test_file_readme_dot_html_file(output_dir):
        this_file = os.path.join(output_dir, 'readme.htm')
        assert os.path.isfile(this_file) and os.path.getsize(this_file) > 0, \
            'readme_dot_html_file {} must be a non-empty file.'.format(this_file)

    @staticmethod
    @pytest.mark.parametrize('subdir', FILES_SUBDIRS)
    def test_tarball_in_files(output_dir, subdir):
        this_file = os.path.join(output_dir, 'files', subdir + '.tar.gz')
        assert os.path.isfile(this_file) and os.path.getsize(this_file) > 0, \
            'tarball {} must be a non-empty file.'.format(this_file)

    @staticmethod
    @pytest.mark.parametrize('subdir', IMAGES_SUBDIRS)
    def test_tarball_in_images(output_dir, subdir):
        this_file = os.path.join(output_dir, 'images', subdir + '.tar.gz')
        assert os.path.isfile(this_file) and os.path.getsize(this_file) > 0, \
            'tarball {} must be a non-empty file.'.format(this_file)
