# Unittest of the process_components_cif_cli.py command line script
# method based on http://dustinrcollins.com/testing-python-command-line-apps
# adapted to use nose
import glob
import os
import shutil

from nose.tools import assert_raises, assert_true, assert_equal

from pdbeccdutils.process_components_cif_cli import create_parser, process_components_cif, file_subdirs, images_subdirs
from pdbeccdutils.utilities import test_components_cif_first_file_comps, file_name_in_tsts_out, cif_filename


def test_with_empty_args():
    """
    User passes no args, should produce a usage statement and then raise SystemExit. Usage statement will appear
    """
    parser = create_parser()
    parser.parse_args
    assert_raises(SystemExit, parser.parse_args, [])


def test_input_file_that_cannot_exist_raises_system_exit():
    parser = create_parser()
    args = parser.parse_args(['/////impossible_to_open_file', '--debug'])
    assert_raises(SystemExit, process_components_cif, args)


def test_with_problematic_cif_7om():
    cif_file = cif_filename('7OM')
    test_output_dir = file_name_in_tsts_out('test_process_components_cif_70M')
    if os.path.isdir(test_output_dir):
        shutil.rmtree(test_output_dir)
    parser = create_parser()
    args = parser.parse_args([cif_file, '-o', test_output_dir, '--debug'])
    process_components_cif(args)
    yield assert_true, os.path.isdir(test_output_dir), 'output directory  {} must be created'.format(test_output_dir)


def test_with_components_cif_first_file_comps():
    parser = create_parser()
    test_components_cif = test_components_cif_first_file_comps
    test_output_dir = file_name_in_tsts_out('test_process_components_cif')
    if os.path.isdir(test_output_dir):
        shutil.rmtree(test_output_dir)
    chem_comp_ids = ('000', '001', '002', '003', '004')
    args = parser.parse_args([test_components_cif, '-o', test_output_dir, '--debug'])
    process_components_cif(args)
    yield assert_true, os.path.isdir(test_output_dir), 'output directory {} must be created'.format(test_output_dir)
    files_dir = os.path.join(test_output_dir, 'files')
    yield assert_true, os.path.isdir(files_dir), 'files sub-directory {} must be created'.format(files_dir)
    for subdir in file_subdirs:
        path = os.path.join(files_dir, subdir)
        yield assert_true, os.path.isdir(path), '{} sub-directory {} must be created'.format(subdir, path)
        for chem_comp_id in chem_comp_ids:
            # simple check that there is a single file starting with the chem_comp_id
            file_for_chem_comp_id = glob.glob1(path, chem_comp_id + '*')
            yield assert_equal, len(file_for_chem_comp_id), 1, \
                'there should be a file matching {}* in {}'.format(chem_comp_id, subdir)
    chem_comp_dot_list_file = os.path.join(test_output_dir, 'chem_comp.list')
    try:
        with open(chem_comp_dot_list_file, 'r') as chem_comp_file:
            lines = chem_comp_file.read().splitlines()
            yield assert_equal, lines, list(chem_comp_ids), \
                'chem_comp.list file should contain list of ccd''s one per line'
    except IOError as message:
        yield assert_true, False, 'problem opening chem_comp.list "{}"'.format(message)
    chem_dot_xml_file_name = os.path.join(test_output_dir, 'chem.xml')
    try:
        with open(chem_dot_xml_file_name, 'r') as chem_dot_xml_file:
            lines = chem_dot_xml_file.read().splitlines()
            strip_lines = [item.strip() for item in lines]
            for test_str in ('</chemCompList>', '<id>000</id>', '<id>004</id>',
                             '<name>(2S)-amino(phenyl)ethanoic acid</name>', '<fragment id="2" name="phenyl">'):
                yield assert_true, test_str in strip_lines, 'chem.xml should contain {}'.format(test_str)
            # yield assert_true, '<id>000</id>' in strip_lines, 'chem.xml should contain <id>000</id>'
            # yield assert_true, '<id>004</id>'in strip_lines, 'chem.xml should contain <id>004</id>'
            # yield assert_true, '<name>(2S)-amino(phenyl)ethanoic acid</name>' in strip_lines, \
            #     'chem.xml should contain name of 004'
            # yield assert_true, '<fragment id="2" name="phenyl"> 'in strip_lines, 'chem.xml should contain '

    except IOError as message:
        yield assert_true, False, 'problem opening chem.xml "{}"'.format(message)
    images_dir = os.path.join(test_output_dir, 'images')
    yield assert_true, os.path.isdir(images_dir), 'images sub-directory {} must be created'.format(images_dir)
    for subdir in images_subdirs:
        yield assert_true, os.path.isdir(path), '{} sub-directory {} must be created'.format(subdir, path)
        for chem_comp_id in chem_comp_ids:
            if subdir[:3] == 'svg':
                file_type = 'svg'
                path = os.path.join(images_dir, subdir, chem_comp_id[:1])
            else:
                path = os.path.join(images_dir, subdir)
                file_type = 'gif'
            # simple check that there is a single file starting with the chem_comp_id
            file_for_chem_comp_id = glob.glob1(path, chem_comp_id + '*.' + file_type)
            yield assert_equal, len(file_for_chem_comp_id), 1, \
                'there should be an {} file matching {}*.{} in {}'.\
                format(file_type, file_type, chem_comp_id, path)
    readme_dot_html_file = os.path.join(test_output_dir, 'readme.htm')
    yield assert_true, os.path.isfile(readme_dot_html_file) and os.path.getsize(readme_dot_html_file) > 0, \
        'readme_dot_html_file {} must be a non-empty file.'.format(readme_dot_html_file)

    for this_dir, this_subdirs in {files_dir: file_subdirs, images_dir: images_subdirs}.items():
        for subdir in this_subdirs:
            tar_ball_file_name = os.path.join(this_dir, subdir + '.tar.gz')
            yield assert_true, os.path.isfile(tar_ball_file_name) and os.path.getsize(tar_ball_file_name) > 0, \
                'tar_ball {} must be a non-empty file.'.format(tar_ball_file_name)