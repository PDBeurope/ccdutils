# Unittest of the process_components_cif_cli.py command line script
# method based on http://dustinrcollins.com/testing-python-command-line-apps
# adapted to use nose
import glob
import os
import shutil

from nose.tools import assert_raises, assert_true, assert_equal

from process_components_cif_cli import create_parser, process_components_cif, file_subdirs
from utilities import test_components_cif_first_file_comps, file_name_in_tsts_out


def test_with_empty_args():
    """
    User passes no args, should produce a usage statement and then raise SystemExit. Usage statement will start
    """
    parser = create_parser()
    assert_raises(SystemExit, parser.parse_args, [])


def test_with_components_cif_first_file_comps():
    parser = create_parser()
    test_components_cif = test_components_cif_first_file_comps
    test_output_dir = file_name_in_tsts_out('test_process_components_cif_cli')
    if os.path.isdir(test_output_dir):
        shutil.rmtree(test_output_dir)
    chem_comp_ids = ('000', '001', '002', '003', '004')
    args = parser.parse_args([test_components_cif, test_output_dir, '--debug'])
    process_components_cif(args.COMPONENTS_CIF, args.OUTPUT_DIR, args.debug)
    yield assert_true, os.path.isdir(test_output_dir), 'output directory {} must be created'.format(test_output_dir)
    files_dir = os.path.join(test_output_dir, 'files')
    yield assert_true, os.path.isdir(test_output_dir), 'files sub-directory {} must be created'.format(files_dir)
    for subdir in file_subdirs:
        path = os.path.join(files_dir, subdir)
        yield assert_true, os.path.isdir(path), '{} sub-directory {} must be created'.format(subdir, path)
        for chem_comp_id in chem_comp_ids:
            # simple check that there is a single file starting with the chem_comp_id
            file_for_chem_comp_id = glob.glob1(path, chem_comp_id + '*')
            yield assert_equal, len(file_for_chem_comp_id), 1, \
                'there should be a file matching {}* in {}'.format(chem_comp_id, subdir)
    chem_comp_dot_list_file = os.path.join(test_output_dir, 'chem_comp.list')
    yield assert_true, os.path.isfile(chem_comp_dot_list_file), \
        'chem_comp.list file: {} must be created'.format(chem_comp_dot_list_file)
    with open(chem_comp_dot_list_file, 'r') as chem_comp_file:
        lines = chem_comp_file.read().splitlines()
        yield assert_equal, lines, list(chem_comp_ids), \
            'chem_comp.list file should contain list of ccd''s one per line'
