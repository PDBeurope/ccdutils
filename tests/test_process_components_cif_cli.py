# Unittest of the process_components_cif_cli.py command line script
# method based on http://dustinrcollins.com/testing-python-command-line-apps

import os
import shutil
from process_components_cif_cli import create_parser, process_components_cif
from unittest import TestCase
from utilities import test_components_cif_first_file_comps, file_name_in_tsts_out


class CommandLineTestCase(TestCase):
    """
    Base TestCase class, sets up a CLI parser
    """
    @classmethod
    def setUpClass(cls):
        parser = create_parser()
        cls.parser = parser


class ProcessComponentsCIFTestCase(CommandLineTestCase):
    def test_with_empty_args(self):
        """
        User passes no args, should fail with SystemExit
        """
        with self.assertRaises(SystemExit):
            self.parser.parse_args([])

    def test_with_components_cif_first_file_comps(self):
        test_components_cif = test_components_cif_first_file_comps
        test_output_dir = file_name_in_tsts_out('test_process_components_cif_cli')
        if os.path.isdir(test_output_dir):
            shutil.rmtree(test_output_dir)
        args = self.parser.parse_args([test_components_cif, test_output_dir, '--debug'])
        process_components_cif(args.COMPONENTS_CIF, args.OUTPUT_DIR,  args.debug)
        self.assertTrue(os.path.isdir(test_output_dir), 'output directory {} must be created'.format(test_output_dir))
        files_dir = os.path.join(test_output_dir, 'files')
        self.assertTrue(os.path.isdir(test_output_dir), 'files sub-directory {} must be created'.format(files_dir))
        mmcif_dir = os.path.join(files_dir, 'mmcif')
        self.assertTrue(os.path.isdir(test_output_dir), 'files sub-directory {} must be created'.format(mmcif_dir))
