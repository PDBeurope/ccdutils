"""Unit tests for checking whether or not depictions are written properly.
"""


import os

import pytest

from pdbeccdutils.core import ccd_reader
from pdbeccdutils.tests.tst_utilities import cif_filename
from pdbeccdutils.core.depictions import DepictionManager


def load_molecule(id):
    depiction = DepictionManager()
    c = ccd_reader.read_pdb_cif_file(cif_filename(id)).component
    c.sanitize()
    c.compute_2d(depiction)
    return c


class TestWriteImg:

    @staticmethod
    def test_file_generated(tmpdir):  # tmpdir is a fixture with temporary directory
        mol = load_molecule('ATP')
        path = str(tmpdir.join('atp_test.svg'))
        mol.export_2d_svg(path)

        assert os.path.isfile(path)

    @staticmethod
    @pytest.mark.parametrize("id,expected,names", [
        ("NAG", 'C8', True),
        ("ATP", 'C5&apos;', True),
        ("08T", 'BE', True),
        ("BCD", 'C66', True),
        ("ATP", '<svg:rect', False),
        ("08T", '<svg:rect', False),
        ("10R", 'Image not', False),
        ("0OD", '<svg:rect', False),
    ])
    def test_image_generation_with_names(tmpdir, id, expected, names):
        mol = load_molecule(id)
        path = str(tmpdir.join('{}_{}.svg'.format(id, 'names' if names else 'no_names')))
        mol.export_2d_svg(path, names=names)

        with open(path, 'r') as f:
            content = f.readlines()
        assert any(expected in i for i in content)
