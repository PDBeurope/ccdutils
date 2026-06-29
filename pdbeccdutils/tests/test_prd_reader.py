import logging

from pdbeccdutils.core import prd_reader
from pdbeccdutils.tests.tst_utilities import prd_cif_filename


def test_read_prd_components_file_parses_all_for_empty_include():
    path = prd_cif_filename("PRDCC_000103")

    for include in (None, []):
        result = prd_reader.read_pdb_components_file(
            path, sanitize=False, include=include
        )

        assert list(result) == ["PRD_000103"]
