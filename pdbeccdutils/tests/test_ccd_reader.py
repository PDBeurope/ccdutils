import logging
from pathlib import Path

from gemmi import cif
from rdkit import Chem

from pdbeccdutils.core import ccd_reader
from pdbeccdutils.tests.tst_utilities import cif_filename


def _components_cif(tmp_path):
    path = tmp_path / "components.cif"
    ccds = [Path(cif_filename("00O")), Path(cif_filename("007"))]
    path.write_text("\n".join(ccd.read_text() for ccd in ccds))
    return path


def test_read_pdb_components_file_parses_all_for_empty_include(tmp_path):
    path = _components_cif(tmp_path)

    for include in (None, []):
        result = ccd_reader.read_pdb_components_file(
            str(path), sanitize=False, include=include
        )

        assert list(result) == ["00O", "007"]


