"""Test fragment library functionality
"""

import os
import xml.etree.ElementTree as ET
from pdbeccdutils.helpers.drawing import svg_namespace


def test__img_crated(library, tmpdir):
    file_path = str(tmpdir.join("library.svg"))

    library.to_image(file_path)

    assert os.path.isfile(file_path)

    xml = ET.parse(file_path)
    assert len(xml.findall("svg:rect", svg_namespace)) > 100


def test_img_png_created(library, tmpdir):
    file_path = str(tmpdir.join("library.png"))
    library.to_image(file_path, source="PDBe")

    assert os.path.isfile(file_path)


def test_generate_conformers(library):
    library.generate_conformers()

    for entry in library.library.values():
        assert entry.mol.GetConformers()
