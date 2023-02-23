"""Test fragment library functionality
"""

import os


def test__img_crated(library, tmpdir):
    file_path = str(tmpdir.join("library.svg"))

    library.to_image(file_path)

    assert os.path.isfile(file_path)


def test_img_png_created(library, tmpdir):
    file_path = str(tmpdir.join("library.png"))
    library.to_image(file_path, source="PDBe")

    assert os.path.isfile(file_path)


def test_generate_conformers(library):
    library.generate_conformers()

    for entry in library.library.values():
        assert entry.mol.GetConformers()
