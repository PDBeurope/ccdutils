import os
import re
import json
import pytest
from gemmi import cif

import xml.etree.ElementTree as ET
from pdbeccdutils.helpers.drawing import svg_namespace
from pdbeccdutils.tests.tst_utilities import updated_mmcif_filename
from pdbeccdutils.scripts.boundmolecule_cli import (
    PDBeBmManager,
    create_parser,
)


class TestBoundmoleculeProcessig:
    """
    run boundmolecule_cli on test file:
    """

    CLC_IDS = ["CLC_1", "CLC_2", "CLC_3", "CLC_4", "CLC_5"]

    @pytest.fixture(scope="class")
    def pipeline_wd(self, tmpdir_factory):
        wd = tmpdir_factory.mktemp("pdbechem_test_bm")
        parser = create_parser()
        args = parser.parse_args(
            ["-i", updated_mmcif_filename("1c4q"), "-o", str(wd), "-id", "1c4q"]
        )
        bm_manager = PDBeBmManager()
        bm_manager.process_entry(args.input_cif, args.pdb_id, args.output_dir)
        return args.output_dir

    @staticmethod
    @pytest.mark.parametrize("clc_id", CLC_IDS)
    def test_subdir_tree_created(pipeline_wd, clc_id):
        path = os.path.join(pipeline_wd, clc_id)
        assert os.path.isdir(path)

    @staticmethod
    def test_boundmolecule_json(pipeline_wd: str):
        file = os.path.join(pipeline_wd, "bound_molecules.json")
        assert os.path.isfile(file)
        assert os.path.getsize(file) > 0

    @staticmethod
    def test_processed_cif(pipeline_wd: str):
        file = os.path.join(pipeline_wd, "1c4q_processed.cif")
        assert os.path.isfile(file)
        assert os.path.getsize(file) > 0

    @staticmethod
    @pytest.mark.parametrize("clc_id", CLC_IDS)
    def test_all_bm_files_created(pipeline_wd: str, clc_id):
        path = os.path.join(pipeline_wd, clc_id)
        files = os.listdir(path)
        assert len(files) == 15

        for f in files:
            file_pointer = os.path.join(pipeline_wd, clc_id, f)
            assert os.path.getsize(file_pointer) > 0

    @staticmethod
    @pytest.mark.parametrize("clc_id", CLC_IDS)
    def test_images_with_names_created(pipeline_wd: str, clc_id):
        """Test if the depictions with names contain certain atom labels
        as expected.
        """
        path = os.path.join(pipeline_wd, clc_id, f"{clc_id}_100_names.svg")
        xml = ET.parse(path)
        path_elements = xml.findall("svg:path", svg_namespace)

        for e in path_elements:
            if "class" in e.attrib:
                if re.fullmatch(r"atom-\d+", e.attrib["class"]):
                    return  # label is where it should be we can quit

        assert False  # this compound should have had label

    @staticmethod
    @pytest.mark.parametrize("clc_id", CLC_IDS)
    def test_cif_files(pipeline_wd: str, clc_id):
        """Test if the CIF file is parsable."""
        path = os.path.join(pipeline_wd, clc_id, f"{clc_id}.cif")
        cif_block = cif.read(path).sole_block()

        assert cif_block
        assert cif_block.name == clc_id

    @staticmethod
    @pytest.mark.parametrize("clc_id", CLC_IDS)
    def test_cml_files(pipeline_wd: str, clc_id):
        """Test if the CML file is parsable. Each of the tested compounds
        should contain element with type C.
        """
        path = os.path.join(pipeline_wd, clc_id, f"{clc_id}.cml")

        xml_root = ET.parse(str(path)).getroot()
        atoms = xml_root.find("molecule").find("atomArray").findall("atom")

        assert any(a.attrib["elementType"] == "C" for a in atoms)

    @staticmethod
    @pytest.mark.parametrize("clc_id", CLC_IDS)
    def test_annotation_file(pipeline_wd: str, clc_id):
        """Test if the annotation.json is parsable. Each of the tested compounds
        should contain some data.
        """

        file = os.path.join(pipeline_wd, clc_id, f"{clc_id}_annotation.json")

        assert os.path.isfile(file)

        with open(file, "r") as fp:
            data = json.load(fp)
            assert data
            assert data["atoms"]
            assert data["bonds"]
            assert data["resolution"]
