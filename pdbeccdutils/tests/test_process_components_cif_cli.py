"""
Unittest of the process_components_cif_cli.py command line script

method based on
http://dustinrcollins.com/testing-python-command-line-apps
adapted to use nose then converted to pytest
"""
import json
import os
import re
import xml.etree.ElementTree as ET
from gemmi import cif

import pytest
from pdbeccdutils.helpers.drawing import svg_namespace
from pdbeccdutils.helpers.cif_tools import get_prd_code
from pdbeccdutils.scripts.process_components_cif_cli import (
    PDBeChemManager,
    create_parser,
)
from pdbeccdutils.tests.tst_utilities import cif_filename, prd_cif_filename

CHEM_COMP_IDS = ["00O"]
CHEM_COMP_NAMES = {
    "00O": {"name": "HB1", "alt_name": "H8"},
    "007": {"name": "H1C1", "alt_name": "1H1C"},
}
PRDCC_IDS = ["PRDCC_000103"]


class TestCommandLineArgs:
    @staticmethod
    def test_with_empty_args():
        """
        User passes no args, should produce a usage statement and then
        raise SystemExit. Usage statement will appear
        """
        parser = create_parser()
        with pytest.raises(SystemExit):
            parser.parse_args()

    @staticmethod
    def test_input_file_that_cannot_exist_raises_system_exit():
        parser = create_parser()
        with pytest.raises(SystemExit):
            parser.parse_args(["-o foo", "/////impossible_to_open_file", "--debug"])


@pytest.fixture(scope="session", params=CHEM_COMP_IDS)
def ccd_prd_pipeline_data(tmpdir_factory, request):
    wd = tmpdir_factory.mktemp("pdbechem_ccd_test")
    parser = create_parser()
    chem_comp_id = request.param
    args = parser.parse_args(["-o", str(wd), "-i", cif_filename(chem_comp_id)])

    m = PDBeChemManager()
    m.run(args.input_cif, args.output_dir)
    return args.output_dir, chem_comp_id


@pytest.fixture(scope="session", params=PRDCC_IDS)
def prd_pipeline_data(tmpdir_factory, request):
    wd = tmpdir_factory.mktemp("pdbechem_prd_test")
    prdcc_id = request.param
    parser = create_parser()
    args = parser.parse_args(["-o", str(wd), "-i", prd_cif_filename(prdcc_id), "--prd"])

    m = PDBeChemManager(procedure=args.procedure)
    m.run(args.input_cif, args.output_dir)

    return args.output_dir, prdcc_id


class TestProcessComponentsCif:
    """
    run process_components_cif_cli on CCD test file

    use pytest fixture to do the run once but then have separate tests
    for the creation of each directory file that is required in the
    ftp area.
    """

    @staticmethod
    def test_subdir_tree_created(ccd_prd_pipeline_data):
        path = (
            ccd_prd_pipeline_data[0]
            / ccd_prd_pipeline_data[1][0]
            / ccd_prd_pipeline_data[1]
        )
        assert path.is_dir()

    @staticmethod
    def test_all_files_created(ccd_prd_pipeline_data):
        path = (
            ccd_prd_pipeline_data[0]
            / ccd_prd_pipeline_data[1][0]
            / ccd_prd_pipeline_data[1]
        )
        files = os.listdir(path)
        assert len(files) == 19

        for f in files:
            file_pointer = path / f
            assert file_pointer.stat().st_size > 0

    @staticmethod
    def test_images_with_names_created(ccd_prd_pipeline_data):
        """Test if the depictions with names contain certain atom labels
        as expected.
        """
        path = (
            ccd_prd_pipeline_data[0]
            / ccd_prd_pipeline_data[1][0]
            / ccd_prd_pipeline_data[1]
            / f"{ccd_prd_pipeline_data[1]}_100_names.svg"
        )
        xml = ET.parse(path)
        path_elements = xml.findall("svg:path", svg_namespace)

        for e in path_elements:
            if "class" in e.attrib:
                if re.fullmatch(r"atom-\d+", e.attrib["class"]):
                    return  # label is where it should be we can quit

        assert False  # this compound should have had label

    @staticmethod
    def test_correct_atom_naming_in_files(ccd_prd_pipeline_data):
        """Test if alternate names are used for exported model/ideal
        pdb files.
        """
        alts = [
            ccd_prd_pipeline_data[0]
            / ccd_prd_pipeline_data[1][0]
            / ccd_prd_pipeline_data[1]
            / f"{ccd_prd_pipeline_data[1]}_ideal_alt.pdb",
            ccd_prd_pipeline_data[0]
            / ccd_prd_pipeline_data[1][0]
            / ccd_prd_pipeline_data[1]
            / f"{ccd_prd_pipeline_data[1]}_model_alt.pdb",
        ]
        regular = [
            ccd_prd_pipeline_data[0]
            / ccd_prd_pipeline_data[1][0]
            / ccd_prd_pipeline_data[1]
            / f"{ccd_prd_pipeline_data[1]}_ideal.pdb",
            ccd_prd_pipeline_data[0]
            / ccd_prd_pipeline_data[1][0]
            / ccd_prd_pipeline_data[1]
            / f"{ccd_prd_pipeline_data[1]}_model.pdb",
        ]

        name = CHEM_COMP_NAMES[ccd_prd_pipeline_data[1]]["name"]
        alt_name = CHEM_COMP_NAMES[ccd_prd_pipeline_data[1]]["alt_name"]

        for i in alts:
            with open(i) as f:
                str_repr = f.read()
                assert alt_name in str_repr and name not in str_repr

        for i in regular:
            with open(i) as f:
                str_repr = f.read()
                assert name in str_repr and alt_name not in str_repr

    @staticmethod
    def test_cif_files(ccd_prd_pipeline_data):
        """Test if the CIF file is parsable."""
        path = (
            ccd_prd_pipeline_data[0]
            / ccd_prd_pipeline_data[1][0]
            / ccd_prd_pipeline_data[1]
            / f"{ccd_prd_pipeline_data[1]}.cif"
        )
        cif_block = cif.read(str(path)).sole_block()
        assert cif_block.name == ccd_prd_pipeline_data[1]

    @staticmethod
    def test_cml_files(ccd_prd_pipeline_data):
        """Test if the CML file is parsable. Each of the tested compounds
        should contain element with type C.
        """
        path = (
            ccd_prd_pipeline_data[0]
            / ccd_prd_pipeline_data[1][0]
            / ccd_prd_pipeline_data[1]
            / f"{ccd_prd_pipeline_data[1]}.cml"
        )

        xml_root = ET.parse(str(path)).getroot()
        atoms = xml_root.find("molecule").find("atomArray").findall("atom")

        assert any(a.attrib["elementType"] == "C" for a in atoms)

    @staticmethod
    def test_annotation_file(ccd_prd_pipeline_data):
        """Test if the annotation.json is parsable. Each of the tested compounds
        should contain some data.
        """

        path = (
            ccd_prd_pipeline_data[0]
            / ccd_prd_pipeline_data[1][0]
            / ccd_prd_pipeline_data[1]
            / f"{ccd_prd_pipeline_data[1]}_annotation.json"
        )

        assert path.is_file()

        with open(path, "r") as fp:
            data = json.load(fp)

            assert data
            assert data["atoms"]
            assert data["bonds"]
            assert data["resolution"]


class TestProcessPRDCif:
    """
    run process_components_cif_cli on PRD test file
    """

    @staticmethod
    def test_subdir_tree_created(prd_pipeline_data):
        path = prd_pipeline_data[0] / prd_pipeline_data[1][-1] / prd_pipeline_data[1]
        assert path.is_dir()

    @staticmethod
    def test_all_files_created(prd_pipeline_data):
        path = prd_pipeline_data[0] / prd_pipeline_data[1][-1] / prd_pipeline_data[1]
        files = os.listdir(path)
        assert len(files) == 19

        for f in files:
            file_pointer = path / f
            assert file_pointer.stat().st_size > 0

    @staticmethod
    def test_cif_files(prd_pipeline_data):
        """Test if the CIF file is parsable."""
        path = (
            prd_pipeline_data[0]
            / prd_pipeline_data[1][-1]
            / prd_pipeline_data[1]
            / f"{prd_pipeline_data[1]}.cif"
        )
        cif_block = cif.read(str(path)).sole_block()
        assert cif_block.name == get_prd_code(prd_pipeline_data[1])

    @staticmethod
    def test_cif_enriched_field_comp_id(prd_pipeline_data):
        """Test if the comp_id of enriched fields in CIF is same as comp_id of file."""
        path = (
            prd_pipeline_data[0]
            / prd_pipeline_data[1][-1]
            / prd_pipeline_data[1]
            / f"{prd_pipeline_data[1]}.cif"
        )
        cif_block = cif.read(str(path)).sole_block()
        enriched_field = "_pdbe_chem_comp_atom_depiction."
        if enriched_field in cif_block.get_mmcif_category_names():
            enriched_field_block = cif_block.find([enriched_field, "comp_id"])
            for row in enriched_field_block:
                assert row["comp_id"] == cif_block.name
