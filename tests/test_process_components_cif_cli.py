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
from pathlib import Path

import pytest
from pdbeccdutils.helpers.drawing import svg_namespace
from pdbeccdutils.scripts.process_components_cif_cli import (
    PDBeChemManager,
    create_parser,
)
from pdbeccdutils.tests.tst_utilities import test_cut_down_components_cif


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


class TestCutDownComponentsCif:
    """
    run process_components_cif_cli on test file:

    cut_down_components.cif

    that is a cutdown components cif with just the first 5 chemical
    components definitions.

    use pytest fixture to do the run once but then have separate tests
    for the creation of each directory file that is required in the
    ftp area.
    """

    CHEM_COMP_IDS = ["000", "001", "002", "003", "004", "ZPN"]

    @pytest.fixture(scope="class")
    def pipeline_wd(self, tmpdir_factory):
        wd = tmpdir_factory.mktemp("pdbechem_test")

        parser = create_parser()
        args = parser.parse_args(["-o", str(wd), test_cut_down_components_cif])

        m = PDBeChemManager()
        m.run(args.components_cif, args.output_dir)

        return args.output_dir

    @staticmethod
    @pytest.mark.parametrize("chem_comp_id", CHEM_COMP_IDS)
    def test_subdir_tree_created(pipeline_wd, chem_comp_id):
        path = pipeline_wd / chem_comp_id[0] / chem_comp_id
        assert path.is_dir()

    @staticmethod
    @pytest.mark.parametrize("chem_comp_id", CHEM_COMP_IDS)
    def test_all_files_created(pipeline_wd: Path, chem_comp_id):
        path = pipeline_wd / chem_comp_id[0] / chem_comp_id
        files = os.listdir(path)
        assert len(files) == 19

        for f in files:
            file_pointer = path / f
            assert file_pointer.stat().st_size > 0

    @staticmethod
    @pytest.mark.parametrize("id_", ["000", "001", "002", "003", "004"])
    def test_images_with_names_created(pipeline_wd: Path, id_):
        """Test if the depictions with names contain certain atom labels
        as expected.
        """
        path = pipeline_wd / id_[0] / id_ / f"{id_}_100_names.svg"
        xml = ET.parse(path)
        path_elements = xml.findall("svg:path", svg_namespace)

        for e in path_elements:
            if re.fullmatch(r"atom-\d+", e.attrib["class"]):
                return  # label is where it should be we can quit

        assert False  # this compound should have had label

    @staticmethod
    @pytest.mark.parametrize(
        "id_,name,alt_name",
        [("001", "H021", "1H02"), ("002", "H121", "1H12"), ("003", "H121", "1H12")],
    )
    def test_correct_atom_naming_in_files(pipeline_wd: Path, id_, name, alt_name):
        """Test if alternate names are used for exported model/ideal
        pdb files.
        """
        alts = [
            pipeline_wd / id_[0] / id_ / f"{id_}_ideal_alt.pdb",
            pipeline_wd / id_[0] / id_ / f"{id_}_model_alt.pdb",
        ]
        regular = [
            pipeline_wd / id_[0] / id_ / f"{id_}_ideal.pdb",
            pipeline_wd / id_[0] / id_ / f"{id_}_model.pdb",
        ]

        for i in alts:
            with open(i) as f:
                str_repr = f.read()
                assert alt_name in str_repr and name not in str_repr

        for i in regular:
            with open(i) as f:
                str_repr = f.read()
                assert name in str_repr and alt_name not in str_repr

    @staticmethod
    @pytest.mark.parametrize("chem_comp_id", CHEM_COMP_IDS)
    def test_cml_files(pipeline_wd: Path, chem_comp_id):
        """Test if the CML file is parsable. Each of the tested compounds
        should contain element with type C.
        """
        path = pipeline_wd / chem_comp_id[0] / chem_comp_id / f"{chem_comp_id}.cml"

        xml_root = ET.parse(str(path)).getroot()
        atoms = xml_root.find("molecule").find("atomArray").findall("atom")

        assert any(a.attrib["elementType"] == "C" for a in atoms)

    @staticmethod
    @pytest.mark.parametrize("chem_comp_id", CHEM_COMP_IDS)
    def test_annotation_file(pipeline_wd: Path, chem_comp_id):
        """Test if the annotation.json is parsable. Each of the tested compounds
        should contain some data.
        """

        path = (
            pipeline_wd
            / chem_comp_id[0]
            / chem_comp_id
            / f"{chem_comp_id}_annotation.json"
        )

        assert path.is_file()

        with open(path, "r") as fp:
            data = json.load(fp)

            assert data
            assert data["atoms"]
            assert data["bonds"]
            assert data["resolution"]
