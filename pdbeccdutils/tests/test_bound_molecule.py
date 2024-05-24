"""Common fixtures shared among all the tests
"""

import json
import pytest

import gemmi
from networkx import MultiDiGraph
from pdbeccdutils.tests import tst_utilities
from pdbeccdutils.core import boundmolecule
from pdbeccdutils.core.models import AssemblyResidue, BoundMolecule


test_residues = {
    "first": [
        AssemblyResidue(
            "NAG",
            "F",
            "1",
            None,
            "1",
            "F",
            "",
        ),
        AssemblyResidue(
            "NAG",
            "F",
            "2",
            None,
            "1",
            "F",
            "",
        ),
    ],
    "second": [
        AssemblyResidue(
            "NAG",
            "F_2",
            "2",
            None,
            "1",
            "F",
            "2",
        ),
        AssemblyResidue(
            "NAG",
            "F_2",
            "1",
            None,
            "1",
            "F",
            "2",
        ),
    ],
    "third": [
        AssemblyResidue(
            "NAG",
            "J_2",
            "2",
            None,
            "1",
            "J",
            "2",
        ),
        AssemblyResidue(
            "NAG",
            "F_2",
            "1",
            None,
            "1",
            "F",
            "2",
        ),
    ],
}

test_inputs = {
    "1c4q": {"au_fallback": False},
    "1tqh": {"au_fallback": False},
}


def create_boundmolecule(residues: list[AssemblyResidue]):
    graph = MultiDiGraph()
    for residue in residues:
        graph.add_node(residue)

    bm = BoundMolecule(graph)
    return bm


def test_boundmolecule_equivalence():
    bm_1 = create_boundmolecule(test_residues["first"])
    bm_2 = create_boundmolecule(test_residues["second"])
    bm_3 = create_boundmolecule(test_residues["third"])
    assert bm_1.is_equivalent(bm_2)
    assert not bm_1.is_equivalent(bm_3)


class TestBoundMolecule:
    @pytest.fixture(autouse=True, params=list(test_inputs.keys()))
    def ligand_sites(self, request):
        self.bio_assembly_path = tst_utilities.bio_assembly_filename(
            request.param, test_inputs[request.param]["au_fallback"]
        )
        self.boundmolecules = tst_utilities.bms_filename(request.param)
        self.bio_assembly_block = gemmi.cif.read(self.bio_assembly_path).sole_block()
        self.structure = gemmi.make_structure_from_block(self.bio_assembly_block)
        self.bio_assembly_cat = self.bio_assembly_block.get_mmcif_category_names()
        self.discarded_ligands = ["HOH"]

    def test_parse_ligands_from_nonpoly_scheme(self):
        if "_pdbx_nonpoly_scheme." not in self.bio_assembly_cat:
            pytest.skip("_pdbx_nonpoly_scheme not present")
        else:
            nonpoly_scheme = self.bio_assembly_block.get_mmcif_category(
                "_pdbx_nonpoly_scheme."
            )
            new_nonpoly_scheme = tst_utilities.remove_discarded_ligands(
                nonpoly_scheme, "pdb_mon_id", self.discarded_ligands
            )
            nonpoly_ligands = boundmolecule.parse_ligands_from_nonpoly_scheme(
                nonpoly_scheme, self.discarded_ligands, assembly=True
            )
            assert len(new_nonpoly_scheme["pdb_mon_id"]) == len(
                list(nonpoly_ligands.nodes())
            )

    def test_parse_ligands_from_branch_scheme(self):
        if "_pdbx_branch_scheme." not in self.bio_assembly_cat:
            pytest.skip("_pdbx_branch_scheme not present")
        else:
            branch_scheme = self.bio_assembly_block.get_mmcif_category(
                "_pdbx_branch_scheme."
            )
            new_branch_scheme = tst_utilities.remove_discarded_ligands(
                branch_scheme, "pdb_mon_id", self.discarded_ligands
            )
            branch_ligands = boundmolecule.parse_ligands_from_branch_scheme(
                branch_scheme, self.discarded_ligands, MultiDiGraph(), assembly=True
            )

            assert len(new_branch_scheme["pdb_mon_id"]) == len(
                list(branch_ligands.nodes())
            )

    def test_infer_bound_molecules(self):
        bound_molecules = boundmolecule.infer_bound_molecules(
            self.bio_assembly_path, self.discarded_ligands
        )
        with open(self.boundmolecules, "r") as fh:
            bms_json = json.load(fh)

        assert len(bound_molecules) == len(bms_json["boundMolecules"])
