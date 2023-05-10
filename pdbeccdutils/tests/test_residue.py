import pytest
from pdbeccdutils.core.models import Residue, AssemblyResidue

test_residue_input = {
    "TRS": {
        "name": "TRS",
        "chain": "A",
        "res_id": "1007",
        "ins_code": None,
        "ent_id": "3",
    },
    "FMT": {
        "name": "FMT",
        "chain": "A_2",
        "res_id": "1008",
        "ins_code": None,
        "ent_id": "4",
    },
}

test_assembly_residue_input = {
    "TRS": {
        "name": "TRS",
        "chain": "A",
        "res_id": "1007",
        "ins_code": None,
        "ent_id": "3",
        "orig_chain": "A",
        "operator": "",
    },
    "FMT": {
        "name": "FMT",
        "chain": "A_2",
        "res_id": "1008",
        "ins_code": None,
        "ent_id": "4",
        "orig_chain": "A",
        "operator": "_2",
    },
}


class TestResidue:
    @pytest.fixture(autouse=True, params=list(test_residue_input.keys()))
    def initialize_residue(self, request):
        residue = Residue(
            test_residue_input[request.param]["name"],
            test_residue_input[request.param]["chain"],
            test_residue_input[request.param]["res_id"],
            test_residue_input[request.param]["ins_code"],
            test_residue_input[request.param]["ent_id"],
        )
        self.residue = residue
        self.resname = request.param

    def test_residue_initialization(self):
        assert self.residue.name == test_residue_input[self.resname]["name"]
        assert self.residue.chain == test_residue_input[self.resname]["chain"]
        assert self.residue.res_id == test_residue_input[self.resname]["res_id"]
        assert self.residue.ins_code == ""
        assert self.residue.ent_id == test_residue_input[self.resname]["ent_id"]
        assert (
            self.residue.id
            == test_residue_input[self.resname]["chain"]
            + test_residue_input[self.resname]["res_id"]
            + self.residue.ins_code
        )

    def test_residue_to_dict(self):
        residue_dict = self.residue.to_dict()
        assert residue_dict["id"] == self.residue.id
        assert residue_dict["label_comp_id"] == self.residue.name
        assert residue_dict["auth_asym_id"] == self.residue.chain
        assert residue_dict["auth_seq_id"] == self.residue.res_id
        assert residue_dict["pdbx_PDB_ins_code"] == " "
        assert residue_dict["entity_id"] == self.residue.ent_id

    def test_residue_to_arpeggio(self):
        arpeggio_str = self.residue.to_arpeggio()
        assert (
            arpeggio_str
            == f"/{self.residue.chain}/{self.residue.res_id}{self.residue.ins_code}/"
        )

    def test_residue_equality(self):
        residue_1 = Residue(
            test_residue_input["TRS"]["name"],
            test_residue_input["TRS"]["chain"],
            test_residue_input["TRS"]["res_id"],
            test_residue_input["TRS"]["ins_code"],
            test_residue_input["TRS"]["ent_id"],
        )
        residue_2 = Residue(
            test_residue_input["FMT"]["name"],
            test_residue_input["FMT"]["chain"],
            test_residue_input["FMT"]["res_id"],
            test_residue_input["FMT"]["ins_code"],
            test_residue_input["FMT"]["ent_id"],
        )
        residue_3 = Residue(
            test_residue_input["TRS"]["name"],
            test_residue_input["TRS"]["chain"],
            test_residue_input["TRS"]["res_id"],
            test_residue_input["TRS"]["ins_code"],
            test_residue_input["TRS"]["ent_id"],
        )

        assert residue_1 == residue_3
        assert residue_1 != residue_2


class TestAssemblyResidue:
    @pytest.fixture(autouse=True, params=list(test_residue_input.keys()))
    def initialize_residue(self, request):
        residue = AssemblyResidue(
            test_assembly_residue_input[request.param]["name"],
            test_assembly_residue_input[request.param]["chain"],
            test_assembly_residue_input[request.param]["res_id"],
            test_assembly_residue_input[request.param]["ins_code"],
            test_assembly_residue_input[request.param]["ent_id"],
            test_assembly_residue_input[request.param]["orig_chain"],
            test_assembly_residue_input[request.param]["operator"],
        )
        self.residue = residue
        self.resname = request.param

    def test_residue_initialization(self):
        assert self.residue.name == test_assembly_residue_input[self.resname]["name"]
        assert self.residue.chain == test_assembly_residue_input[self.resname]["chain"]
        assert (
            self.residue.res_id == test_assembly_residue_input[self.resname]["res_id"]
        )
        assert self.residue.ins_code == ""
        assert (
            self.residue.ent_id == test_assembly_residue_input[self.resname]["ent_id"]
        )
        assert (
            self.residue.id
            == test_assembly_residue_input[self.resname]["chain"]
            + test_assembly_residue_input[self.resname]["res_id"]
            + self.residue.ins_code
        )
        assert (
            self.residue.orig_chain
            == test_assembly_residue_input[self.resname]["orig_chain"]
        )
        assert (
            self.residue.operator
            == test_assembly_residue_input[self.resname]["operator"]
        )

    def test_residue_to_dict(self):
        residue_dict = self.residue.to_dict()
        assert residue_dict["id"] == self.residue.id
        assert residue_dict["label_comp_id"] == self.residue.name
        assert residue_dict["auth_asym_id"] == self.residue.chain
        assert residue_dict["auth_seq_id"] == self.residue.res_id
        assert residue_dict["pdbx_PDB_ins_code"] == " "
        assert residue_dict["entity_id"] == self.residue.ent_id
        assert residue_dict["orig_auth_asym_id"] == self.residue.orig_chain
        assert residue_dict["operator"] == self.residue.operator

    def test_residue_to_arpeggio(self):
        arpeggio_str = self.residue.to_arpeggio()
        assert (
            arpeggio_str
            == f"/{self.residue.chain}/{self.residue.res_id}{self.residue.ins_code}/"
        )

    def test_residue_equality(self):
        residue_1 = AssemblyResidue(
            test_assembly_residue_input["TRS"]["name"],
            test_assembly_residue_input["TRS"]["chain"],
            test_assembly_residue_input["TRS"]["res_id"],
            test_assembly_residue_input["TRS"]["ins_code"],
            test_assembly_residue_input["TRS"]["ent_id"],
            test_assembly_residue_input["TRS"]["orig_chain"],
            test_assembly_residue_input["TRS"]["operator"],
        )
        residue_2 = AssemblyResidue(
            test_assembly_residue_input["FMT"]["name"],
            test_assembly_residue_input["FMT"]["chain"],
            test_assembly_residue_input["FMT"]["res_id"],
            test_assembly_residue_input["FMT"]["ins_code"],
            test_assembly_residue_input["FMT"]["ent_id"],
            test_assembly_residue_input["FMT"]["orig_chain"],
            test_assembly_residue_input["FMT"]["operator"],
        )
        residue_3 = AssemblyResidue(
            test_assembly_residue_input["TRS"]["name"],
            test_assembly_residue_input["TRS"]["chain"],
            test_assembly_residue_input["TRS"]["res_id"],
            test_assembly_residue_input["TRS"]["ins_code"],
            test_assembly_residue_input["TRS"]["ent_id"],
            test_assembly_residue_input["TRS"]["orig_chain"],
            test_assembly_residue_input["TRS"]["operator"],
        )

        assert residue_1 == residue_3
        assert residue_1 != residue_2
