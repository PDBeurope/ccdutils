import os
import pytest

from gemmi import cif
from pdbeccdutils.helpers import cif_tools
from pdbeccdutils.tests import tst_utilities


test_inputs = ["1c4q", "1tqh"]


class TestFixUpdatedMmcif:
    @pytest.fixture(autouse=True, params=test_inputs)
    def run_fix_updated_mmcif(self, request, tmpdir_factory):
        wd = tmpdir_factory.mktemp("fix_updated_mmcif_test")
        clean_mmcif = tst_utilities.updated_mmcif_filename(request.param)
        self.fixed_mmcif = os.path.join(wd, f"{request.param}_processed.cif")
        cif_tools.fix_updated_mmcif(clean_mmcif, self.fixed_mmcif)

    def test_fixed_mmcif_file_generated(self):
        assert os.path.isfile(self.fixed_mmcif)
        assert os.path.getsize(self.fixed_mmcif) > 0

    def test_multiple_models_removed(self):
        cif_block = cif.read(self.fixed_mmcif).sole_block()
        atom_site = cif_block.get_mmcif_category("_atom_site.")
        assert len(set(atom_site["pdbx_PDB_model_num"])) == 1

    def test_alt_loc_removed(self):
        cif_block = cif.read(self.fixed_mmcif).sole_block()
        atom_site = cif_block.get_mmcif_category("_atom_site.")
        for alt_id in atom_site["label_alt_id"]:
            assert not alt_id

    def test_nonpoly_scheme_updated(self):
        cif_block = cif.read(self.fixed_mmcif).sole_block()
        cif_block_cat = cif_block.get_mmcif_category_names()
        if "_pdbx_nonpoly_scheme" in cif_block_cat:
            nonpoly_scheme = cif_block_cat.get_mmcif_category("_pdbx_nonpoly_scheme.")
            nonpoly_entity_ids = set(nonpoly_scheme["entity_id"])
            output_st_residues = [
                (res.entity_id, res.subchain, res.name, str(res.seqid.num))
                for chain in self.output_st[0]
                for res in chain
                if res.entity_id in nonpoly_entity_ids
            ]
            for i in range(len(nonpoly_scheme["entity_id"])):
                nonpoly_scheme_residue = (
                    nonpoly_scheme["entity_id"][i],
                    nonpoly_scheme["asym_id"][
                        i
                    ],  # _pdbx_nonpoly_scheme.asym_id maps to _atom_site.label_asym_id
                    nonpoly_scheme["pdb_mon_id"][i],
                    nonpoly_scheme["pdb_seq_num"][i],
                )
                assert nonpoly_scheme_residue in output_st_residues

    def test_branch_scheme_updated(self):
        cif_block = cif.read(self.fixed_mmcif).sole_block()
        cif_block_cat = cif_block.get_mmcif_category_names()
        if "_pdbx_branch_scheme" in cif_block_cat:
            branch_scheme = cif_block_cat.get_mmcif_category("_pdbx_branch_scheme.")
            branch_entity_ids = set(branch_scheme["entity_id"])
            branch_chains = set(branch_scheme["pdb_asym_id"])
            # chain number from gemmi maps to _atom_site.auth_asym_id
            output_st_residues = [
                (res.entity_id, chain, res.name, str(res.seqid.num))
                for chain in branch_chains
                for res in self.output_st[0][chain]
                if res.entity_id in branch_entity_ids
            ]

            for i in range(len(branch_scheme["entity_id"])):
                branch_scheme_residue = (
                    branch_scheme["entity_id"][i],
                    branch_scheme["pdb_asym_id"][
                        i
                    ],  # _pdbx_branch_scheme.pdb_asym_id maps to _atom_site.auth_asym_id
                    branch_scheme["pdb_mon_id"][i],
                    branch_scheme["num"][i],
                )
                assert branch_scheme_residue in output_st_residues
