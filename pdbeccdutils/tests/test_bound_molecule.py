"""Common fixtures shared among all the tests
"""
import os
import json
import pytest

import gemmi
import rdkit
import numpy as np
from rdkit import Chem
from networkx import MultiDiGraph

from pdbeccdutils.core import ccd_writer
from pdbeccdutils.core.models import ConformerType
from pdbeccdutils.helpers.mol_tools import fix_conformer
from pdbeccdutils.tests import tst_utilities
from pdbeccdutils.core import boundmolecule


test_inputs = {
    "1c4q": {"au_fallback": False},
    "1tqh": {"au_fallback": False},
    "4xkw": {"au_fallback": False},
}


def test_write_bound_molecule(boundMolecules, tmpdir_factory):
    wd = tmpdir_factory.mktemp("bm_test")
    for bm in boundMolecules:
        for remove_hs in True, False:
            suffix = f'{("" if remove_hs else "H")}'
            sdf_file = os.path.join(wd, f"{bm.id}_{suffix}.sdf")
            conf_type = ConformerType.Model
            ccd_writer.write_molecule(
                path=sdf_file,
                component=bm,
                conf_type=conf_type,
                remove_hs=remove_hs,
            )
            rdkit_mol = bm.mol_no_h if remove_hs else bm.mol
            assert os.path.isfile(sdf_file)
            assert os.path.getsize(sdf_file) > 0
            mol = Chem.MolFromMolFile(sdf_file, sanitize=False)
            assert isinstance(mol, Chem.rdchem.Mol)
            assert mol.GetNumAtoms() == rdkit_mol.GetNumAtoms()


def test_bound_molecule_has_degenerate_conformer(boundMolecules):
    for bm in boundMolecules:
        assert not bm.has_degenerated_conformer(ConformerType.Model)


def test_bound_molecule_conformer_is_broken_ion():
    mol = rdkit.Chem.RWMol()
    atom = rdkit.Chem.Atom("H")
    mol.AddAtom(atom)
    conformer = rdkit.Chem.Conformer(1)
    atom_position = rdkit.Chem.rdGeometry.Point3D(np.NaN, np.NaN, np.NaN)
    conformer.SetAtomPosition(0, atom_position)
    mol.AddConformer(conformer, assignId=True)
    m = mol.GetMol()
    c = m.GetConformer(0)
    fix_conformer(c)
    assert c.GetAtomPosition(0).x == 0.0
    assert c.GetAtomPosition(0).y == 0.0
    assert c.GetAtomPosition(0).z == 0.0


def test_bound_molecule_conformer_has_broken_atom():
    mol = rdkit.Chem.RWMol()
    o = rdkit.Chem.Atom("O")
    h = rdkit.Chem.Atom("H")
    mol.AddAtom(o)
    mol.AddAtom(h)
    mol.AddBond(0, 1, rdkit.Chem.rdchem.BondType(1))
    conformer = rdkit.Chem.Conformer(1)
    o_position = rdkit.Chem.rdGeometry.Point3D(1, 2, 3)
    h_position = rdkit.Chem.rdGeometry.Point3D(np.NaN, np.NaN, np.NaN)
    conformer.SetAtomPosition(0, o_position)
    conformer.SetAtomPosition(1, h_position)
    mol.AddConformer(conformer, assignId=True)
    m = mol.GetMol()
    c = m.GetConformer(0)
    fix_conformer(c)
    assert c.GetAtomPosition(0).x != 0.0
    assert c.GetAtomPosition(0).y != 0.0
    assert c.GetAtomPosition(0).z != 0.0
    assert c.GetAtomPosition(1).x == 0.0
    assert c.GetAtomPosition(1).y == 0.0
    assert c.GetAtomPosition(1).z == 0.0


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
