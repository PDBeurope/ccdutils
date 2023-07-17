import pytest
import os
from pdbeccdutils.core import clc_writer
from pdbeccdutils.core.component import Component
from rdkit import Chem
from gemmi import cif

must_have_categories = [
    "_chem_comp.",
    "_chem_comp_atom.",
    "_chem_comp_bond.",
    "_pdbx_chem_comp_descriptor.",
]


class TestBMWriter:
    @staticmethod
    @pytest.mark.parametrize("rem_hs", [True, False])
    def test_cif_write(component_globotriose: Component, tmpdir, rem_hs):
        path = tmpdir.join(f"{component_globotriose.id}.cif")
        to_check = must_have_categories.copy()

        clc_writer.write_molecule(str(path), component_globotriose, remove_hs=rem_hs)
        cif_block = cif.read(str(path)).sole_block()

        assert cif_block
        assert component_globotriose.id == cif_block.name
        for c in to_check:
            assert c in cif_block.get_mmcif_category_names()

    @staticmethod
    @pytest.mark.parametrize("rem_hs", [True, False])
    def test_pdb_write(component_globotriose, tmpdir_factory, rem_hs):
        wd = tmpdir_factory.mktemp("pdb_test")
        suffix = "" if rem_hs else "H"
        pdb_file = os.path.join(wd, f"{component_globotriose.id}{suffix}.pdb")
        clc_writer.write_molecule(
            pdb_file,
            component_globotriose,
            remove_hs=rem_hs,
        )
        rdkit_mol = (
            component_globotriose.mol_no_h if rem_hs else component_globotriose.mol
        )

        assert os.path.isfile(pdb_file)
        assert os.path.getsize(pdb_file) > 0

        mol = Chem.MolFromPDBFile(
            pdb_file,
            removeHs=False,
            sanitize=False,
        )
        assert isinstance(mol, Chem.rdchem.Mol)
        assert mol.GetNumAtoms() == rdkit_mol.GetNumAtoms()
