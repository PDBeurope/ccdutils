import logging
from pathlib import Path

from gemmi import cif
from rdkit import Chem

from pdbeccdutils.core import ccd_reader
from pdbeccdutils.tests.tst_utilities import cif_filename


def _components_cif(tmp_path):
    path = tmp_path / "components.cif"
    ccds = [Path(cif_filename("00O")), Path(cif_filename("007"))]
    path.write_text("\n".join(ccd.read_text() for ccd in ccds))
    return path


def test_read_pdb_components_file_parses_all_for_empty_include(tmp_path):
    path = _components_cif(tmp_path)

    for include in (None, []):
        result = ccd_reader.read_pdb_components_file(
            str(path), sanitize=False, include=include
        )

        assert list(result) == ["00O", "007"]


def test_read_pdb_components_file_uses_include_order_and_skips_missing(
    tmp_path, caplog
):
    path = _components_cif(tmp_path)

    with caplog.at_level(logging.WARNING):
        result = ccd_reader.read_pdb_components_file(
            str(path), sanitize=False, include=["007", "MISSING", "00O", "007"]
        )

    assert list(result) == ["007", "00O"]
    assert "Data block MISSING not found" in caplog.text


def test_parse_pdb_bonds_uses_atom_id_lookup():
    block = cif.read_string(
        """
        data_TST
        loop_
        _chem_comp_atom.atom_id
        A
        B
        loop_
        _chem_comp_bond.atom_id_1
        _chem_comp_bond.atom_id_2
        _chem_comp_bond.value_order
        A B SING
        """
    ).sole_block()
    mol = Chem.RWMol()
    mol.AddAtom(Chem.Atom("C"))
    mol.AddAtom(Chem.Atom("O"))
    errors = []

    ccd_reader._parse_pdb_bonds(mol, block, errors)

    assert errors == []
    assert mol.GetNumBonds() == 1
    assert mol.GetBondBetweenAtoms(0, 1) is not None


def test_parse_pdb_bonds_missing_atom_1_does_not_mask_error():
    block = cif.read_string(
        """
        data_TST
        loop_
        _chem_comp_atom.atom_id
        A
        loop_
        _chem_comp_bond.atom_id_1
        _chem_comp_bond.atom_id_2
        _chem_comp_bond.value_order
        ? A SING
        """
    ).sole_block()
    mol = Chem.RWMol()
    mol.AddAtom(Chem.Atom("C"))
    errors = []

    ccd_reader._parse_pdb_bonds(mol, block, errors)

    assert errors == [
        f"Missing atom in 0 entry in _chem_comp_bond"
    ]
    assert mol.GetNumBonds() == 0