"""Common fixtures shared among all the tests
"""
from pdbeccdutils.core import ccd_writer
from pdbeccdutils.core.models import ConformerType
from pdbeccdutils.helpers.mol_tools import fix_conformer
import rdkit
from rdkit import Chem
import numpy as np
import os


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
