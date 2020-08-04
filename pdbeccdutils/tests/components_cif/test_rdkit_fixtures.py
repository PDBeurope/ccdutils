import rdkit
import numpy as np
from pdbeccdutils.helpers.rdkit_fixtures import fix_conformer

class TestRDKitFixtures:
    @staticmethod
    def test_conformer_is_broken_ion():
        mol = rdkit.Chem.RWMol()
        atom =  rdkit.Chem.Atom('H')

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

    @staticmethod
    def test_conformer_has_broken_atom():
        mol = rdkit.Chem.RWMol()
        o =  rdkit.Chem.Atom('O')
        h =  rdkit.Chem.Atom('H')

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
