import os
from collections import namedtuple


class PdbChemicalComponents(object):

    def __init__(self, file_name=None):
        self.chem_comp_id = None
        self.chem_comp_name = None
        self.Atom = namedtuple('Atom', 'atom_id pdbx_stereo_config xyz_ideal')
        self.atoms = []
        if file_name is not None:
            self.read_ccd_from_cif_file(file_name)

    @property
    def atom_ids(self):
        atom_ids = []
        for atom in self.atoms:
            atom_ids.append(atom.atom_id)
        return atom_ids

    @property
    def number_atoms(self):
        return len(self.atoms)


    def read_ccd_from_cif_file(self, file_name):
        """
        reads the ccd from a cif fileC

        Args:
            file_name (str): the filename

        Returns:
            None
        """
        if not os.path.isfile(file_name):
            raise ValueError('cannot read chemical compoents from %s as file not found' % file_name)
        pass  # TODO write method

    def test_hard_code_CMO(self):
        """
        stub to produce a hard coded carbon monoxide ccd object for development idea/testing

        Returns:
            None
        """
        # _chem_comp.id                                    CMO
        self.chem_comp_id = 'CMO'
        # _chem_comp.name                                  "CARBON MONOXIDE"
        self.chem_comp_name = 'CARBON MONOXIDE'
        # CMO C C C -1 1 N N N -0.296 8.526 17.112 0.607  0.000 0.000 C CMO 1
        # CMO O O O 1  1 N N N 0.023  7.997 18.053 -0.600 0.000 0.000 O CMO 2
        this_atom = self.Atom(atom_id='C', pdbx_stereo_config='N', xyz_ideal=(0.607, 0.000, 0.000))
        self.atoms.append(this_atom)
        this_atom = self.Atom(atom_id='O', pdbx_stereo_config='N', xyz_ideal=(-0.600, 0.000, 0.000))
        self.atoms.append(this_atom)
