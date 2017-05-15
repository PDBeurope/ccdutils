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
        """
        tuple of the atom_id's (aka atom names) in the chem_comp

        Returns:
            (str): the atom_id's
        """
        atom_ids = []
        for atom in self.atoms:
            atom_ids.append(atom.atom_id)
        return tuple(atom_ids)

    @property
    def number_atoms(self):
        """
        The number of atoms in the chem_comp

        Returns:
            int: the number of atoms
        """
        return len(self.atoms)

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
        self.read_ccd_from_file_mmcifIO(file_name)

    def read_ccd_from_file_mmcifIO(self, file_name):
        import mmCif.mmcifIO as mmcif
        cif_parser = mmcif.CifFileReader(input='data', preserve_order=True)
        cifObj = cif_parser.read(file_name, output='cif_wrapper')
        chem_comp = list(cifObj.values())[0]
        self.chem_comp_id = chem_comp._chem_comp['id'][0]
        self.chem_comp_name = chem_comp._chem_comp['name'][0]
        atoms=chem_comp._chem_comp_atom
        self.atoms = []
        for atom in atoms:
            atom_id= atom['atom_id']
            pdbx_stereo_config = atom['pdbx_stereo_config']
            ideal_x = float(atom['pdbx_model_Cartn_x_ideal'])
            ideal_y = float(atom['pdbx_model_Cartn_y_ideal'])
            ideal_z = float(atom['pdbx_model_Cartn_z_ideal'])
            this_atom = self.Atom(atom_id=atom_id,
                                  pdbx_stereo_config=pdbx_stereo_config,
                                  xyz_ideal=(ideal_x, ideal_y, ideal_z))
            self.atoms.append(this_atom)
