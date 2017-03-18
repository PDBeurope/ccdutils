import os

class PdbChemicalComponents(object):
    def __init__(self, file_name=None):
        self.chem_comp_id = None
        self.chem_comp_name = None
        self.number_of_atoms = 0
        self.atom_id = []
        if file_name is not None:
            self.read_ccd_from_cif_file(file_name)

    def read_ccd_from_cif_file(self, file_name):
        """
        reads the ccd from a cif file

        Args:
            file_name (str): the filename

        Returns:
            None
        """
        if not os.path.isfile(file_name):
            raise ValueError('cannot read chemical compoents from %s as file not found' % file_name)
        pass  # TODO write method
