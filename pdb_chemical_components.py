import os

class PdbChemicalComponents(object):
    def __init__(self):
        pass

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

    def number_atoms(self):
        """
        gives the number of atoms in the ccd

        Returns:
            int: the number of atoms
        """
        return 0  # TODO write method
