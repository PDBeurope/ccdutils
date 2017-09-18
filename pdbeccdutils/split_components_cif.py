# software from PDBe: Protein Data Bank in Europe; http://pdbe.org
#
# Copyright 2017 EMBL - European Bioinformatics Institute
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing,
# software distributed under the License is distributed on
# an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
# KIND, either express or implied. See the License for the
# specific language governing permissions and limitations
# under the License.
#
import mmCif.mmcifIO as mmcifIO
import pprint


class SplitComponentsCif(object):
    """deals with splitting the wwPDB Chemical Component Dictionary into individual PdbChemicalComponentsRDKit
    objects for each chem_comp definition.
    """
    def __init__(self, file_name, logger=None):
        cif_parser = mmcifIO.CifFileReader(input='data', preserve_order=True)
        self.cif_dictionary = cif_parser.read(file_name, output='cif_dictionary')
        if self.cif_dictionary is None:
            raise IOError('problem reading file {}'.format(file_name))
        if logger is None:
            self.logger_or_print = pprint.pprint
        else:
            self.logger_or_print = logger.error

    def individual_cif_dictionary(self):
        """
        gives each datablock as an individual (python) dictionary

        Yields:
            each data in the cif file in turn as an individual cif python dictionary.
        """
        for data_block_id, data_block in self.cif_dictionary.items():
            individual_cif_dictionary = {data_block_id: data_block}
            yield individual_cif_dictionary

    @staticmethod
    def write_individual_cif_dictionary(individual_cif_dictionary, file_name):
        """
        writes an individual_cif_dictionary to a file

        Args:
            individual_cif_dictionary: an individual cif dictionary containing  {data_block_id: data_block}
            file_name (str): Name for file
        """
        cfd = mmcifIO.CifFileWriter(file_name, preserve_order=True)
        cfd.write(individual_cif_dictionary)
