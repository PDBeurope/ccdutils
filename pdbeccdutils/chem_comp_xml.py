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
from lxml import etree
from pdbeccdutils.fragment_library import FragmentLibrary


class ChemCompXMl(object):
    """deals with creating records/file with information for the file
    http://ftp.ebi.ac.uk/pub/databases/msd/pdbechem/chem.xml
    and also
    http://ftp.ebi.ac.uk/pub/databases/msd/pdbechem/chem_comp.list
    """

    def __init__(self, override_fragment_library_file_path=None):
        self.top = etree.Element('chemCompList')
        self.fragment_library = FragmentLibrary(override_fragment_library_file_path)
        self.chem_comp_id_list = []

    def store_ccd(self, pdb_ccd):
        """
        stores a chemical component in the chem.xml and the chem_comp_id in the list.

        Args:
            pdb_ccd: a PdbChemicalComponents for a chemical component
        """
        chem_comp_id = pdb_ccd.chem_comp_id
        self.chem_comp_id_list.append(chem_comp_id)
        chem_comp = etree.SubElement(self.top, 'chemComp')
        this_id = etree.SubElement(chem_comp, 'id')
        this_id.text = chem_comp_id
        name = etree.SubElement(chem_comp, 'name')
        name.text = pdb_ccd.chem_comp_name
        formula = etree.SubElement(chem_comp, 'formula')
        formula.text = pdb_ccd.chem_comp_formula
        systematic_name_text = pdb_ccd.systematic_name_openeye
        if systematic_name_text is not None:
            systematic_name = etree.SubElement(chem_comp, 'systematicName')
            systematic_name.text = systematic_name_text
        stereo_smiles = etree.SubElement(chem_comp, 'stereoSmiles')
        stereo_smiles.text = pdb_ccd.smiles_canonical_cactvs
        non_stereo_smiles = etree.SubElement(chem_comp, 'nonStereoSmiles')
        non_stereo_smiles.text = pdb_ccd.smiles_cactvs
        inchi = etree.SubElement(chem_comp, 'InChi')
        inchi.text = pdb_ccd.inchi
        inchikey = etree.SubElement(chem_comp, 'InChIKey')
        inchikey.text = pdb_ccd.inchikey
        fragments = self.fragment_library.fragments_for_pdb_chemical_components_rdkit(pdb_ccd)
        if len(fragments.keys()) > 0:
            for fragment_name, matches in sorted(fragments.items()):
                id_numb = 1
                for match in matches:
                    fragment = etree.SubElement(chem_comp, 'fragment', name=fragment_name, id=str(id_numb))
                    id_numb += 1
                    for atom_name in match:
                        atom_id = etree.SubElement(fragment, 'atom_id')
                        atom_id.text = atom_name

    def to_string(self, pretty_print=True):
        """
        produces a string representation of the chem_comp.xml

        Returns:
            str: a string with the chem_comp.xml pretty printed
        """
        xml_string = etree.tostring(self.top, pretty_print=pretty_print)
        xml_string = xml_string.decode('utf-8')  # needed for python3
        return xml_string

    def to_file(self, chem_dot_xml_file_name):
        """
        writes the chem_comp.xml aka chem.xml file

        Args:
            chem_dot_xml_file_name (str): file name to be written

        """
        with open(chem_dot_xml_file_name, 'w') as chem_dot_xml_file:
            chem_dot_xml_file.write(self.to_string())

    def chem_comp_id_list_to_string(self):
        """
        provides the list of chem_comp_ids as a string with new lines

        Returns:
            str: the list of chem_comp_ids with a \n after each one

        """
        return '\n'.join(self.chem_comp_id_list + [''])

    def chem_comp_id_list_to_file(self, file_name):
        """
        writes out the chem_comp.list file

        Args:
            file_name (str): the name for the file.
        """
        with open(file_name, 'w') as this_file:
            this_file.write(self.chem_comp_id_list_to_string())
