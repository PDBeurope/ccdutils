#!/usr/bin/env python
# software from PDBe: Protein Data Bank in Europe; https://pdbe.org
#
# Copyright 2018 EMBL - European Bioinformatics Institute
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


import os

import numpy
import requests
from pdbeccdutils.core import ccd_reader
from rdkit import Chem
from rdkit.Chem import AllChem

import logging

logging.getLogger("urllib3").setLevel(logging.WARNING)


def rescale_molecule(path, factor):
    """
    Rescale molecule coords to a given factor

    Args:
        path (str) Path to the molecule to be rescaled.
        factor (float): rescaling factor
    """
    mol = Chem.MolFromMolFile(path, sanitize=True)
    matrix = numpy.zeros((4, 4), float)

    for i in range(3):
        matrix[i, i] = factor
    matrix[3, 3] = 1

    AllChem.TransformMol(mol, matrix)
    Chem.MolToMolFile(mol, path)


def download_template(destination: str, template_id: str, inchikey: str) -> bool:
    """Download 2D layout from the PubChem FTP.

    Args:
        destination (str): Path to the pubchem template
        template_id (str): CCD id of a pubchem template
        inchikey (str): CCD's INCHIKey

    Returns:
        bool: If the download was successful or no.
    """
    pubchem_api = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound"

    if os.path.isfile(destination):
        return False

    inchi_url = f"{pubchem_api}/inchikey/{inchikey}/cids/json"
    response = requests.get(inchi_url)

    if response.status_code != 200:
        return False

    json_file = response.json()
    cid = json_file["IdentifierList"]["CID"][0]

    structure_url = f"{pubchem_api}/cid/{cid}/record/SDF/?record_type=2d&response_type=save&response_basename={template_id}.sdf"
    response = requests.get(structure_url)

    if response.status_code != 200:
        return False

    with open(destination, "wb") as fp:
        fp.write(response.content)

    return True


class PubChemDownloader:
    """
    Toolkit to retrieve 2D layouts from the PubChem database.
    """

    def __init__(self, pubchem_templates: str) -> None:
        if not os.path.isdir(pubchem_templates):
            raise ValueError(f"{pubchem_templates} is not a valid path")

        self.pubchem_templates = pubchem_templates

    def update_ccd_dir(self, components: str):
        """Update 2D images of PDBeChem components which are available
        in the pubchem database

        Args:
            components (str): Path to the directory with components in
                .cif format
        """

        for f in os.listdir(components):
            c = ccd_reader.read_pdb_cif_file(os.path.join(components, f)).component
            self.process_template(c)

    def update_ccd_file(self, ccd: str) -> None:
        """Update 2d images of pdbechem components which are available
        in the pubchem database from CCD files.

        Args:
            ccd (str): Path to the the `.cif` CCD file
        """
        components = ccd_reader.read_pdb_components_file(ccd)

        for i in components.values():
            self.process_template(i.component)

    def process_template(self, component):
        """Process template for a given component. First the component
        is attempted to be downloaded and re-scaled. Since the RDKit
        default depiction has 1.5A single bond size whereas templates
        from pubchem are 1.0A.

        Args:
            component (Component): Component
            destination (str): Path to the pubchem 2D template dir

        Returns:
            bool: whether or not the new template has been processed
        """
        destination = os.path.join(self.pubchem_templates, f"{component.id}.sdf")
        downloaded = download_template(destination, component.id, component.inchikey)

        if downloaded:
            rescale_molecule(destination, 1.5)

        return downloaded
