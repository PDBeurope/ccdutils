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

import json
import os
import urllib.request

from pdbeccdutils.core import ccd_reader as sr


class PubChemDownloader:
    """
    Toolkit to retrieve 2D layouts from the PubChem database.
    """

    def __init__(self, pubchem_templates):
        if not os.path.isdir(pubchem_templates):
            raise ValueError(pubchem_templates + ' is not a valid path')

        self.pubchem_templates = pubchem_templates

    def update_ccd_dir(self, components):
        """Update 2D images of pdbechem components which are available
        in the pubchem database

        Args:
            components (str): Path to the directory with components in
                .cif format
        """

        for f in os.listdir(components):
            c = sr.read_pdb_cif_file(os.path.join(components, f)).component
            self.download_template(c)

    def update_ccd_file(self, ccd):
        """Update 2d images of pdbechem components which are available
        in the pubchem database from CCD files.

        Args:
            ccd (str): Path to the the `.cif` CCD file
        """
        components = sr.read_pdb_components_file(ccd)

        for k, v in components.items():
            self.download_template(v.component)

    def download_template(self, component):
        """Downloads 2D layout of a given component

        Args:
            component (pdbeccdutils.core.component.Component): Component
            destination (str): Path to the pubchem 2D template dir

        Returns:
            bool: whether or not the new template has been downloaded
        """
        pubchem_api = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound'

        template = os.path.join(self.pubchem_templates, component.id + '.sdf')
        if os.path.isfile(template):
            return False

        inchikey = component.inchikey

        try:
            inchi_url = '{}/inchikey/{}/cids/json'.format(pubchem_api, inchikey)
            response = urllib.request.urlopen(inchi_url).read().decode('utf-8')
            jsonFile = json.loads(response)
            cid = jsonFile['IdentifierList']['CID'][0]

            structure_url = '{}/cid/{}/record/SDF/?record_type=2d&response_type=save&response_basename={}'.format(pubchem_api, cid, component.id + '.sdf')
            urllib.request.urlretrieve(structure_url, template)
            return True

        except urllib.request.HTTPError:
            return False
        except urllib.error.URLError:
            return False
