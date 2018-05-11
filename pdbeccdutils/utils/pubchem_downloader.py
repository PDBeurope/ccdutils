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
import sys
import urllib.request

from rdkit import Chem

from pdbeccdutils.core import structure_reader as sr
from pdbeccdutils.core import ReleaseStatus


class PubChemDownloader:
    """
    Toolkit to retrieve pubchem 2D depictions from the database
    """

    def __init__(self, pubchem_templates):
        if not os.path.isdir(pubchem_templates):
            raise ValueError(pubchem_templates + ' is not a valid path')

        self.pubchem_templates = pubchem_templates

    def update_ccd_dir(self, components):
        """Update 2d images of pdbechem components which are available
        in the pubchem database

        Args:
            components (str): Path to the directory with components in
                .cif format
        """

        print('Querying pubchem database...')
        counter = 0
        downloaded = 0

        for f in os.listdir(components):
            id = f.split('.')[0]
            c = sr.read_pdb_cif_file(os.path.join(components, f)).component
            destination = os.path.join(self.pubchem_templates, id + '.sdf')
            success = self.download_template(c)

            if success:
                downloaded += 1
            counter += 1

            print('{} | new {}'.format(counter, downloaded), end='\r')

        print('Downloaded {} new structures.'.format(downloaded))

    def update_ccd_file(self, ccd):
        """Update 2d images of pdbechem components which are available
        in the pubchem database from CCD files. Only released components
        are downloaded

        Args:
            ccd (str): Path to the the z`.cif CCD file
        """

        print('Querying pubchem database...')
        counter = 0
        downloaded = 0
        components = sr.read_pdb_components_file(ccd)

        for k, v in components.items():
            destination = os.path.join(self.pubchem_templates, k + '.sdf')
            if v.component.released != ReleaseStatus.OBS:
                continue

            success = self.download_template(v.component)

            if success:
                downloaded += 1
            counter += 1

            print('{} | new {}'.format(counter, downloaded), end='\r')

        print('Downloaded {} new structures.'.format(downloaded))

    def download_template(self, component):
        """Downloads 2D structure of a given component

        Args:
            component (pdbeccdutils.core.Component): Component
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
