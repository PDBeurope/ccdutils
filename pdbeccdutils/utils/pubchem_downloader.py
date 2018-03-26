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
import sys
import json
import urllib.request
from rdkit import Chem

from pdbeccdutils.core import structure_reader


class PubChemDownloader:
    """Toolkit to retrieve pubchem 2D depictions from the database
    """

    def __init__(self, components, pubchem_templates):
        self.components = components
        self.pubchem_templates = pubchem_templates
        self.blacklist = list()

    def run(self):
        """Update 2d images of pdbechem components which are available
        in the pubchem database
        """

        print('Querying pubchem database...')
        downloaded = self._download()
        print('Downloaded {} new structures.'.format(downloaded))

    def _download(self):
        """
        Downloads 2D structures of the components and returns a number
        of new structures
        """
        counter = 0
        pubchem_api = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound'
        i = 0
        for file in os.listdir(self.components):
            id = os.path.basename(file).split('.')[0]
            destination = os.path.join(self.pubchem_templates, id + '.sdf')
            counter += 1
            print('{} | new {}'.format(counter, i), end='\r')

            if os.path.isfile(destination):
                continue
            inchikey = structure_reader.read_pdb_cif_file(os.path.join(self.components, file)).component.inchikey

            try:
                inchi_url = '{}/inchikey/{}/cids/json'.format(pubchem_api, inchikey)
                response = urllib.request.urlopen(inchi_url).read().decode('utf-8')
                jsonFile = json.loads(response)
                cid = jsonFile['IdentifierList']['CID'][0]

                structure_url = '{}/cid/{}/record/SDF/?record_type=2d&response_type=save&response_basename={}'.format(pubchem_api, cid, id + '.sdf')
                urllib.request.urlretrieve(structure_url, destination)
                i += 1
            except urllib.request.HTTPError:
                pass
            except urllib.error.URLError:
                pass
        return i
