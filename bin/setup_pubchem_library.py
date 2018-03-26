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

import argparse
from pdbeccdutils.core import structure_reader
from pdbeccdutils.utils import PubChemDownloader

"""Modul to create and maintain up-to-date pubchem library
"""

parser = argparse.ArgumentParser(description='PDBe downloader of pubchem depictions')
parser.add_argument('-components', type=str, help='Path to the component library', required=True)
parser.add_argument('-pubchem_templates', type=str, help='Path to the pubchem templates.',
                    required=True)

config = parser.parse_args()
PubChemDownloader(config.components, config.pubchem_templates).run()
