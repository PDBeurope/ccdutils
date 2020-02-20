#!/usr/bin/env python
# software from PDBe: Protein Data Bank in Europe; https://pdbe.org
#
# Copyright 2018 EMBL - European Bioinformatics Institute
# Licensed under the Apache License, Version 2.0 (the 'License');
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing,
# software distributed under the License is distributed on
# an 'AS IS' BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
# KIND, either express or implied. See the License for the
# specific language governing permissions and limitations
# under the License.

"""module handling UniChem mapping
"""

import requests

# these are the resources we agreed on providing as part of the PDBeChem
agreed_resources = {
    '1': ('ChEMBL', 'ChEMBL'),
    '2': ('DrugBank', 'DrugBank'),
    '4': ('Guide to Pharmacology', 'Guide to Pharmacology'),
    '6': ('KEGG LIGAND', 'KEGG (Kyoto Encyclopedia of Genes and Genomes) Ligand'),
    '7': ('ChEBI', 'ChEBI (Chemical Entities of Biological Interest).'),
    '9': ('ZINC', 'ZINC'),
    '22': ('PubChem', 'PubChem Compounds'),
    '24': ('NMRShiftDB', 'NMRShiftDB'),
    '31': ('BindingDb', 'BindingDb'),
    '33': ('LipidMaps', 'LipidMaps'),
    '36': ('MetaboLights', 'MetaboLights'),
    '37': ('BRENDA', 'BRENDA'),
    '38': ('Rhea', 'Rhea')
}

# list taken from https://www.ebi.ac.uk/unichem/ucquery/listSources
all_resource = {
    '1': ('ChEMBL', 'ChEMBL'),
    '2': ('DrugBank', 'DrugBank'),
    '3': ('PDBe', 'Protein Data Bank in Europe'),
    '4': ('Guide to Pharmacology', 'Guide to Pharmacology'),
    '5': ('PubChem DOTF', 'PubChem ("Drugs of the Future" subset'),
    '6': ('KEGG LIGAND', 'KEGG (Kyoto Encyclopedia of Genes and Genomes) Ligand'),
    '7': ('ChEBI', 'ChEBI (Chemical Entities of Biological Interest).'),
    '8': ('NIH', 'NIH Clinical Collection'),
    '9': ('ZINC', 'ZINC'),
    '10': ('eMolecules', 'eMolecules'),
    '11': ('IBM', 'IBM strategic IP insight platform and the National Institutes of Health'),
    '12': ('atlas', 'Gene Expression Atlas'),
    '14': ('fdasrs', 'FDA/USP Substance Registration System (SRS)'),
    '15': ('SureChEMBL', 'SureChEMBL'),
    '17': ('PharmGKB', 'PharmGKB'),
    '18': ('HMDB', 'Human Metabolome Database (HMDB)'),
    '20': ('Selleck', 'Selleck'),
    '21': ('PubChem TPHARMA', 'PubChem ("Thomson Pharma" subset)'),
    '22': ('PubChem', 'PubChem Compounds'),
    '23': ('Mcule', 'Mcule'),
    '24': ('NMRShiftDB', 'NMRShiftDB'),
    '25': ('LINCS', 'Library of Integrated Network-based Cellular Signatures'),
    '26': ('ACTor', 'ACTor'),
    '27': ('Recon', 'Recon'),
    '28': ('MolPort', 'MolPort'),
    '29': ('Nikkaji', 'Nikkaji'),
    '31': ('BindingDb', 'BindingDb'),
    '32': ('EPA CompTox Dashboard', 'EPA (Environmental Protection Agency) CompTox Dashboard'),
    '33': ('LipidMaps', 'LipidMaps'),
    '34': ('DrugCentral', 'DrugCentral'),
    '35': ('Carotenoid Database', 'Carotenoid Database'),
    '36': ('MetaboLights', 'MetaboLights'),
    '37': ('BRENDA', 'BRENDA'),
    '38': ('Rhea', 'Rhea'),
    '39': ('ChemicalBook', 'ChemicalBook'),
    '41': ('SwissLipids', 'SwissLipids')
}

url_prefix = f'https://www.ebi.ac.uk/unichem/rest/inchikey/'


def get_agreed_unichem_mapping(inchikey):
    """Get the agreed mappings for PDBeChem:
        * ChEMBL, DrugBank, KEGG LIGAND, ChEBI, ZINC, NMRShiftDB, BindingDB
          MetaboLights, BRENDA, Rhea

    Args:
        inchikey (str): InchiKey (preferably from CIF CCD)

    Returns:
        dict of str: str: Resource mapping. Key is a resource name,
            value is the internal identifier.
    """

    return _download_unichem_mapping(inchikey, agreed_resources)


def get_all_unichem_mapping(inchikey):
    """Get all mappings for PDBeChem

    Args:
        inchikey (str): InchiKey (preferably from CIF CCD)

    Returns:
        dict of str: str: Resource mapping. Key is a resource name,
            value is the internal identifier.
    """

    return _download_unichem_mapping(inchikey, all_resource)


def _download_unichem_mapping(inchikey, resources):
    """Internal function to query unichem server.

    Args:
        inchikey (str): InchiKey (preferably from CIF CCD)
        matched_dict (dict of str: str): Internal UniChem mapping of
            external servers.

    Returns:
        dict of str: [str]: Resource mapping. Key is a resource name,
            value is the list of internal identifiers.
    """
    mapping = []
    url = f'{url_prefix}/{inchikey}'

    try:
        headers = {
            'Content-Type': 'application/json'
        }
        r = requests.get(url, headers)
        jsonFile = r.json()

        for entity in jsonFile:
            if entity["src_id"] not in resources:
                continue

            resource = resources[entity["src_id"]]
            key = resource[0]
            mapping.append((key, entity["src_compound_id"]))

        return mapping

    except Exception:
        return mapping
