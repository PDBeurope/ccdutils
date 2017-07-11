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
import logging
import os
import tempfile
from ccdc import io, conformer
from pdb_chemical_components_rdkit import PdbChemicalComponentsRDKit

class PdbCCDMogul(object):
    """ run Mogul on PDB CCD"""
    def __init__(self, file_name=None):
        logging.debug('initialize PdbCCD with cif file {}'.format(file_name))
        self.pdb_ccd_rdkit = PdbChemicalComponentsRDKit(file_name=file_name)
        # write out the molecule as an sdf file to a temporary directory
        __, sdf_temp = tempfile.mkstemp(suffix='.sdf')
        logging.debug('load into PdbChemicalComponentsRDKit and write out temporary sdf file "{}"'.format(sdf_temp))
        self.pdb_ccd_rdkit.sdf_file_or_string(file_name=sdf_temp)
        if not os.path.isfile(sdf_temp) or os.path.getsize(sdf_temp) == 0:
            raise RuntimeError('cannot write out sdf file')
        mol_reader = io.MoleculeReader(sdf_temp)
        molecule = mol_reader[0]
        molecule.standardise_aromatic_bonds()
        molecule.standardise_delocalised_bonds()
        logging.debug('CSD smiles string {}'.format(molecule.smiles))
        engine = conformer.GeometryAnalyser()
        geometry_analysed_molecule = engine.analyse_molecule(molecule)
        logging.debug('number of Mogul analysed bonds={}'.format(len(geometry_analysed_molecule.analysed_bonds)))
        self.analysed_bonds = geometry_analysed_molecule.analysed_bonds
        self.analysed_angles =  geometry_analysed_molecule.analysed_angles
        self.analysed_torsions = geometry_analysed_molecule.analysed_torsions
        self.analysed_rings = geometry_analysed_molecule.analysed_rings
        os.remove(sdf_temp)






