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
import collections
import logging
import os
import tempfile
from ccdc import io, conformer
from pdb_chemical_components_rdkit import PdbChemicalComponentsRDKit


class PdbCCDMogul(object):
    """ run Mogul on PDB CCD"""
    def __init__(self, file_name=None):
        self.settings_bond_few_hits_threshold = 5
        self.settings_rfactor_filter = '<5%'
        self.settings_generalisation = False
        self.settings_summary_from_mogul = None
        logging.debug('initialize PdbCCD with cif file {}'.format(file_name))
        self.pdb_ccd_rdkit = PdbChemicalComponentsRDKit(file_name=file_name)
        self.store_bonds = []
        self.store_angles = []
        self.store_torsions = []
        self.store_rings = []
        self.classify_bonds = []

    def run_mogul(self):
        """
        Runs the Mogul analysis, storing results in self.store_* lists

        Returns:
            None
        """
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
        engine.settings.generalisation = self.settings_generalisation
        engine.settings.rfactor_filter = self.settings_rfactor_filter
        engine.settings.bond.few_hits_threshold = self.settings_bond_few_hits_threshold
        self.settings_summary_from_mogul = engine.settings.summary()
        logging.debug('engine.settings.summary()=\n{}'.format(self.settings_summary_from_mogul))
        geometry_analysed_molecule = engine.analyse_molecule(molecule)
        logging.debug('number of Mogul analysed bonds={}'.format(len(geometry_analysed_molecule.analysed_bonds)))
        self.store_observation(geometry_analysed_molecule, 'bond')
        self.store_observation(geometry_analysed_molecule, 'angle')
        self.store_observation(geometry_analysed_molecule, 'torsion')
        self.store_observation(geometry_analysed_molecule, 'ring')
        os.remove(sdf_temp)


    def store_observation(self, geometry_analysed_molecule, observation_type):
        if observation_type == 'bond':
            analysed = geometry_analysed_molecule.analysed_bonds
            place_in = self.store_bonds
            hist_max = 4.0
        elif observation_type == 'angle':
            analysed = geometry_analysed_molecule.analysed_angles
            place_in = self.store_angles
            hist_max = 180
        elif observation_type == 'torsion':
            analysed = geometry_analysed_molecule.analysed_torsions
            place_in = self.store_torsions
            hist_max = 180
        elif observation_type == 'ring':
            analysed = geometry_analysed_molecule.analysed_rings
            place_in = self.store_rings
            hist_max = 180
        else:
            raise RuntimeError('unrecognized observation_type={}'.format(observation_type))
        for thing in analysed:
            store = collections.OrderedDict()
            store['indices'] = thing.atom_indices
            atom_ids = []
            for index in thing.atom_indices:
                atom_id = self.pdb_ccd_rdkit.atom_ids[index]
                atom_ids.append(atom_id)
            store['atoms_ids'] = atom_ids
            for key in ['classification', 'd_min', 'local_density', 'lower_quartile', 'maximum',
                        'mean', 'median', 'minimum', 'nhits', 'standard_deviation', 'type',
                        'unusual', 'upper_quartile', 'value', 'z_score']:
                store[key] = getattr(thing, key)
            store['histogram'] = thing.histogram(minimum=0.0, maximum=hist_max)
            store['hist_max'] = hist_max
            store_nt = collections.namedtuple('stored_mogul_' + observation_type, store.keys())(**store)
            logging.debug(store_nt)
            place_in.append(store_nt)

    def classify_observation(self, observation_type):
        """
        classifies Mogul results using own analysis.

        Args:
            observation_type (str): one of 'bond', 'angle', 'torsion' or 'ring'

        Returns:
           None

        Notes:
            initial idea take self.store_bonds

            recalculate Z using a minimum s.d. of 0.010.

            reorder according to z-score and add classification for outlier
            "outlier": Z > 5 purple
            "unusual":  2 >= Z < 5 violet
            "ok" Z < 2: green
            "too few hits": not enough Mogul hits to classify
            store result in self.classify_bonds

            for bonds

        """
        pass  # TODO to be written!
