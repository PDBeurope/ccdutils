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
from yattag import Doc

ANGSTROM = '&Aring;'
SIGMA = '&sigma;'
CLASSIFICATION_COLOR = {'outlier': (215. / 255., 48. / 255., 39. / 255.),  # blood orange
                        'very-unusual': (252. / 255., 141. / 255., 89. / 255.),  # mid orange
                        'unusual': (254. / 255., 224. / 255., 144. / 255.),  # yellow/orange
                        'common': (145. / 255., 191. / 255., 219. / 255.),  # mid blue
                        'very-common': (69. / 255., 117. / 255., 180. / 255.)  # blue
                        }
JQUERY_SCRIPT = '''
<script src="https://ajax.googleapis.com/ajax/libs/jquery/3.2.1/jquery.min.js"></script>
<script>
$(document).ready(function(){
    $("#bond_details").hide();
    $("#angle_details").hide();
    $("#torsion_details").hide();
    $("#ring_details").hide();
    $(".toggle").click(function(){
        var details = "#".concat($(this).val(), "_details")
        var show_button = "#".concat($(this).val(), "_show_button")
        $(details).toggle();
        $(show_button).toggle();
    });
});
</script>'''
PIXELS_X = 500
PIXELS_Y = 300


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
        self.classify_observation('bond')


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

            recalculate Z* using a minimum s.d. of 0.010 - store in Zstar
            reorder according to z-score and add classification for outlier
            "outlier": Z > 5 purple
            "unusual":  2 >= Z < 5 violet
            "ok" Z < 2: green
            "too few hits": not enough Mogul hits to classify
            store result in self.classify_bonds.Z
        """
        if observation_type == 'bond':
            work_from = self.store_bonds
            place_in = self.classify_bonds
            sd_min = 0.010
            few_hits_threshold = self.settings_bond_few_hits_threshold
        else:
            raise RuntimeError('unrecognized observation_type={}'.format(observation_type))

        for thing in work_from:
            if thing.nhits == 0:
                zstar = None
            else:
                zstar = (thing.value - thing.mean)/max(thing.standard_deviation,sd_min)
            if thing.nhits >= few_hits_threshold:
                zorder = abs(zstar)
            elif thing.nhits != 0:
                zorder = -100. + zstar
            else:
                zorder = -101.
            logging.debug('zstar={} zorder={}'.format(zstar, zorder))
            if zorder > 5.:
                classification = 'outlier'
            elif zorder > 3.5:
                classification = 'very-unusual'
            elif zorder > 2.0:
                classification = 'unusual'
            elif zorder > 1.:
                classification = 'common'
            elif zorder > 0.:
                classification = 'very-common'
            else:
                classification = 'too few hits'
            store = thing._asdict()
            store['zstar'] = zstar
            store['zorder'] = zorder
            store['classification'] = classification
            store_nt = collections.namedtuple('stored_mogul_' + observation_type, store.keys())(**store)
            logging.debug(store_nt)
            place_in.append(store_nt)


    def prepare_html(self):
        doc, tag, text, line = Doc().ttl()

        chem_comp_id = self.pdb_ccd_rdkit.chem_comp_id
        chem_comp_name = self.pdb_ccd_rdkit.chem_comp_name
        svg_diagram = self.pdb_ccd_rdkit.image_file_or_string(atom_labels=True, pixels_x=PIXELS_X, pixels_y=PIXELS_Y)
        title = 'proof of concept - Mogul analysis of PDB-CCD coordinates for {}'.format(chem_comp_id)
        bond_title, bond_rows, bond_svg = self.prepare_bond_table()
        logging.debug(bond_title)

        with tag('html'):
            with tag('head'):
                with tag('title'):
                    text(title)
                with tag('style'):
                    text('table, th, td {border: 2px solid black; border-collapse: collapse;}')
                    text('th, td { padding: 5px; text-align: center }')
                doc.asis(JQUERY_SCRIPT)
            with tag('body'):
                with tag('h1'):
                    text(title)
                with tag('ul', ):
                    line('li', 'chem_comp_id =' + chem_comp_id)
                    line('li', 'chem_comp_name = ' + chem_comp_name)
                doc.asis(svg_diagram)
                with tag('h2'):
                    text('bond lengths')
                if len(self.store_bonds) == 0:
                    line('p', 'no bonds found')
                else:
                    doc.asis(bond_svg)
                    with tag('div', id="bond_show_button"):
                        with tag('button', klass='toggle', value='bond'):
                            text('Show detailed table showing results for each bond')
                    with tag('div', id="bond_details"):
                        with tag('button', klass='toggle', value='bond'):
                            text('Hide hide detailed table')
                        with tag('table'):
                            with tag('tr'):
                                for item in bond_title:
                                    with tag('th'):
                                        doc.asis(item)
                            for row in bond_rows:
                                with tag('tr'):
                                    for item in row:
                                        with tag('td'):
                                            text(item)
                        with tag('button', klass='toggle', value='bond'):
                            text('Hide hide detailed table')
        result = doc.getvalue()
        return result

    def prepare_bond_table(self):
        title_row = ('atoms', 'actual in ' + ANGSTROM, 'Mogul mean in ' + ANGSTROM, 'difference in ' + ANGSTROM,
                     'Mogul ' + SIGMA + ' in ' + ANGSTROM, ' Mogul # hits', 'Z*-score', 'classification')
        rows = []
        for bond in sorted(self.classify_bonds, key=lambda b: b.zorder, reverse=True):
            atoms = '-'.join(bond.atoms_ids)
            actual = '{:.3f}'.format(bond.value)
            mean = '{:.3f}'.format(bond.mean)
            difference = '{:.3f}'.format(bond.value - bond.mean)
            sigma = '{:.3f}'.format(bond.standard_deviation)
            nhits = '{}'.format(bond.nhits)
            try:
                z_score = '{:.2f}'.format(bond.zstar)
            except ValueError:
                z_score = ' '
            classification = bond.classification
            rows.append((atoms, actual, mean, difference, sigma, nhits, z_score, classification))

        highlight_bonds = collections.OrderedDict()
        for bond in sorted(self.classify_bonds, key=lambda b: b.zorder):
            classification = bond.classification
            if classification == 'too few hits':
                pass
            else:
                highlight_bonds[(min(bond.indices), max(bond.indices))] =  CLASSIFICATION_COLOR[classification]
        svg_string = self.pdb_ccd_rdkit.image_file_or_string( hydrogen=False, atom_labels=False, wedge=False,
                                                              highlight_bonds=highlight_bonds, black=True,
                                                              pixels_x=PIXELS_X, pixels_y=PIXELS_Y)
        return title_row, rows, svg_string

