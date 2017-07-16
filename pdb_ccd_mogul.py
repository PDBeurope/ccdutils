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
from collections import OrderedDict

ANGSTROM = '&Aring;'
SIGMA = '&sigma;'
CLASSIFICATION_NAME = OrderedDict([(5, 'outlier'),
                                   (4, 'very-unusual'), 
                                   (3, 'unusual'), 
                                   (2, 'common'), 
                                   (1, 'very-common'), 
                                   (0, 'too few hits')])
CLASSIFICATION_ZLIMIT = OrderedDict([(5, 5.0), 
                                     (4, 3.5),
                                     (3, 2.0),
                                     (2, 1.0),
                                     (1, 0.0),
                                     (0, -9999999.0)])
CLASSIFICATION_COLOR = {5: (215. / 255., 48. / 255., 39. / 255.),  # blood orange
                        4: (252. / 255., 141. / 255., 89. / 255.),  # mid orange
                        3: (254. / 255., 224. / 255., 144. / 255.),  # yellow/orange
                        2: (145. / 255., 191. / 255., 219. / 255.),  # mid blue
                        1: (69. / 255., 117. / 255., 180. / 255.),  # blue
                        0: None }
CLASSIFICATION_HTML_COLOR = {5:'#D73027',
                             4:'#FC8D59', 
                             3:'#FEE090',
                             2:'#91BFDB',
                             1:'#4575B4',
                             0:'#FFFFFF'}

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
MOGUL_OBSERVATION_TYPES = ('bond', 'angle', 'torsion', 'ring')

class PdbCCDMogul(object):
    """ run Mogul on PDB CCD"""
    def __init__(self, file_name=None):
        self.settings_bond_few_hits_threshold = 5
        self.settings_angle_few_hits_threshold = 5
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
        self.classify_angles = []
        self.classify_torsions = []
        self.classify_rings = []
        self.detailed_html_table = {}
        """dictionary for each observation, giving detailed html table: 1st row is the header, rest the data"""
        self.svg_coloured_diagram = {}
        """for each observation type"""

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
        engine.settings.angle.few_hits_threshold = self.settings_angle_few_hits_threshold
        self.settings_summary_from_mogul = engine.settings.summary()
        logging.debug('engine.settings.summary()=\n{}'.format(self.settings_summary_from_mogul))
        geometry_analysed_molecule = engine.analyse_molecule(molecule)
        logging.debug('number of Mogul analysed bonds={}'.format(len(geometry_analysed_molecule.analysed_bonds)))
        for observation_type in MOGUL_OBSERVATION_TYPES:
            self.store_observation(geometry_analysed_molecule, observation_type)
        os.remove(sdf_temp)
        for observation_type in MOGUL_OBSERVATION_TYPES:
            self.classify_observation(observation_type)
            self.prepare_html_table(observation_type)
            self.prepare_svg_coloured_diagram(observation_type)

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
        elif observation_type == 'angle':
            work_from = self.store_angles
            place_in = self.classify_angles
            sd_min = 1.0
            few_hits_threshold = self.settings_angle_few_hits_threshold
        elif observation_type == 'torsion':
            return  # TODO code up!
        elif observation_type == 'ring':
            return  # TODO code up!
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
                zorder = -100. + abs(zstar)
            else:
                zorder = -101.
            logging.debug('zstar={} zorder={}'.format(zstar, zorder))
            for class_num, limit in CLASSIFICATION_ZLIMIT.items():
                if zorder > limit:
                    classification = class_num
                    break
            store = thing._asdict()
            store['zstar'] = zstar
            store['zorder'] = zorder
            store['classification'] = classification
            store_nt = collections.namedtuple('stored_mogul_' + observation_type, store.keys())(**store)
            logging.debug(store_nt)
            place_in.append(store_nt)

    def prepare_html_table(self, observation_type):
        if observation_type == 'bond':
            work_from = self.classify_bonds
            units = ANGSTROM
            sf_format = '{:.3f}'
        elif observation_type == 'angle':
            work_from = self.classify_angles
            units = 'degrees'
            sf_format = '{:.1f}'
        elif observation_type == 'torsion':
            return  # TODO code up!
        elif observation_type == 'ring':
            return  # TODO code up!
        else:
            raise RuntimeError('unrecognized observation_type={}'.format(observation_type))
        rows = []
        title_row = ('atoms', 'actual in ' + units, 'Mogul mean in ' + units, 'difference in ' + units,
                     'Mogul ' + SIGMA + ' in ' + units, ' Mogul # hits', 'Z*-score', 'classification', 'color')
        rows.append(title_row)
        for thing in sorted(work_from, key=lambda t: t.zorder, reverse=True):
            atoms = '-'.join(thing.atoms_ids)
            actual = sf_format.format(thing.value)
            mean = sf_format.format(thing.mean)
            difference = sf_format.format(thing.value - thing.mean)
            sigma = sf_format.format(thing.standard_deviation)
            nhits = '{}'.format(thing.nhits)
            try:
                z_score = '{:.2f}'.format(thing.zstar)
            except ValueError:
                z_score = ' '
            classification = CLASSIFICATION_NAME[thing.classification]
            html_color = CLASSIFICATION_HTML_COLOR[thing.classification]
            rows.append((atoms, actual, mean, difference, sigma, nhits, z_score, classification, html_color))
        self.detailed_html_table[observation_type] = rows
 
    def prepare_svg_coloured_diagram(self, observation_type):
        if observation_type == 'bond':
            work_from = self.classify_bonds
        elif observation_type == 'angle':
            work_from = self.classify_angles
        elif observation_type == 'torsion':
            return  # TODO code up!
        elif observation_type == 'ring':
            return  # TODO code up!
        else:
            raise RuntimeError('unrecognized observation_type={}'.format(observation_type))
        highlight_bonds = OrderedDict()
        for thing in sorted(work_from, key=lambda b: b.zorder):
            classification = thing.classification
            if classification == 0:  # too few hits
                pass
            elif observation_type == 'bond':
                highlight_bonds[(min(thing.indices), max(thing.indices))] =  CLASSIFICATION_COLOR[classification]
            elif observation_type == 'angle':
                logging.debug('thing.indices is {} thing.atom_ids={} classification={}'.format(thing.indices, thing.atoms_ids, classification))
                first_bond_in_angle = (min(thing.indices[0:2]), max(thing.indices[0:2]))
                second_bond_in_angle = (min(thing.indices[1:3]), max(thing.indices[1:3]))
                logging.debug(' first_bond_in_angle is {}'.format(first_bond_in_angle))
                logging.debug(' second_bond_in_angle is {}'.format(second_bond_in_angle))
                for this_bond in first_bond_in_angle, second_bond_in_angle:
                    if this_bond in highlight_bonds:
                        del highlight_bonds[this_bond]
                    highlight_bonds[this_bond] =  CLASSIFICATION_COLOR[classification]
        logging.debug('hightlight_bonds={}'.format(highlight_bonds))
        svg_string = self.pdb_ccd_rdkit.image_file_or_string( hydrogen=False, atom_labels=False, wedge=False,
                                                              highlight_bonds=highlight_bonds, black=True,
                                                              pixels_x=PIXELS_X, pixels_y=PIXELS_Y)
        self.svg_coloured_diagram[observation_type] = svg_string

    def prepare_html(self):
        doc, tag, text, line = Doc().ttl()

        chem_comp_id = self.pdb_ccd_rdkit.chem_comp_id
        chem_comp_name = self.pdb_ccd_rdkit.chem_comp_name
        svg_diagram = self.pdb_ccd_rdkit.image_file_or_string(atom_labels=True, pixels_x=PIXELS_X, pixels_y=PIXELS_Y)
        title = 'proof of concept - Mogul analysis of PDB-CCD coordinates for {}'.format(chem_comp_id)
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
                for observation_type in MOGUL_OBSERVATION_TYPES:
                   some_html = self.prepare_html_section( observation_type, doc, tag, text, line)
        result = doc.getvalue()
        return result

    def prepare_html_section( self, observation_type, doc, tag, text, line):
        if observation_type == 'bond':
            section_title = 'bond lengths'
            this_classify = self.classify_bonds
        elif observation_type == 'angle':
            section_title = 'bond angles'
            this_classify = self.classify_angles
        elif observation_type == 'torsion':
            section_title = 'torsions'
            this_classify = self.classify_angles
        elif observation_type == 'ring':
            section_title = 'rings'
            this_classify = self.classify_angles
        else:
            raise RuntimeError('unrecognized observation_type={}'.format(observation_type))
        with tag('h2'):
            text(section_title)
        if observation_type not in self.detailed_html_table:
            line('p', 'no ' + observation_type + 's found')
        else:
            doc.asis(self.svg_coloured_diagram[observation_type])
            with tag('div', id=observation_type + '_show_button'):
                with tag('button', klass='toggle', value=observation_type):
                    text('Show detailed table showing results for each ' + observation_type)
            with tag('div', id=observation_type + "_details"):
                with tag('button', klass='toggle', value=observation_type):
                    text('Hide detailed table')
                with tag('table'):
                    with tag('tr'):
                        for item in self.detailed_html_table[observation_type][0][:-1]:
                            with tag('th'):
                                doc.asis(item)
                    for row in self.detailed_html_table[observation_type][1:]:
                        with tag('tr'):
                            for item in row[:-2]:  # all but the last two
                                with tag('td'):
                                    text(item)
                            with tag('td', bgcolor=row[-1]):
                                text(row[-2])
                with tag('button', klass='toggle', value='bond'):
                    text('Hide detailed table')
 
