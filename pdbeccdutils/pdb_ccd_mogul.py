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
from math import sqrt
from ccdc import conformer
from ccdc.molecule import Molecule as ccdcMolecule
from ccdc.descriptors import MolecularDescriptors as MD
from pdbeccdutils.pdb_chemical_components_rdkit import PdbChemicalComponentsRDKit
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
                        0: None}
CLASSIFICATION_HTML_COLOR = {5: '#D73027',
                             4: '#FC8D59',
                             3: '#FEE090',
                             2: '#91BFDB',
                             1: '#4575B4',
                             0: '#FFFFFF'}
STYLE = '''
<style>
table, th, td {border: 2px solid black; border-collapse: collapse;}
th, td { padding: 5px; text-align: center }
table.no_border, th.no_border, td.no_border { border: 0px ; border-collapse: collapse;}
th.key, td.key { border: 4px solid white; text-align: left}
</style>
'''
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
        self.svg_coloured_diagram_labels = {}
        """for each observation type"""

    def run_mogul(self):
        """
        Runs the Mogul analysis, storing results in self.store_* lists

        Returns:
            None
        """
        # write out the molecule as an sdf file to a temporary directory
        override_xyz = None
        if self.pdb_ccd_rdkit.ideal_xyz[0][2] == 0.000:
            logging.debug('ideal xyz first z coordinate is 0.000')
            override_xyz = []
            for xyz in self.pdb_ccd_rdkit.ideal_xyz:
                override_xyz.append((xyz[0]+0.001, xyz[1]+0.001, xyz[2]+0.001))
            logging.debug(override_xyz)
        sdf_string = self.pdb_ccd_rdkit.sdf_file_or_string(xyz=override_xyz)
        logging.debug('load ccd into PdbChemicalComponentsRDKit sdf string=\n{}'.format(sdf_string))
        molecule = ccdcMolecule.from_string(sdf_string, format='sdf')
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
            if observation_type == 'ring':
                store['sybyl_atom_types'] = self.query_ring_sybyl_types(geometry_analysed_molecule, this_ring=thing)
                store['ring_torsions_labels'] = self.query_ring_tors_labels(atom_ids)
                store['ring_hits'] = self._ring_hit_list(this_ring=thing)
            for key in ['classification', 'd_min', 'local_density', 'lower_quartile', 'maximum',
                        'mean', 'median', 'minimum', 'nhits', 'standard_deviation', 'type',
                        'unusual', 'upper_quartile', 'value', 'z_score']:
                store[key] = getattr(thing, key)
            store['histogram'] = thing.histogram(minimum=0.0, maximum=hist_max)
            store['hist_max'] = hist_max
            store_nt = collections.namedtuple('stored_mogul_' + observation_type, store.keys())(**store)
            place_in.append(store_nt)
            logging.debug('store {}:'.format(observation_type))
            for name, value in store.items():
                if name == 'ring_hits':
                    logging.debug('\t\tring_hits:')
                    for hit in value:
                        logging.debug('\t\t\t{}'.format(hit))
                else:
                    logging.debug('\t\t{}\t{}'.format(name, value))

    @staticmethod
    def query_ring_sybyl_types(geometry_analysed_molecule, this_ring):
        number_atoms_in_ring = len(this_ring.atom_indices)
        query_atoms = geometry_analysed_molecule.atoms
        query_ring_sybyl_types = []
        for i0 in range(number_atoms_in_ring):
            csd_atom = query_atoms[this_ring.atom_indices[i0]]
            query_ring_sybyl_types.append(str(csd_atom.sybyl_type))
        return query_ring_sybyl_types

    @staticmethod
    def query_ring_tors_labels(atom_ids):
        number_atoms_in_ring = len(atom_ids)
        query_ring_torsions_labels = []
        for i0 in range(number_atoms_in_ring):
            i1 = (i0 + 1) % number_atoms_in_ring
            i2 = (i0 + 2) % number_atoms_in_ring
            i3 = (i0 + 3) % number_atoms_in_ring
            tors_label = '{}-{}-{}-{}'.format(atom_ids[i0], atom_ids[i1], atom_ids[i2], atom_ids[i3])
            query_ring_torsions_labels.append(tors_label)
        return query_ring_torsions_labels

    def _ring_hit_list(self, this_ring):
        hit_list = []
        for hit in this_ring.hits:
            this_hit = collections.OrderedDict()
            this_hit['csd_identifier'] = str(hit.identifier)
            this_hit['csd_atom_ids'], this_hit['sybyl_types'], this_hit['ring_torsions'] = \
                self.ring_hit_information_to_store(hit.atoms)
            this_hit['rms_ring_torsion'] = self.rms(this_hit['ring_torsions'])
            hit_list.append(this_hit)
        return hit_list

    @staticmethod
    def ring_hit_information_to_store(hit_atoms):
        """
        process a ring hit - finding the atom_ids, sybyl_atom_types and ring torsions

        Args:
            hit_atoms: list of csd python api atoms

        Returns:
            tuple (list of atom ids, list of  sybyl_atom_types, list of ring torsions angles in degrees)
        """
        number_atoms_in_ring = len(hit_atoms)
        csd_atom_ids = []
        sybyl_types = []
        hit_ring_torsions = []
        for i0 in range(number_atoms_in_ring):
            csd_atom_ids.append(str(hit_atoms[i0].label))
            sybyl_types.append(str(hit_atoms[i0].sybyl_type))
            i1 = (i0 + 1) % number_atoms_in_ring
            i2 = (i0 + 2) % number_atoms_in_ring
            i3 = (i0 + 3) % number_atoms_in_ring
            tors = MD.atom_torsion_angle(hit_atoms[i0], hit_atoms[i1], hit_atoms[i2], hit_atoms[i3])
            hit_ring_torsions.append(tors)
        return csd_atom_ids, sybyl_types, hit_ring_torsions

    @staticmethod
    def rms(float_list):
        """
        works out the root mean squared for a list of floats
        
        Args:
            float_list: a list of floats. 

        Returns:

        """
        this_sum = 0.
        for element in float_list:
            this_sum += element*element
        return sqrt(this_sum/float(len(float_list)))

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
        if observation_type == 'bond' or observation_type == 'angle':
            self.classify_observation_bonds_or_angles(observation_type)
        elif observation_type == 'torsion':
            return  # TODO code up!
        elif observation_type == 'ring':
            self.score_and_classify_rings()
            return
        else:
            raise RuntimeError('unrecognized observation_type={}'.format(observation_type))

    def classify_observation_bonds_or_angles(self, observation_type):
        """
        classifies Mogul results using own analysis.

        Args:
            observation_type (str): one of 'bond', 'angle'

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
        else:
            raise RuntimeError('call to classify_observation_bonds_or_angles with bad observation')
        for thing in work_from:
            if thing.nhits == 0:
                zstar = None
            else:
                zstar = (thing.value - thing.mean)/max(thing.standard_deviation, sd_min)
            if thing.nhits >= few_hits_threshold:
                zorder = abs(zstar)
            elif thing.nhits != 0:
                zorder = -100. + abs(zstar)
            else:
                zorder = -101.
            logging.debug('zstar={} zorder={}'.format(zstar, zorder))
            classification = None
            for class_num, limit in CLASSIFICATION_ZLIMIT.items():
                if zorder > limit:
                    classification = class_num
                    break
            classify = thing._asdict()
            classify['zstar'] = zstar
            classify['zorder'] = zorder
            classify['classification'] = classification
            store_nt = collections.namedtuple('classify_mogul_' + observation_type, classify.keys())(**classify)
            logging.debug(store_nt)
            place_in.append(store_nt)

    def score_and_classify_rings(self):
        """
        Scores ring using RDKit measurements of ring torsion angles and then classifies them into self.classify_rings
        list.
        """
        for store_ring in self.store_rings:
            # find the query ring torsion angles in degrees - using rdkit
            query_ring_torsions = []
            indices = store_ring.indices
            number_atoms_in_ring = len(store_ring.indices)
            for i0 in range(number_atoms_in_ring):
                i1 = (i0 + 1) % number_atoms_in_ring
                i2 = (i0 + 2) % number_atoms_in_ring
                i3 = (i0 + 3) % number_atoms_in_ring
                torsion_indices = (indices[i0], indices[i1], indices[i2], indices[i3])
                torsion = self.pdb_ccd_rdkit.calculate_torsion(atom_indices=torsion_indices)
                query_ring_torsions.append(torsion)
            scored_hits = self.score_ring(query_ring_torsions=query_ring_torsions, store_ring=store_ring)
            logging.debug('scored_hits = {}'.format(scored_hits))
            # TODO classify ring

    def score_ring(self, query_ring_torsions, store_ring):
        """
        'Scores' a ring working out the strangeness (rmsd of query ring torsions to each of the CSD hits stored).

        Args:
            query_ring_torsions: list of ring torsion angles in degrees for the query molecule.
            store_ring: a stored_mogul_ring named tuple object containing stored information to classify ring.

        Returns:
            list of scored hits

        """
        logging.debug('call to score_ring')
        logging.debug('(ideal) query_ring_torsions {}'.format(query_ring_torsions))
        logging.debug('store_ring={}'.format(store_ring))
        scored_hits = []
        for hit in store_ring.ring_hits:
            scored_hits.append(self.score_ring_hit(query_ring_torsions=query_ring_torsions,
                                                   query_sybyl_atom_types=store_ring.sybyl_atom_types,
                                                   hit=hit))
        return scored_hits
    @staticmethod
    def score_ring_hit(query_ring_torsions, query_sybyl_atom_types, hit):
        """
        Find the ring strangeness of the query molecule to a hit from the stored_mogul_ring named tuple rings_hit list
        of information about a hit to a particular ring from a CSD enty.

        Args:
            query_ring_torsions: list of ring torsion angles in degrees for the query molecule.
            query_sybyl_atom_types: list of the sybyl atoms types of the query ring
            hit: an ordered dictionary created by the self._ring_hit_list method

        Returns:
            OrderDict with results for this hit

        Notes:
            to score the ring have to match up the query ring sybyl atom types with the hits. Have to try all
            combinations in case the ring has identical types. Ring has to be reversed (going the other way around)

            working out the ring torsions on reversal could be done in a neater way - currently a dictionary is used.
        """
        hit_sybyl_atom_types = hit['sybyl_types']
        number_atoms_in_ring = len(query_ring_torsions)
        csd_atom_ids = hit['csd_atom_ids']
        hit_ring_tors_dict = {}
        for i0 in range(number_atoms_in_ring):
            i1 = (i0 + 1) % number_atoms_in_ring
            i2 = (i0 + 2) % number_atoms_in_ring
            i3 = (i0 + 3) % number_atoms_in_ring
            hit_ring_tors_label = '{}-{}-{}-{}'.format(csd_atom_ids[i0], csd_atom_ids[i1],
                                                       csd_atom_ids[i2], csd_atom_ids[i3])
            hit_ring_tors_reversed = '{}-{}-{}-{}'.format(csd_atom_ids[i3], csd_atom_ids[i2],
                                                          csd_atom_ids[i1], csd_atom_ids[i0])
            hit_ring_tors_dict[hit_ring_tors_label] = hit['ring_torsions'][i0]
            hit_ring_tors_dict[hit_ring_tors_reversed] = hit['ring_torsions'][i0]
        strangeness = 99999.00
        match_invert = None
        matched_atoms_labels = None
        matched_torsions = None
        for reverse in False, True:
            for offset_by in range(number_atoms_in_ring):
                offset_indices = []
                for ia in range(number_atoms_in_ring):
                    offset_indices.append((ia + offset_by) % number_atoms_in_ring)
                if reverse:
                    offset_indices.reverse()
                sybyl_types_match = True
                for ia in range(number_atoms_in_ring):
                    if query_sybyl_atom_types[ia] != hit_sybyl_atom_types[offset_indices[ia]]:
                        sybyl_types_match = False
                if not sybyl_types_match:
                    continue
                sum_delta_squared = 0.
                sum_delta_squared_invert = 0.
                offset_atoms_labels = []
                hit_ring_torsions = []
                hit_ring_torsions_invert = []
                for ia in range(number_atoms_in_ring):
                    query_ring_tor = query_ring_torsions[ia]
                    i0 = offset_indices[ia]
                    i1 = offset_indices[(ia + 1) % number_atoms_in_ring]
                    i2 = offset_indices[(ia + 2) % number_atoms_in_ring]
                    i3 = offset_indices[(ia + 3) % number_atoms_in_ring]
                    hit_ring_tors_label = '{}-{}-{}-{}'.format(csd_atom_ids[i0], csd_atom_ids[i1],
                                                               csd_atom_ids[i2], csd_atom_ids[i3])
                    tors = hit_ring_tors_dict[hit_ring_tors_label]
                    offset_atoms_labels.append(csd_atom_ids[i0])
                    delta = query_ring_tor - tors
                    delta_invert = query_ring_tor + tors
                    sum_delta_squared += delta * delta
                    sum_delta_squared_invert += delta_invert * delta_invert
                    hit_ring_torsions.append(tors)
                    hit_ring_torsions_invert.append(-tors)
                my_ring_rmsd = sqrt(sum_delta_squared / float(number_atoms_in_ring))
                my_ring_rmsd_invert = sqrt(sum_delta_squared_invert / float(number_atoms_in_ring))
                if my_ring_rmsd < strangeness:
                    strangeness = my_ring_rmsd
                    match_invert = False
                    matched_atoms_labels = offset_atoms_labels
                    matched_torsions = hit_ring_torsions
                if my_ring_rmsd_invert < strangeness:
                    strangeness = my_ring_rmsd_invert
                    match_invert = True
                    matched_atoms_labels = offset_atoms_labels
                    matched_torsions = hit_ring_torsions_invert
        if strangeness == 99999.00:
            raise RuntimeError('impossible error failed to find any match to hit {}'.format(hit))
        logging.debug('{} strangeness={} {} invert={} {}'.
                      format(hit['csd_identifier'], strangeness, matched_atoms_labels, match_invert, matched_torsions))
        hit_score = collections.OrderedDict()
        hit_score['csd_identifier'] = hit['csd_identifier']
        hit_score['strangeness'] = strangeness
        hit_score['matched_atoms_labels'] = matched_atoms_labels
        hit_score['match_invert'] = match_invert
        hit_score['match_invert'] = matched_torsions
        return hit_score

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
        if len(work_from) > 0:
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
            logging.debug('thing={}'.format(thing))
            classification = thing.classification
            if classification == 0:  # too few hits
                pass
            elif observation_type == 'bond':
                highlight_bonds[(thing.atoms_ids[0], thing.atoms_ids[1])] = CLASSIFICATION_COLOR[classification]
            elif observation_type == 'angle':
                logging.debug('thing.indices is {} thing.atom_ids={} classification={}'.
                              format(thing.indices, thing.atoms_ids, classification))
                first_bond_in_angle = (thing.atoms_ids[0], thing.atoms_ids[1])
                second_bond_in_angle = (thing.atoms_ids[1], thing.atoms_ids[2])
                logging.debug(' first_bond_in_angle is {}'.format(first_bond_in_angle))
                logging.debug(' second_bond_in_angle is {}'.format(second_bond_in_angle))
                for this_bond in first_bond_in_angle, second_bond_in_angle:
                    if this_bond in highlight_bonds:
                        del highlight_bonds[this_bond]
                    highlight_bonds[this_bond] = CLASSIFICATION_COLOR[classification]
        logging.debug('highlight_bonds={}'.format(highlight_bonds))
        for atom_labels in False, True:
            svg_string = self.pdb_ccd_rdkit.image_file_or_string(hydrogen=False, atom_labels=atom_labels, wedge=False,
                                                                 highlight_bonds=highlight_bonds, black=True,
                                                                 pixels_x=PIXELS_X, pixels_y=PIXELS_Y)
            svg_string = svg_string.replace('svg:', '')
            if not atom_labels:
                self.svg_coloured_diagram[observation_type] = svg_string
            else:
                self.svg_coloured_diagram_labels[observation_type] = svg_string

    def prepare_file_html(self, html_file_name):
        html_text = self.prepare_html()
        with open(html_file_name, "w") as text_file:
            text_file.write(html_text)

    def prepare_html(self):
        doc, tag, text, line = Doc().ttl()

        chem_comp_id = self.pdb_ccd_rdkit.chem_comp_id
        chem_comp_name = self.pdb_ccd_rdkit.chem_comp_name
        title = 'proof of concept - Mogul analysis of PDB-CCD coordinates for {}'.format(chem_comp_id)
        with tag('html'):
            with tag('head'):
                with tag('title'):
                    text(title)
                doc.asis(STYLE)
                doc.asis(JQUERY_SCRIPT)
            with tag('body'):
                with tag('h1'):
                    text(title)
                with tag('ul', ):
                    with tag('li'):
                        text('code: ')
                        with tag('a', href='https://gitlab.com/pdbe/ccd_utils/'):
                            text('https://gitlab.com/pdbe/ccd_utils/')
                    with tag('li'):
                        text('development notes: ')
                        with tag('a', href='https://gitlab.com/pdbe/ccd_utils/issues/20'):
                            text('https://gitlab.com/pdbe/ccd_utils/issues/20')
                    line('li', 'chem_comp_id ' + chem_comp_id)
                    line('li', 'chem_comp_name ' + chem_comp_name)
                    line('li', "This analysis is of the wwPDB chemical component definition 'ideal' coordinates.")
                for observation_type in MOGUL_OBSERVATION_TYPES:
                    self.prepare_html_section(observation_type, doc, tag, text, line)
                with tag('h3'):
                    text('Mogul run conditions')
                with tag('i'):
                    text('TODO: add Mogul run conditions')

        result = doc.getvalue()
        return result

    def prepare_html_section(self, observation_type, doc, tag, text, line):
        if observation_type == 'bond':
            section_title = 'bond lengths'
        elif observation_type == 'angle':
            section_title = 'bond angles'
        elif observation_type == 'torsion':
            section_title = 'torsions'
        elif observation_type == 'ring':
            section_title = 'rings'
        else:
            raise RuntimeError('unrecognized observation_type={}'.format(observation_type))
        with tag('h3'):
            text(section_title)
        if observation_type not in self.detailed_html_table:
            line('p', 'no ' + observation_type + 's found or NOT YET CODED!')
        else:
            with tag('div', id=observation_type + '_show_button'):
                # svg and key in little table
                with tag('table', klass='no_border'):
                    with tag('tr', klass='no_border'):
                        with tag('td', klass='no_border'):
                            doc.asis(self.svg_coloured_diagram[observation_type])
                        with tag('td', klass='no_border'):
                            self.bond_angle_key(doc, tag, text)
                with tag('button', klass='toggle', value=observation_type):
                    text('Show {} details'.format(observation_type))
            with tag('div', id=observation_type + "_details"):
                with tag('table', klass='no_border'):
                    with tag('tr', klass='no_border'):
                        with tag('td', klass='no_border'):
                            doc.asis(self.svg_coloured_diagram_labels[observation_type])
                        with tag('td', klass='no_border'):
                            self.bond_angle_key(doc, tag, text)
                with tag('button', klass='toggle', value=observation_type):
                    text('Hide {} details'.format(observation_type))
                with tag('p'):
                    with tag('i'):
                        text('TODO: add table with metrics - number and % for outliers, very-unusual, usual plus rmsZ*')
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
                if len(self.detailed_html_table[observation_type]) > 10:
                    with tag('button', klass='toggle', value=observation_type):
                        text('Hide {} details'.format(observation_type))

    @staticmethod
    def bond_angle_key( doc, tag, text):
        z_next = ''
        with tag('table'):
            for class_num, name in CLASSIFICATION_NAME.items():
                if class_num == 0:
                    break
                with tag('tr',  klass='key'):
                    with tag('td', bgcolor=CLASSIFICATION_HTML_COLOR[class_num], klass='key'):
                        text(name)
                    with tag('td', klass='key'):
                        z_limit = '{} &gt; Z* {}'.format(CLASSIFICATION_ZLIMIT[class_num], z_next)
                        z_next = '&ge; {}'.format(CLASSIFICATION_ZLIMIT[class_num])
                        doc.asis(z_limit)
