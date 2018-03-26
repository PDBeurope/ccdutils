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
import logging
import os
import sys

from pdbeccdutils.core import ConformerType, structure_reader
from pdbeccdutils.utils import FlatteningManager


def setup_debug(name, log_file, level=logging.INFO):
    handler = logging.FileHandler(log_file)
    logger = logging.getLogger(name)
    logger.setLevel(level)
    logger.addHandler(handler)

    formatter = logging.Formatter('[%(asctime)s]  %(message)s')
    handler.setFormatter(formatter)

    return logger


def _parse_arguments():
    parser = argparse.ArgumentParser(description="Sample script for 2d depiction generation from CCD files.")
    parser.add_argument('-component_file', type=str, help='mmcif component file with all ligands.', required=False)
    parser.add_argument('-pubchem_templates', type=str, help='Path to the pubchem templates.', required=False)
    parser.add_argument('-generic_templates', type=str, help='Path to generic templates.', required=False)
    parser.add_argument('-sdf', type=str, help='Directory to store ideal sdf files.', required=False)
    parser.add_argument('-pdb', type=str, help='Directory to store ideal pdb files.', required=False)
    parser.add_argument('-svg', type=str, help='Directory to store ligand depictions.', required=False)
    parser.add_argument('-logging', action='store_true', help='Turn on logging information.', required=False, default=False)
    parser.add_argument('-with_names', action='store_true', help='Generate figures with atom names too.', required=False, default=False)

    return parser.parse_args()


def check_arguments(config):
    config.component_file = config.component_file or ''

    if not os.path.isfile(config.component_file):
        print('Component file: {} does not exist'.format(config.component_file), file=sys.stderr)
        sys.exit()

    print('Starting mock application of pdbeccdutils...\n')
    _check_path_param(config.svg, 'depictions')
    _check_path_param(config.sdf, 'sdf files')
    _check_path_param(config.pdb, 'pdb files')

    os.makedirs(os.path.join(config.svg, 'with_names'), exist_ok=True)
    os.makedirs(os.path.join(config.svg, 'without_names'), exist_ok=True)


def _check_path_param(param, label):
    if param is None:
        print('INFO | No {} are going to be generated. Argument missing.'.format(label))
    else:
        print('INFO | {} will be saved here: {}'.format(label, param))


def _retrieve_log_file(create_log_file):
    if create_log_file:
        return '2d_depiction.log'

    return os.devnull


def _check_structure(reader_result, log):
    component = reader_result.component

    if len(reader_result.warnings) > 0:
        log.info('{} | warning for structure parsing: {}'
                 .format(component.id, ';'.join(reader_result.warnings)))

    if len(reader_result.errors) > 0:
        log.info('{} | errors for structure parsing: {}'
                 .format(component.id, ';'.join(reader_result.warnings)))

    if component.is_degenerated_conformer(ConformerType.Model):
        log.info('{} | contains degenerated model coordinates.'.format(component.id))

    if component.is_degenerated_conformer(ConformerType.Ideal):
        log.info('{} | contains degenerated ideal coordinates.'.format(component.id))

    sanitization_success = component.sanitize(fast=False)

    if not sanitization_success:
        log.info('{} | sanitization failed!'.format(component.id))


def _generate_depiction(component, flattening, config, log):
    if not os.path.isdir(config.svg):
        return

    depiction_result = component.compute_2d(flattening)
    if depiction_result is not None:
        log.info('{} | 2D depiction generated using {} method with {} template and score {}.'
                 .format(component.id, str(depiction_result.source).split('.')[1],
                         depiction_result.template_name, depiction_result.score))

    file_name = os.path.join(config.svg, 'without_names', component.id + '.svg')
    component.export_2d_svg(file_name, width=500)

    if config.with_names:
        file_name = os.path.join(config.svg, 'with_names', component.id + '.svg')
        component.export_2d_svg(file_name, width=500)


def _save_structures(component, config, log):
    ideal_ok = component.is_degenerated_conformer(ConformerType.Ideal)
    conformer_to_save = ConformerType.Ideal

    if not ideal_ok:
        computed_3d = component.compute_3d()
        if computed_3d:
            conformer_to_save = ConformerType.Computed
        else:
            log.info('{} | cannot generate RDKit conformer. Nothing saved'.format(component.id))
            return

    if os.path.isdir(config.sdf):
        component.export_mol_representation(str_format='sdf', conf_type=conformer_to_save)
        log.info('{} | sdf file with ideal coords generated.'.format(component.id))

    if os.path.isdir(config.pdb):
        component.export_mol_representation(str_format='sdf', conf_type=conformer_to_save)
        log.info('{} | pdb file with ideal coords generated.'.format(component.id))


def main():
    config = _parse_arguments()
    check_arguments(config)
    log_file = _retrieve_log_file(config.logging)
    log = setup_debug('2d_depiction', log_file)

    print('Reading in components.cif file...')
    mmcifs = structure_reader.read_pdb_components_file(config.component_file)
    flattening = FlatteningManager(config.generic_templates, config.pubchem_templates)

    for ccd_id, reader_result in mmcifs.items():
        try:
            component = reader_result.component

            _check_structure(reader_result, log)
            _generate_depiction(component, flattening, config, log)
            _save_structures(component, config, log)
        except Exception:
            log.error('{} FATALITY :)'.format(ccd_id))

    print('We are done in here.')


if __name__ == '__main__':
    main()
