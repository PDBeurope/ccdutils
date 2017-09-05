#!/usr/bin/env python
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
"""
Script for PDBeChem backend infrastructure.
Processes the wwPDB Chemical Components Dictionary file components.cif
producing files for

http://ftp.ebi.ac.uk/pub/databases/msd/pdbechem/

To do this components.cif is split into individual PDB chemical component
definitions cif files, sdf files, pdb files and image files.
In addition creates chem.xml and chem_comp.list for all components.
"""
import argparse
import logging
import os
from argparse import RawTextHelpFormatter
from collections import OrderedDict

import cairosvg

from PIL import Image
from yattag import Doc, indent

from split_components_cif import SplitComponentsCif
from utilities import create_directory_using_mkdir_unless_it_exists

clean_existing = True  # might want an update run mode later but for now remove existing directories/files
file_subdirs = 'mmcif', 'sdf', 'sdf_nh', 'sdf_r', 'sdf_r_nh', 'pdb', 'pdb_r', 'cml', 'xyz', 'xyz_r'
images_subdirs = 'large', 'small', 'hydrogen'


def create_parser():
    """
    Sets up parse the command line options.

    Returns:
         ArgumentParser parser
    """
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=RawTextHelpFormatter)
    parser.add_argument('COMPONENTS_CIF', help='input PDB-CCD components.cif file (must be specified)')
    parser.add_argument('OUTPUT_DIR', help='directory for output (must be specified)')
    parser.add_argument('--debug', action='store_true', help='turn on debug message logging output')
    return parser


def process_components_cif(components_cif, output_dir, debug):
    """
    processes components.cif for PDBeChem type output

    Args:
        components_cif (str): file name/path for components.cif or a test version
        output_dir (str): path for the directory where output will be written
        debug (bool): produce debug type logging

    Returns:

    """
    logger = logging.getLogger(' ')
    if debug:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.WARNING)
    logger.debug('components_cif={} output_dir={}'.format(components_cif, output_dir))
    create_directory_using_mkdir_unless_it_exists(output_dir)
    chem_comp_dot_list_file_name = os.path.join(output_dir, 'chem_comp.list')
    chem_dot_xml_file_name = os.path.join(output_dir, 'chem.xml')
    with open(chem_comp_dot_list_file_name, 'w') as chem_comp_dot_list_file, \
            open(chem_dot_xml_file_name, 'w') as chem_dot_xml_file:
        chem_dot_xml_file.write('<chemCompList>\n')
        files_subdirs_path = _create_files_or_images_subdirs(logger, output_dir, 'files', file_subdirs)
        images_subdirs_path = _create_files_or_images_subdirs(logger, output_dir, 'images', images_subdirs)
        split_cc = SplitComponentsCif(components_cif)
        logger.debug('have opened {} and it contains {} individual CCD cif definitions '.
                     format(components_cif, len(split_cc.cif_dictionary)))
        for pdb_cc_rdkit in split_cc.individual_pdb_ccd_rdkit():
            chem_comp_id = pdb_cc_rdkit.chem_comp_id
            logger.debug('chem_comp_id={}'.format(chem_comp_id))
            chem_comp_dot_list_file.write('{}\n'.format(chem_comp_id))
            chem_dot_xml_file.write(pdb_cc_rdkit.chem_comp_xml())
            _write_coordinate_files_for_ccd(logger, files_subdirs_path, pdb_cc_rdkit, chem_comp_id)
            _write_image_files_for_ccd(logger, images_subdirs_path, pdb_cc_rdkit, chem_comp_id)
        chem_dot_xml_file.write('</chemCompList>\n')
    _create_readme_dot_html(logger, output_dir)
    _create_tar_balls(logger, files_subdirs_path)
    _create_tar_balls(logger, images_subdirs_path)


def _create_readme_dot_html(logger, output_dir):
    """
    writes file to become http://ftp.ebi.ac.uk/pub/databases/msd/pdbechem/readme.htm

    Args:
        logger: logging object
        output_dir (str):  the output directory


    Returns:
        None
    """
    doc, tag, text, line = Doc().ttl()
    title = 'wwPDB ligand dictionary resources'
    with tag('html'):
        with tag('head'):
            with tag('title'):
                text(title)
        with tag('body'):
            with tag('h1'):
                text(title)
            with tag('p'):
                text('This area provides various data files and images for the wwPDB ligand dictionary ')
                ccd_url = 'http://www.wwpdb.org/ccd.html'
                with tag('a', href=ccd_url):
                    text(ccd_url)
                doc.stag('br')
                text('and is part of the PDBeChem web application ')
                with tag('a', href='http://pdbe.org/chem'):
                    text('PDBeChem')
            with tag('p'):
                text('It is updated weekly on the same schedule of the PDB release and PDBe database update.')
            with tag('h2'):
                text('Contents:')
            with tag('h3'):
                text('General')
            with tag('ul', ):
                with tag('li'):
                    text('List of 3 letter codes in the ligand dictionary: ')
                    with tag('a', href='chem_comp.list'):
                        text('chem_comp.list')
                with tag('li'):
                    text('XML file with summary information for each ligand: ')
                    with tag('a', href='chem.xml'):
                        text('chem.xml')
            for section in 'Images', 'Files':
                with tag('h3'):
                    text(section)
                descriptions = OrderedDict()
                if section == 'Images':
                    descriptions['images/large'] = 'Large images with atom names but without hydrogen atoms: '
                    descriptions['images/small'] = 'Small images without hydrogen atoms: '
                    descriptions['images/hydrogen'] = 'Large images with atom names and hydrogen atoms: '
                else:
                    descriptions['files/sdf'] = 'Molfile (SDF) with ideal coordinates and hydrogen atoms: '
                    descriptions['files/sdf_r'] = 'Molfile (SDF) with representative coordinates and hydrogen atoms:'
                    descriptions['files/sdf_nh'] = 'Molfile (SDF) with ideal coordinates without hydrogen atoms: '
                    descriptions['files/sdf_r_nh'] = \
                        'Molfile (SDF) with representative coordinates without hydrogen atoms: '
                    descriptions['files/pdb'] = 'PDB with ideal coordinates: '
                    descriptions['files/pdb_r'] = 'PDB with representative coordinates:'
                    descriptions['files/cml'] = 'CML: '
                    descriptions['files/mmcif'] = 'mmcif individual PDB chemical component definitions: '
                with tag('ul', ):
                    for key, value in descriptions.items():
                        with tag('li'):
                            text(value)
                            with tag('tt'):
                                text(key + ' ')
                            with tag('a', href=key):
                                text('FTP')
                            text(' - ')
                            with tag('a', href=key + '.tar.gz'):
                                text('gzipped tar ball')
    html = indent(doc.getvalue())
    readme_dot_html_file_name = os.path.join(output_dir, 'readme.htm')
    with open(readme_dot_html_file_name, 'w') as readme_dot_html_file:
        readme_dot_html_file.write(html)
        logger.debug('Have written {}'.format(readme_dot_html_file_name))


def _create_files_or_images_subdirs(logger, output_dir, files_or_images, subdirs_list):
    """
    creates the 'files' or 'images' directory and the required subdirectories in it

    Args:
        logger: logging object
        output_dir (str):  the output directory
        files_or_images (str): either 'files' or 'images'
        subdirs_list: list of subdirectories to be created.

    Returns:
        dictionary giving the path to each subdir type

    """
    files_or_images_dir = os.path.join(output_dir, files_or_images)
    create_directory_using_mkdir_unless_it_exists(files_or_images_dir, clean_existing)
    subdirs_path = {}
    for subdir in subdirs_list:
        subdirs_path[subdir] = os.path.join(files_or_images_dir, subdir)
        create_directory_using_mkdir_unless_it_exists(subdirs_path[subdir])
        logger.debug('have created {} subdir {}'.format(files_or_images, subdirs_path[subdir]))
    return subdirs_path


def _write_coordinate_files_for_ccd(logger, subdirs_path, pdb_cc_rdkit, chem_comp_id):
    """
    writes the coordinate files for a particular ccd

    Args:
        logger: logging object
        subdirs_path: dictionary giving the path to each subdir type
        pdb_cc_rdkit (PdbChemicalComponentsRDKit): object for ccd to be written
        chem_comp_id (str): the chem comp id aka 3 letter code for the ccd (eg ATP)

    Returns:
        None
    """
    for subdir in file_subdirs:
        if subdir == 'mmcif':
            file_type = '.cif'
        else:
            file_type = '.' + subdir[:3]
        output_file = os.path.join(subdirs_path[subdir], chem_comp_id + file_type)
        if subdir == 'mmcif':
            pdb_cc_rdkit.write_ccd_cif(output_file)
        elif subdir == 'sdf':
            pdb_cc_rdkit.sdf_file_or_string(file_name=output_file, ideal=True, hydrogen=True)
        elif subdir == 'sdf_nh':
            pdb_cc_rdkit.sdf_file_or_string(file_name=output_file, ideal=True, hydrogen=False)
        elif subdir == 'sdf_r':
            pdb_cc_rdkit.sdf_file_or_string(file_name=output_file, ideal=False, hydrogen=True)
        elif subdir == 'sdf_r_nh':
            pdb_cc_rdkit.sdf_file_or_string(file_name=output_file, ideal=False, hydrogen=False)
        elif subdir == 'pdb':
            pdb_cc_rdkit.pdb_file_or_string(file_name=output_file, ideal=True)
        elif subdir == 'pdb_r':
            pdb_cc_rdkit.pdb_file_or_string(file_name=output_file, ideal=False)
        elif subdir == 'cml':
            pdb_cc_rdkit.cml_file_or_string(file_name=output_file)
        elif subdir == 'xyz':
            pdb_cc_rdkit.xyz_file_or_string(file_name=output_file, ideal=True)
        elif subdir == 'xyz_r':
            pdb_cc_rdkit.xyz_file_or_string(file_name=output_file, ideal=False)
        else:
            raise NotImplementedError('unrecognized subdir {}'.format(subdir))
        if os.path.isfile(output_file):
            logger.debug('written file {}'.format(output_file))
        else:
            logger.warn('failed to write {}'.format(output_file))


def _write_image_files_for_ccd(logger, subdirs_path, pdb_cc_rdkit, chem_comp_id):
    """
    writes the image files for a particular ccd

    Args:
        logger: logging object
        subdirs_path: dictionary giving the path to each subdir type
        pdb_cc_rdkit (PdbChemicalComponentsRDKit): object for ccd to be written
        chem_comp_id (str): the chem comp id aka 3 letter code for the ccd (eg ATP)

    Returns:
        None
    """
    for subdir in subdirs_path.keys():
        output_svg = os.path.join(subdirs_path[subdir], chem_comp_id + '.svg')
        output_png = os.path.join(subdirs_path[subdir], chem_comp_id + '.png')
        output_gif = os.path.join(subdirs_path[subdir], chem_comp_id + '.gif')
        if subdir == 'large':
            pdb_cc_rdkit.image_file_or_string(file_name=output_svg, pixels_x=400, pixels_y=400,
                                              wedge=True, atom_labels=True, hydrogen=False)
        elif subdir == 'small':
            pdb_cc_rdkit.image_file_or_string(file_name=output_svg, pixels_x=100, pixels_y=100,
                                              wedge=True, atom_labels=False, hydrogen=False)
        elif subdir == 'hydrogen':
            pdb_cc_rdkit.image_file_or_string(file_name=output_svg, pixels_x=600, pixels_y=600,
                                              wedge=True, atom_labels=True, hydrogen=True)
        else:
            raise NotImplementedError('unrecognized subdir {}'.format(subdir))
        if os.path.isfile(output_svg):
            logger.debug('written file {}'.format(output_svg))
        else:
            logger.warn('failed to write {}'.format(output_svg))
            continue

        try:
            cairosvg.svg2png(file_obj=open(output_svg, "rb"), write_to=output_png)
        except Exception as ex:
            logging.error('cairosvg.svg2png raised exception on file {}'.format(output_svg))
            logging.error('... exception type: {} message: {}'.format(type(ex).__name__, ex))
            import traceback
            print(traceback.format_exc())
        if os.path.isfile(output_png):
            logger.debug('written file {}'.format(output_png))
        else:
            logger.warn('failed to write {}'.format(output_png))
            continue


        img = Image.open(output_png)
        img.save(output_gif)
        if os.path.isfile(output_gif):
            logger.debug('written file {}'.format(output_gif))
        else:
            logger.warn('failed to write {}'.format(output_gif))
            continue


def _create_tar_balls(logger, subdirs_path):
    """
    creates .tar.gz tarballs for each subdirectory in either files or images

    Args:
        logger: logging object
        subdirs_path: dictionary giving the path to each subdir type

    Returns:
        None
    """
    logger.debug('_create_tar_balls called for {}'.format(subdirs_path))
    for subdir, path in subdirs_path.items():
        path_for_tar_ball = os.path.dirname(path)
        tar_ball_file = subdir + '.tar.gz'
        command = 'cd {}; tar czf {} {}'.format(path_for_tar_ball, tar_ball_file, subdir)
        logger.debug('command {}'.format(command))
        os.system(command)
        tar_ball_file = os.path.join(path_for_tar_ball, tar_ball_file)
        if os.path.isfile(tar_ball_file):
            logger.debug('written tarball {}'.format(tar_ball_file))
        else:
            logger.warn('failed to write {}'.format(tar_ball_file))


def main():
    parser = create_parser()
    args = parser.parse_args()
    process_components_cif(args.COMPONENTS_CIF, args.OUTPUT_DIR, args.debug)

if __name__ == "__main__":
    main()
