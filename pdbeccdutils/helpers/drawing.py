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

"""Module helper providing drawinf functionality based on RDKit
"""

import os
import re
from sys import platform

import rdkit
import xml.etree.ElementTree as ET
from PIL import Image, ImageDraw, ImageFont


def supply_font():
    """
    Platform non-specific function to locate sans-serif font in the environment.

    Returns:
        str: path to the font
    """
    font = ''
    if platform == "linux" or platform == "linux2":
        font = '/usr/share/fonts/gnu-free/FreeSans.ttf'
    elif platform == "darwin":
        font = '/Library/Fonts/arial.ttf'
    elif platform == "win32":
        font = 'c:\\windows\\font\\arial.ttf'

    if os.path.isfile(font):
        return font
    else:
        return None


def save_no_image(path_to_image, width=200):
    """
    Generate pretty image with 'No image available' message in case
    the 2D depiction cannot be created.

    Args:
        path_to_image (str): path to the image
        width (int, optional): Defaults to 200. width of the image
    """
    if path_to_image.split('.')[-1] == "svg":
        svg = _svg_no_image(width)
        with open(path_to_image, 'w') as f:
            f.write(svg)
    else:
        _png_no_image(path_to_image, width)


def _png_no_image(path_to_image, width):
    """
    Save image with the text 'No image available' as a png.

    Args:
        path_to_image (str): path to save the image
        width (int): Width of an image
    """

    font = None
    font_path = supply_font()

    if font is not None:
        font_path = ImageFont.truetype(font_path, size=(int(width / 8)))
    else:
        font = ImageFont.load_default()

    white = (255, 255, 255)
    black = (0, 0, 0)
    img = Image.new("RGBA", (width, width), white)
    draw = ImageDraw.Draw(img)
    draw.multiline_text((width / 4, width / 3), "No image\n available",
                        font=font, align='center', fill=black)
    draw = ImageDraw.Draw(img)
    img.save(path_to_image)


def _svg_no_image(width=200):
    """Get svg representation
        width (int, optional): Defaults to 200. width of the image

    Returns:
        str: string representation of an svg image.
    """

    svg = """<?xml version='1.0' encoding='iso-8859-1'?>
            <svg version='1.1' baseProfile='full'
              xmlns='http://www.w3.org/2000/svg'
                      xmlns:rdkit='http://www.rdkit.org/xml'
                      xmlns:xlink='http://www.w3.org/1999/xlink'
                  xml:space='preserve' width='{width}px' height='{width}px' >
            <rect style='opacity:1.0;fill:#FFFFFF;stroke:none' width='{width}' height='{width}' x='0' y='0'> </rect>
            <text alignment-baseline="middle" text-anchor="middle" x="25%" y="25%" style='font-size:{font}px;font-family:sans-serif;text-anchor:start;fill:#000000'>
                Image not
            </text>
            <text alignment-baseline="middle" text-anchor="middle" x="30%" y="50%" style='font-size:{font}px;font-family:sans-serif;text-anchor:start;fill:#000000'>
                available
            </text>
            </svg>
          """
    return svg.format(width=width, font=width / 8)


def draw_molecule(mol, drawer, file_name, width, atom_highlight, bond_highlight):
    """Draw SVG image from the RDKit molecule.

    Args:
        mol (rdkit.Chem.rdchem.Mol): Rdkit mol object to be depicted.
        drawer ([type]): [description]
        file_name (str): Path where the depiction will be saved.
        width (int): Size of the image. Final image will be in NxN
            resolution.
        atom_highlight (`obj`: Dict): Dictionary with atom id and RGB
            mapping color mapping.
        bond_highlight (`obj`: Dict): Dictionary with mapping of atom
            ids and RGB colors.
    """
    try:
        copy = rdkit.Chem.Draw.rdMolDraw2D.PrepareMolForDrawing(mol, wedgeBonds=True,
                                                                kekulize=True, addChiralHs=True)
    except (RuntimeError, ValueError):
        copy = rdkit.Chem.Draw.rdMolDraw2D.PrepareMolForDrawing(mol, wedgeBonds=False,
                                                                kekulize=True, addChiralHs=True)

    if bond_highlight is None:
        drawer.DrawMolecule(copy, highlightAtoms=atom_highlight.keys(),
                            highlightAtomColors=atom_highlight)
    else:
        drawer.DrawMolecule(copy, highlightAtoms=atom_highlight.keys(),
                            highlightAtomColors=atom_highlight,
                            highlightBonds=bond_highlight.keys(), highlightBondColors=bond_highlight)
    drawer.FinishDrawing()

    with open(file_name, 'w') as f:
        svg = drawer.GetDrawingText()

        if width < 201:
            svg = re.sub('stroke-width:2px', 'stroke-width:1px', svg)
        f.write(svg)


def get_drawing_scale(mol):
    """Calculate molecule resolution given 50points per Angstroom

    Args:
        mol (rdkit.Chem.rdchem.Mol): Rdkit mol object.

    Returns:
        [:obj:`tuple` of :obj:`int`]: Dimension of the depictions.
    """
    cnf = mol.GetConformer()

    a = rdkit.Geometry.Point2D(0, 0)
    b = rdkit.Geometry.Point2D(0, 0)
    a.x = b.x = cnf.GetAtomPosition(0).x
    a.y = b.y = cnf.GetAtomPosition(0).y

    for i in range(1, cnf.GetNumAtoms()):
        a.x = min(a.x, cnf.GetAtomPosition(i).x)
        a.y = min(a.y, cnf.GetAtomPosition(i).y)
        b.x = max(b.x, cnf.GetAtomPosition(i).x)
        b.y = max(b.y, cnf.GetAtomPosition(i).y)

    w = 50 * (b.x - a.x) + 1
    h = 50 * (b.y - a.y) + 1

    return (int(w), int(h))


def parse_svg(svg_string, mol: rdkit.Chem.Mol):
    """Parse information from SVG depiction into object.

    Args:
        svg_string (str): SVG as string.
        mol (rdkit.Chem.Mol): RDKit mol object used for depiction.

    Returns:
        :obj:`dict` of :obj:`dict`: object with all the details for
        json serialization.
    """
    svg_string = _fix_svg(svg_string)
    depiction = {}

    svg = ET.fromstring(svg_string)
    atoms_svg = svg.findall('{http://www.rdkit.org/xml}atom')
    bonds_svg = svg.findall('{http://www.w3.org/2000/svg}path')
    dimensions_svg = svg.find('{http://www.w3.org/2000/svg}rect')
    labels_svg = svg.findall('{http://www.w3.org/2000/svg}text')

    depiction['atoms'] = list(map(lambda atom_svg:
                                  {'name': mol.GetAtomWithIdx(int(atom_svg.attrib.get('idx')) - 1).GetProp('name'),
                                   'x': atom_svg.attrib.get('x'),
                                   'y': atom_svg.attrib.get('y')
                                   }, atoms_svg))

    depiction['bonds'] = list(map(lambda bond_svg:
                                  {'coords': bond_svg.attrib.get('d'),
                                   'style': bond_svg.attrib.get('style')
                                   }, bonds_svg))

    depiction['labels'] = list(map(lambda label_svg:
                                   {
                                       'x': label_svg.attrib.get('x'),
                                       'y': label_svg.attrib.get('y'),
                                       'style': label_svg.attrib.get('style'),
                                       'tspans': [{
                                           'value': tspan.text,
                                           'style': '' if tspan.attrib.get('style') is None else tspan.attrib.get('style')
                                       }
                                           for tspan in filter(lambda x: x.text is not None, label_svg.findall('{http://www.w3.org/2000/svg}tspan'))]
                                   }, labels_svg))

    depiction['dimensions'] = {
        'x': dimensions_svg.attrib.get('width'),
        'y': dimensions_svg.attrib.get('height')
    }

    return depiction


def _fix_svg(svg_string):
    
    svg_string = re.sub('<sub>', '&lt;sub&gt;', svg_string)
    svg_string = re.sub('</sub>', '&lt;/sub&gt;', svg_string)
    svg_string = re.sub('<sup>', '&lt;sup&gt;', svg_string)
    svg_string = re.sub('</sup>', '&lt;/sup&gt;', svg_string)
    
    return svg_string
