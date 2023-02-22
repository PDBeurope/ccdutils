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

"""Module helper providing drawing functionality based on RDKit
"""

import os
import re
import xml.etree.ElementTree as ET
from collections import OrderedDict
from sys import platform

import rdkit
from PIL import Image, ImageDraw, ImageFont

MIN_IMG_DIMENSION = 300
svg_namespace = {"svg": "http://www.w3.org/2000/svg"}


def save_no_image(path_to_image, default_msg=None, width=200):
    """
    Generate pretty image with 'No image available' message in case
    the 2D depiction cannot be created.

    Args:
        path_to_image (str): path to the image
        width (int, optional): Defaults to 200. width of the image
    """
    if path_to_image.split(".")[-1] == "svg":
        svg = (
            _svg_no_image_with_id(default_msg, width)
            if default_msg
            else _svg_no_image(width)
        )
        with open(path_to_image, "w") as f:
            f.write(svg)
    else:
        _png_no_image(path_to_image, width)


def draw_molecule(mol, drawer, file_name, wedge_bonds, atom_highlight, bond_highlight):
    """Draw SVG image from the RDKit molecule.

    Args:
        mol (rdkit.Chem.rdchem.Mol): Rdkit mol object to be depicted.
        drawer (rdkit.Chem.Draw.MolDrawing.DrawingOptions): RDKit object
            with parameters for drawing depiction.
        file_name (str): Path where the depiction will be saved.
        wedge_bonds (bool): Whether or not to wedge bonds.
        atom_highlight (dict): Dictionary with atom id and RGB
            mapping color mapping.
        bond_highlight (dict): Dictionary with mapping of atom
            ids and RGB colors.
    """
    try:
        copy = rdkit.Chem.Draw.rdMolDraw2D.PrepareMolForDrawing(
            mol, wedgeBonds=wedge_bonds, kekulize=True, addChiralHs=True
        )
    except (RuntimeError, ValueError):
        try:
            copy = rdkit.Chem.Draw.rdMolDraw2D.PrepareMolForDrawing(
                mol, wedgeBonds=False, kekulize=True, addChiralHs=True
            )
        except (RuntimeError, ValueError):
            copy = rdkit.Chem.Draw.rdMolDraw2D.PrepareMolForDrawing(
                mol, wedgeBonds=False, kekulize=True, addChiralHs=False
            )

    if bond_highlight is None:
        drawer.DrawMolecule(
            copy,
            highlightAtoms=atom_highlight.keys(),
            highlightAtomColors=atom_highlight,
        )
    else:
        drawer.DrawMolecule(
            copy,
            highlightAtoms=atom_highlight.keys(),
            highlightAtomColors=atom_highlight,
            highlightBonds=bond_highlight.keys(),
            highlightBondColors=bond_highlight,
        )
    drawer.FinishDrawing()

    with open(file_name, "w") as f:
        svg = drawer.GetDrawingText()
        f.write(svg)


def get_drawing_scale(mol):
    """Calculate molecule resolution given 50points per Angstroom

    Args:
        mol (rdkit.Chem.rdchem.Mol): Rdkit mol object.

    Returns:
        [tuple[int,int]]: Dimension of the depictions (x, y).
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

    w = int(50 * (b.x - a.x) + 1)
    h = int(50 * (b.y - a.y) + 1)

    return (max(w, MIN_IMG_DIMENSION), max(h, MIN_IMG_DIMENSION))


def convert_svg(svg_string, ccd_id, mol: rdkit.Chem.Mol):
    """Parse information from SVG depiction into object.

    Args:
        svg_string (str): SVG as string.
        ccd_id (str): CCD ID.
        mol (rdkit.Chem.Mol): RDKit mol object used for depiction.

    Returns:
        dict: JSON-style information for the 2D depiction.
    """

    result_bag = OrderedDict(
        [("ccd_id", ccd_id), ("resolution", {}), ("atoms", []), ("bonds", [])]
    )
    svg = ET.fromstring(svg_string)

    atom_elems = svg.findall("svg:circle", svg_namespace)
    path_elems = svg.findall("svg:path", svg_namespace)
    dimensions_svg = svg.find("svg:rect", svg_namespace)

    atoms = _parse_atoms_from_svg(atom_elems, mol)
    bonds = _parse_bonds_from_svg(path_elems, mol)
    _parse_labels_from_svg(path_elems, atoms)

    result_bag["atoms"] = atoms
    result_bag["bonds"] = bonds

    result_bag["resolution"] = {
        "x": float(dimensions_svg.attrib.get("width")),
        "y": float(dimensions_svg.attrib.get("height")),
    }

    return result_bag


def _parse_atoms_from_svg(atom_elements, mol: rdkit.Chem.Mol):
    """Extract atoms from the SVG atom elements

    Args:
        atom_elements (list[xml.etree.ElementTree.Element]): List of extracted XML elements
        mol (rdkit.Chem.rdchem.Mol): RDkit molecule

    Returns:
        list[dict]: list of JSON-style atom representation.
    """
    result = []
    for atom_svg in atom_elements:
        try:
            atom_id_str = re.search(r"\d+", atom_svg.attrib.get("class")).group(0)
            atom_id = int(atom_id_str)

            if atom_id >= mol.GetNumAtoms():
                continue

            temp = {
                "name": mol.GetAtomWithIdx(atom_id).GetProp("name"),
                "labels": [],
                "x": float(atom_svg.attrib.get("cx")),
                "y": float(atom_svg.attrib.get("cy")),
            }
            result.append(temp)
        except RuntimeError:
            pass  # we do not care for H atoms

    return result


def _parse_labels_from_svg(path_elements, atoms):
    """Parse atom label information from the SVG.

    Args:
        path_elements (list[xml.etree.ElementTree.Element]):
            List of path elements that encode bonds and text.
        atoms (list[dict]): JSON-style representation of atoms.
    """
    atom_id_re = r"atom-\d+"
    for label_svg in path_elements:
        try:
            match = re.fullmatch(atom_id_re, label_svg.attrib["class"])
            if not match:
                continue

            atom_id = int(match.group(0)[5:])
            atoms[atom_id]["labels"].append(
                {"d": label_svg.attrib["d"], "fill": label_svg.attrib["fill"]}
            )
        except (IndexError, KeyError):
            pass  # we do not care for H labels and radicals


def _parse_bonds_from_svg(bond_elements, mol):
    """Extract bonding information from SVG elements

    Args:
        bond_elements (list[xml.etree.ElementTree.Element]):
            List of SVG path elements.
        mol (rdkit.Chem.rdchem.Mol): RDKit mol object]

    Returns:
        list[dict]: JSON-style formated bond informations
    """
    result = []
    re_bond_regex = r"bond-\d+"
    for bond_svg in bond_elements:
        try:
            if not re.search(re_bond_regex, bond_svg.attrib["class"]):
                continue

            atoms = re.findall(r"atom-\d+", bond_svg.attrib["class"])
            atom_id_a = int(atoms[0][5:])  # to get int part of "atom-0"
            atom_id_b = int(atoms[1][5:])

            temp = {
                "bgn": mol.GetAtomWithIdx(atom_id_a).GetProp("name"),
                "end": mol.GetAtomWithIdx(atom_id_b).GetProp("name"),
                "coords": bond_svg.attrib["d"],
                "style": bond_svg.attrib["style"],
            }

            result.append(temp)
        except (RuntimeError, KeyError):
            pass  # we do not care about bonded Hydrogens
    return result


def _png_no_image(path_to_image, width):
    """
    Save image with the text 'No image available' as a png.

    Args:
        path_to_image (str): path to save the image
        width (int): Width of an image
    """

    font = None
    font_path = _supply_font()

    if font is not None:
        font_path = ImageFont.truetype(font_path, size=(int(width / 8)))
    else:
        font = ImageFont.load_default()

    white = (255, 255, 255)
    black = (0, 0, 0)
    img = Image.new("RGBA", (width, width), white)
    draw = ImageDraw.Draw(img)
    draw.multiline_text(
        (width / 4, width / 3),
        "No image\n available",
        font=font,
        align="center",
        fill=black,
    )
    draw = ImageDraw.Draw(img)
    img.save(path_to_image)


def _svg_no_image(width=200):
    """Get svg representation
        width (int, optional): Defaults to 200. width of the image

    Returns:
        str: string representation of an svg image.
    """

    font = width / 8
    svg = f"""<?xml version='1.0' encoding='iso-8859-1'?>
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
    return svg


def _svg_no_image_with_id(name, width=200):
    """Get svg representation
        name (str): Name of the component to be displayed, in case
            the image cannot be generated.
        width (int, optional): Defaults to 200. width of the image

    Returns:
        str: string representation of an svg image.
    """
    font = width / 8
    svg = f"""<?xml version='1.0' encoding='iso-8859-1'?>
            <svg version='1.1' baseProfile='full'
              xmlns='http://www.w3.org/2000/svg'
                      xmlns:rdkit='http://www.rdkit.org/xml'
                      xmlns:xlink='http://www.w3.org/1999/xlink'
                  xml:space='preserve' width='{width}px' height='{width}px' >
            <rect style='opacity:1.0;fill:#FFFFFF;stroke:none' width='{width}' height='{width}' x='0' y='0'> </rect>
            <text  dominant-baseline="middle" text-anchor="middle" x="40%" y="50%" style='font-size:{font}px;font-family:sans-serif;text-anchor:start;fill:#000000'>
                {name}
            </text>
            </svg>
          """
    return svg


def _supply_font():
    """
    Platform non-specific function to locate sans-serif font in the environment.

    Returns:
        str: path to the font
    """
    font = ""
    if platform == "linux" or platform == "linux2":
        font = "/usr/share/fonts/gnu-free/FreeSans.ttf"
    elif platform == "darwin":
        font = "/Library/Fonts/arial.ttf"
    elif platform == "win32":
        font = "c:\\windows\\font\\arial.ttf"

    if os.path.isfile(font):
        return font

    return None
