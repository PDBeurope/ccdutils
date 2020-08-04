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


"""Module to aid generation of 2D depictions and evaluation of their
quality

"""
import math
import os
import sys
from collections import OrderedDict
from typing import Dict

import pdbeccdutils.helpers.collection_ext as ext
import rdkit
from pdbeccdutils.core.models import DepictionResult, DepictionSource
from pdbeccdutils.utils import config
from rdkit import Chem, Geometry
from rdkit.Chem import AllChem, rdCoordGen
from scipy.spatial import KDTree
from pdbeccdutils.helpers.rdkit_fixtures import fix_conformer


class DepictionManager:
    """
    Toolkit for depicting ligand's structure using RDKit.
    One can supply either templates or 2D depictions by pubchem.
    PubChem templates can be downloaded using PubChemDownloader class.
    """

    def __init__(self, pubchem_templates_path: str = '', general_templates_path: str = config.general_templates) -> None:
        """
        Initialize component which does the ligand depiction.
        If Nones is provided as parameters just the defalt RDKit
        functionality is going to be used.

        Args:
            pubchem_templates_path (str, optional): Defaults to ''.
                Path to the library with 2D structures downloaded from
                PubChem. Use `setup_pubchem_library` for this task.
            general_templates_path (str, optional): Defaults to
                config.general_templates (supplied with the pdbeccdutils).
                Path to the library with general templates to be used
                for depicting ligand e.g. porphyrin rings.

        """
        self.coordgen_params = rdCoordGen.CoordGenParams()
        self.coordgen_params.coordgenScaling = 50 / 1.5
        self.coordgen_params.templateFileDir = config.coordgen_templates

        self.pubchem_templates = pubchem_templates_path if os.path.isdir(pubchem_templates_path) else ''
        self.templates: Dict[str, rdkit.Chem.rdchem.Mol] = OrderedDict()

        if os.path.isdir(general_templates_path):
            for k in sorted(os.listdir(general_templates_path)):
                template = self._load_template(os.path.join(general_templates_path, k))
                template_name = k.split('.')[0]
                self.templates[template_name] = template

    def depict_molecule(self, het_id: str, mol: rdkit.Chem.rdchem.Mol) -> DepictionResult:
        """
        Given input molecule tries to generate its depictions.

        Presently 3 methods are used:
                Pubchem template - find 2d depiction in pubchem db
                User-provided templates - try to use general templates
                From 3D conformer - just apply default RDKit functionality

        Arguments:
            id (str): id of the ligand
            mol (rdkit.Chem.rdchem.Mol): molecule to be depicted

        Returns:
            DepictionResult: Summary of the ligand depiction process.
        """
        temp_mol = Chem.RWMol(mol)

        templateMol = Chem.RWMol(temp_mol).GetMol()
        pubchemMol = Chem.RWMol(temp_mol).GetMol()
        rdkitMol = Chem.RWMol(temp_mol).GetMol()
        results = []

        pubchem_res = self._get_2D_by_pubchem(het_id, pubchemMol) if self.pubchem_templates else None
        template_res = self._get_2D_by_template(templateMol) if self.templates else []
        rdkit_res = self._get_2D_by_rdkit(rdkitMol)

        if pubchem_res is not None:
            results.append(pubchem_res)
        if rdkit_res is not None:
            results.append(rdkit_res)

        results = results + template_res

        results.sort(key=lambda l: (l.score, l.source))

        if results:
            to_return = results[0]
            fix_conformer(to_return.mol.GetConformer(0))
            
            return to_return
            

        return DepictionResult(source=DepictionSource.Failed, template_name='', mol=None, score=1000)

    def _get_pubchem_template_path(self, het_id):
        """Get path to the PubChem template if it exists.

        Args:
            het_id (str): Ligand in.

        Returns:
            str: Path to the PubChem layout.
        """
        path = os.path.join(self.pubchem_templates, f'{het_id}.sdf')

        return path if os.path.isfile(path) else ''

    def _load_template(self, path):
        """
        Loads a template molecule with 2D coordinates

        Args:
            path (str): path to the model molecule in *.sdf,
                or *.pdb format

        Raises:
            ValueError: if unsupported format is used: sdf|pdb

        Returns:
            rdkit.Chem.rdchem.Mol: RDKit representation of the template
        """
        mol = Chem.RWMol()
        extension = os.path.basename(path).split('.')[1]

        if extension == 'sdf':
            mol = Chem.MolFromMolFile(path, sanitize=True, removeHs=True)
        elif extension == 'pdb':
            mol = Chem.MolFromPDBFile(path, sanitize=True, removeHs=True)
        else:
            raise ValueError('Unsupported molecule type \'{}\''.format(extension))

        p = Chem.AdjustQueryParameters()
        p.makeAtomsGeneric = True
        p.makeBondsGeneric = True

        mol = Chem.AdjustQueryProperties(mol, p)

        return mol

    def _get_2D_by_rdkit(self, mol):
        """
        Get depiction done using solely the default RDKit functionality.

        Args:
            mol (rdkit.Chem.rdchem.Mol): Mol to be depicted

        Returns:
            DepictionResult: Depiction with some usefull metadata
        """
        try:
            rdCoordGen.AddCoords(mol, self.coordgen_params)
            flaws = DepictionValidator(mol).depiction_score()
            return DepictionResult(source=DepictionSource.RDKit, template_name=None, mol=mol, score=flaws)
        except Exception:
            return DepictionResult(source=DepictionSource.Failed, template_name=None, mol=None, score=1000)

    def _get_2D_by_pubchem(self, code, mol):
        """
        Depict ligand using Pubchem templates.

        Args:
            code (str): id of the molecule to be processed
            mol (rdkit.Chem.rdchem.Mol): Mol to be depicted

        Returns:
            DepictionResult: Depiction with some usefull metadata
        """
        try:
            template_path = self._get_pubchem_template_path(code)
            if template_path:
                template = self._load_template(template_path)

                if mol.HasSubstructMatch(template):
                    AllChem.GenerateDepictionMatching2DStructure(mol, template)
                    flaws = DepictionValidator(mol).depiction_score()
                    return DepictionResult(source=DepictionSource.PubChem, template_name=code, mol=mol, score=flaws)
        except Exception as e:
            print(str(e), file=sys.stderr)

        return DepictionResult(source=DepictionSource.Failed, template_name=None, mol=None, score=1000)

    def _get_2D_by_template(self, mol):
        """
        Depict ligand using user-provided templates

        Args:
            mol (rdkit.Chem.rchem.Mol): Mol to be depicted

        Returns:
            :obj:`list` of :obj:`DepictionResult`: Depictions with their
            quality and metadata.
        """
        results = list()
        try:
            for key, template in self.templates.items():
                temp_mol = Chem.RWMol(mol)
                if temp_mol.HasSubstructMatch(template):
                    AllChem.GenerateDepictionMatching2DStructure(temp_mol, template)
                    flaws = DepictionValidator(temp_mol).depiction_score()
                    results.append(DepictionResult(source=DepictionSource.Template,
                                                   template_name=key, mol=temp_mol, score=flaws))
        except Exception:
            pass  # if it fails it fails, but generally it wont

        return results


class DepictionValidator:
    """
    Toolkit for estimation of depiction quality
    """

    def __init__(self, mol):
        self.mol = mol
        self.conformer = mol.GetConformer()
        self.bonds = self.mol.GetBonds()

        atoms = [self.conformer.GetAtomPosition(i) for i in range(0, self.conformer.GetNumAtoms())]
        atom_centers = [[atom.x, atom.y, atom.z] for atom in atoms]

        self.kd_tree = KDTree(atom_centers)

    def _intersection(self, bondA, bondB):
        """
        True if two bonds collide, false otherwise. Note that false is
        retrieved even in case the bonds share common atom, as this is
        not a problem case. Cramer's rule is used for the linear
        equations system.

        Args:
            bondA (rdkit.Chem.rdchem.Bond): this bond
            bondB (rdkit.Chem.rdchem.Bond): other bond

        Returns:
            bool: true if bonds share collide, false otherwise.
        """
        atoms = [bondA.GetBeginAtom(), bondA.GetEndAtom(), bondB.GetBeginAtom(), bondB.GetEndAtom()]
        names = [a.GetProp('name') for a in atoms]
        points = [self.conformer.GetAtomPosition(a.GetIdx()) for a in atoms]

        vecA = Geometry.Point2D(points[1].x - points[0].x, points[1].y - points[0].y)
        vecB = Geometry.Point2D(points[3].x - points[2].x, points[3].y - points[2].y)

        # we need to set up directions of the vectors properly in case
        # there is a common atom. So we identify angles correctly
        # e.g. B -> A; B -> C and not A -> B; C -> B.
        if len(set(names)) == 3:
            angle = self.__get_angle(names, vecA, vecB)
            return angle < 10.0

        # Cramer's rule to identify intersection
        det = vecA.x * -vecB.y + vecA.y * vecB.x
        if round(det, 2) == 0.00:
            return False

        a = points[2].x - points[0].x
        b = points[2].y - points[0].y

        detP = (a * -vecB.y) - (b * -vecB.x)
        p = round(detP / det, 3)

        if (p < 0 or p > 1):
            return False

        detR = (vecA.x * b) - (vecA.y * a)
        r = round(detR / det, 3)

        if 0 <= r <= 1:
            return True

        return False

    def __get_angle(self, names, vecA, vecB):
        """Get the size of the angle formed by two bonds which share
        common atom.

        Args:
            names (list of str): List of atom names forming bonds 
                [A, B, C, D] for AB and CD.
            vecA (Geometry.Point2D): Vector representing AB bond.
            vecB (Geometry.Point2D): Vector representing CD bond.

        Returns:
            float: Size of the angle in degrees.
        """
        pivot = ext.find_element_with_max_occurrence(names)

        if names[0] != pivot:  # Atoms needs to be order to pick vectors correctly
            vecA = vecA * -1

        if names[2] != pivot:
            vecB = vecB * -1

        radians = vecA.AngleTo(vecB)
        angle = 180 / math.pi * radians

        return angle

    def has_degenerated_atom_positions(self, threshold):
        """
        Detects whether the structure has a pair or atoms closer to each
        other than threshold. This can detect structures which may need
        a template as they can be handled by RDKit correctly.

        Arguments:
            threshold (float): Bottom line to use for spatial search.

        Returns:
            (bool): if such atomic pair is found
        """

        for i in range(0, len(self.conformer.GetNumAtoms())):
            center = self.conformer.GetAtomPosition(i)
            point = [center.x, center.y, center.z]
            surrounding = self.kd_tree.query_ball_point(point, threshold)

            if len(surrounding) > 1:
                return True

        return False

    def count_suboptimal_atom_positions(self, lowerBound, upperBound):
        """
        Detects whether the structure has a pair or atoms in the range
        <lowerBound, upperBound> meaning that the depiction could
        be improved.

        Arguments:
            lowerBound (float): lower bound
            upperBound (float): upper bound

        Returns:
            bool: indication whether or not the atoms are not in
            optimal coordinates
        """
        counter = 0
        for i in range(self.conformer.GetNumAtoms()):
            center = self.conformer.GetAtomPosition(i)
            point = [center.x, center.y, center.z]
            surroundingLow = self.kd_tree.query_ball_point(point, lowerBound)
            surroundingHigh = self.kd_tree.query_ball_point(point, upperBound)

            if len(surroundingHigh) - len(surroundingLow) > 0:
                counter += 1

        return counter / 2

    def count_bond_collisions(self):
        """
        Counts number of collisions among all bonds. Can be used for estimations of how 'wrong'
        the depiction is.

        Returns:
            int: number of bond collisions per molecule
        """

        errors = 0

        for i in range(0, len(self.bonds)):
            for a in range(i + 1, len(self.bonds)):
                result = self._intersection(self.bonds[i], self.bonds[a])

                if result:
                    errors += 1
        return errors

    def has_bond_crossing(self):
        """
        Tells if the structure contains collisions

        Returns:
            bool: Indication about bond collisions
        """
        return self.count_bond_collisions() > 0

    def depiction_score(self):
        """
        Calculate quality of the ligand depiction. The higher the worse.
        Ideally that should be 0.

        Returns:
            float: Penalty score.
        """

        collision_penalty = 1
        degenerated_penalty = 0.4

        bond_collisions = self.count_bond_collisions()
        degenerated_atoms = self.count_suboptimal_atom_positions(0.0, 0.5)

        score = collision_penalty * bond_collisions + degenerated_penalty * degenerated_atoms

        return round(score, 1)
