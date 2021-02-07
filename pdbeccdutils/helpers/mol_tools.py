#!/usr/bin/env python
# software from PDBe: Protein Data Bank in Europe; https://pdbe.org
#
# Copyright 2021 EMBL - European Bioinformatics Institute
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

"""
Set of methods for molecular sanitization and work with conformers
"""

import re
import sys
from io import StringIO

import numpy as np
import rdkit

empty_coords = rdkit.Chem.rdGeometry.Point3D(0, 0, 0)
METALS_SMART = (
    "[Li,Na,K,Rb,Cs,F,Be,Mg,Ca,Sr,Ba,Ra,Sc,Ti,V,Cr,Mn,Fe,Co,Ni,Cu,Zn,Al,Ga,Y,Zr,Nb,Mo,"
    "Tc,Ru,Rh,Pd,Ag,Cd,In,Sn,Hf,Ta,W,Re,Os,Ir,Pt,Au,Hg,Tl,Pb,Bi]"
)


def is_degenerate_conformer(conformer):
    """
    Determine if given conformer has missing coordinates or is
    missing completelly from the rdkit.Mol object. This can
    be used to determine, whether or not the coordinates should be
    regenerated.

    Args:
        conformer (rdkit.Chem.rdchem.Conformer): conformer to check

    Returns:
        bool: true if more then 1 atom has coordinates [0, 0, 0]
    """
    counter = 0

    try:
        for i in range(conformer.GetNumAtoms()):
            pos = conformer.GetAtomPosition(i)
            if pos.Distance(empty_coords) == 0.0:
                counter += 1

            if counter > 1:
                return True
    except ValueError:  # Conformer does not exist
        return False

    return False


def sanitize(rwmol):
    """
    Attempts to sanitize mol in place. RDKit's standard error can be
    processed in order to find out what went wrong with sanitization
    to fix the molecule.

    Args:
        rwmol (rdkit.Chem.rdchem.RWMol): rdkit molecule to be sanitized

    Returns:
        bool: Result of the sanitization process.
    """
    success = False

    try:
        success = fix_molecule(rwmol)

        if not success:
            return False

        rdkit.Chem.Kekulize(rwmol)
        # rdkit.Chem.rdmolops.AssignAtomChiralTagsFromStructure(rwmol, confId=0)

        # find correct conformer to assign stereochemistry
        # ideal conformer comes first
        conformers = rwmol.GetConformers()
        conformer_id = -1

        if is_degenerate_conformer(conformers[0]):
            conformer_id = conformers[1].GetId()
        else:
            conformer_id = conformers[0].GetId()

        rdkit.Chem.rdmolops.AssignStereochemistryFrom3D(rwmol, conformer_id)

    except Exception as e:
        print(e, file=sys.stderr)
        return False

    return success


def fix_molecule(rwmol: rdkit.Chem.rdchem.RWMol):
    """
    Single molecule sanitization process. Presently, only valence
    errors are taken care are of.

    Args:
        rwmol (rdkit.Chem.rdchem.RWMol): rdkit molecule to be
            sanitized

    Returns:
        bool: Whether or not sanitization succeeded
    """
    attempts = 10
    success = False
    saved_std_err = sys.stderr
    log = sys.stderr = StringIO()
    rdkit.Chem.WrapLogs()

    while (not success) and attempts >= 0:
        sanitization_result = rdkit.Chem.SanitizeMol(rwmol, catchErrors=True)

        if sanitization_result == 0:
            sys.stderr = saved_std_err
            return True

        sanitization_failures = re.findall("[a-zA-Z]{1,2}, \\d+", log.getvalue())
        if not sanitization_failures:
            sys.stderr = saved_std_err
            return False

        for sanitization_failure in sanitization_failures:

            split_object = sanitization_failure.split(",")  # [0] element [1] valency
            element = split_object[0]
            valency = int(split_object[1].strip())

            smarts_metal_check = rdkit.Chem.MolFromSmarts(
                METALS_SMART + "~[{}]".format(element)
            )
            metal_atom_bonds = rwmol.GetSubstructMatches(smarts_metal_check)
            rdkit.Chem.SanitizeMol(
                rwmol, sanitizeOps=rdkit.Chem.SanitizeFlags.SANITIZE_CLEANUP
            )

            for (metal_index, atom_index) in metal_atom_bonds:
                metal_atom = rwmol.GetAtomWithIdx(metal_index)
                erroneous_atom = rwmol.GetAtomWithIdx(atom_index)

                # change the bond type to dative
                bond = rwmol.GetBondBetweenAtoms(
                    metal_atom.GetIdx(), erroneous_atom.GetIdx()
                )
                bond.SetBondType(rdkit.Chem.BondType.SINGLE)

                if erroneous_atom.GetExplicitValence() == valency:
                    erroneous_atom.SetFormalCharge(erroneous_atom.GetFormalCharge() + 1)
                    metal_atom.SetFormalCharge(metal_atom.GetFormalCharge() - 1)

        attempts -= 1

    sys.stderr = saved_std_err

    return False


def fix_conformer(conformer):
    """In place fixing of rdkit conformer.
    In certain cases the resulting conformer (mainly from depiction process)
    can contain not valid atom coordinatesp [NaN, NaN, NaN]. This
    results in errors in downstream processes so that it is easier
    to fix it when the problem occurs

    Args:
        conformer (rdkit.Chem.rdchem.Conformer): RDKit conformer
    """
    positions = conformer.GetPositions()

    for index, pos in enumerate(positions):
        if all(np.isnan(pos)):
            new_pos = rdkit.Chem.rdGeometry.Point3D(0, 0, 0)
            conformer.SetAtomPosition(index, new_pos)
