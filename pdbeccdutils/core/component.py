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

import io
import re
import sys
from collections import OrderedDict
from datetime import date

import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem.rdchem import BondType
from rdkit.Chem.rdmolops import AssignAtomChiralTagsFromStructure

from pdbeccdutils.core import ConformerType, ReleaseStatus
from pdbeccdutils.helpers import IOGrabber, drawing

METALS_SMART = '[Li,Na,K,Rb,Cs,F,Be,Mg,Ca,Sr,Ba,Ra,Sc,Ti,V,Cr,Mn,Fe,Co,Ni,Cu,Zn,Al,Ga,Y,Zr,Nb,Mo,'\
               'Tc,Ru,Rh,Pd,Ag,Cd,In,Sn,Hf,Ta,W,Re,Os,Ir,Pt,Au,Hg,Tl,Pb,Bi]'


class Component:
    """
    Wrapper for the rdkit.Chem.rdchem.Mol object enabling some of its
    functionality and handling possible erroneous situations.

    Returns:
        pdbeccdutils.core.Component: instance object
    """

    def __init__(self, mol, ccd_cif_dict=None, properties=None, descriptors=None):

        self.mol = mol
        self.ccd_cif_dict = ccd_cif_dict
        self.fragments = {}
        self._2dmol = None
        self._id = ''
        self._name = ''
        self._formula = ''
        self._modified_date = ''
        self._pdbx_release_status = ReleaseStatus.NOT_SET
        self._descriptors = []
        self._inchi_from_rdkit = None
        self._inchikey_from_rdkit = None

        self.conformers_mapping = \
            {ConformerType.AllConformers: - 1,
             ConformerType.Ideal: 0,
             ConformerType.Model: 1 if len(mol.GetConformers()) == 2 else 1000,
             ConformerType.Computed: 2000}

        if properties is not None:
            mod_date = properties.modified_date.split('-')
            self._id = properties.id
            self._name = properties.name
            self._formula = properties.formula
            self._pdbx_release_status = ReleaseStatus[properties.pdbx_release_status]
            self._modified_date = date(int(mod_date[0]), int(mod_date[1]), int(mod_date[2]))

        if descriptors is not None:
            self._descriptors = descriptors

    # region properties
    @property
    def id(self):
        """
        Supply the unique identifier for the PDB-CCD,
        for example 'ATP'.
        Obtained from CCD's _chem_comp.id:

        http://mmcif.wwpdb.org/dictionaries/mmcif_std.dic/Items/_chem_comp.id.html

        If not defined then the empty string '' will be returned.

        Returns:
            str: the _chem_comp.id or ''.
        """
        return self._id

    @property
    def name(self):
        """
        Supply the 'full name' of the PDB-CCD, for example 'ETHANOL'.
        Obtained from PDB-CCD's _chem_comp.name:

        http://mmcif.wwpdb.org/dictionaries/mmcif_std.dic/Items/_chem_comp.name.html

        If not defined then the empty string '' will be returned.

        Returns:
            str: the _chem_comp.name or ''.
        """
        return self._name

    @property
    def formula(self):
        """
        Supply the chemical formula for the PDB-CCD,
        for example 'C2 H6 O'.
        Obtained from PDB-CCD's _chem_comp.formula:

        http://mmcif.wwpdb.org/dictionaries/mmcif_std.dic/Items/_chem_comp.formula.html

        If not defined then the empty string '' will be returned.

        Returns:
            str: the _chem_comp.formula or ''.
        """
        return self._formula

    @property
    def pdbx_release_status(self):
        """
        Supply the pdbx_release_status for the PDB-CCD.
        Obtained from PDB-CCD's _chem_comp.pdbx_rel_status:

        http://mmcif.wwpdb.org/dictionaries/mmcif_pdbx.dic/Items/_chem_comp.pdbx_release_status.html

        Returns:
            pdbeccdutils.core.enums.ReleaseStatus: enum of the release
            status (this includes NOT_SET if no value is defined).
        """
        return self._pdbx_release_status

    @property
    def modified_date(self):
        return self._modified_date

    @property
    def descriptors(self):
        return self._descriptors

    @property
    def inchikey(self):
        """
        Supply the InChIKey for the PDB-CCD.
        Obtained from `PDB-CCD's _pdbx_chem_comp_descriptor` table line
        with `_pdbx_chem_comp_descriptor.type=InChIKey`, see:

        http://mmcif.wwpdb.org/dictionaries/mmcif_pdbx.dic/Items/_pdbx_chem_comp_descriptor.type.html

        If not defined then the empty string '' will be returned.

        Returns:
            str: the InChIKey or ''.
        """
        return next((x.value for x in self._descriptors if x.type == 'InChIKey'), '')

    @property
    def inchi(self):
        return next((x.value for x in self._descriptors if x.type == 'InChI'), '')

    @property
    def inchi_from_rdkit(self):
        """
        provides the InChI worked out by rkdit

        Returns:
            str: the InChI or emptry '' if there was an error finding it.
        """
        if self._inchi_from_rdkit is None:
            try:
                self._inchi_from_rdkit = Chem.inchi.MolToInchi(self.mol)
            except ValueError:
                self._inchi_from_rdkit = ''
        return self._inchi_from_rdkit

    @property
    def inchikey_from_rdkit(self):
        """
        provides the InChIKey worked out by rdkit

        Returns:
            str: the InChIKey or '' if there was an error finding it.
        """
        if self._inchikey_from_rdkit is None:
            inchi = self.inchi_from_rdkit
            if inchi != 'ERROR':
                self._inchikey_from_rdkit = Chem.inchi.InchiToInchiKey(inchi)
            else:
                self._inchikey_from_rdkit = ''
            if self._inchikey_from_rdkit is None:
                self._inchikey_from_rdkit = ''
        return self._inchikey_from_rdkit

    @property
    def released(self):
        """ returns True if PDB-CCD has been released.
        Tests pdbx_release_status is REL"""
        return self._pdbx_release_status == ReleaseStatus.REL

    @property
    def mol_no_h(self):
        no_h = Chem.RemoveHs(self.mol, sanitize=False)
        Chem.SanitizeMol(no_h, catchErrors=True)
        return no_h

    @property
    def number_atoms(self):
        """
        Supplies the number of atoms in the _chem_comp_atom table

        Returns:
            int: the number of atoms in the PDB-CCD
        """
        return self.mol.GetNumAtoms()

    @property
    def atoms_ids(self):
        """
        Supplies a list of the atom_ids obtained from
        `_chem_comp_atom.atom_id`, see:

        http://mmcif.wwpdb.org/dictionaries/mmcif_pdbx.dic/Categories/chem_comp_atom.html

        The order will reflect the order in the input PDB-CCD.

        The atom_id is also also know as 'atom_name', standard amino
        acids have main chain atom names 'N CA C O'

        Returns:
            (:obj:`tuple` of :obj:`str`): `atom_id's` for the PDB-CCD
        """
        return tuple(atom.GetProp('name') for
                     atom in self.mol.GetAtoms())
    # endregion properties

    def inchikey_from_rdkit_matches_ccd(self, connectivity_only=False):
        """
        checks whether inchikey matches between ccd and rdkit

        Args:
            connectivity_only (bool): restrict to the first 14 character - the connectivity information.

        Returns:
            bool: True for match
        """
        if self.inchikey is None or self.inchikey_from_rdkit == 'ERROR':
            return False
        if connectivity_only:
            if len(self.inchikey) < 14 or len(self.inchikey_from_rdkit) < 14:
                return False
            elif self.inchikey[:14] != self.inchikey_from_rdkit[:14]:
                return False
        elif self.inchikey != self.inchikey_from_rdkit:
            return False
        return True

    def compute_2d(self, manager, remove_hs=True):
        """
        Compute 2d depiction of the component using DepictionManager
        instance

        Args:
            manager (pdbeccdutils.utils.DepictionManager):
                Instance of the ligand depiction class.
            remove_hs (bool, optional): Defaults to True. Remove
                hydrogens prior to depiction.

        Returns:
            rdkit.Chem.rdchem.Mol: 2D depiction of the ligand.
        """
        mol_copy = Chem.RWMol(self.mol)
        if remove_hs:
            mol_copy = Chem.RemoveHs(mol_copy, updateExplicitCount=True, sanitize=False)
            Chem.SanitizeMol(mol_copy, catchErrors=True)
        result_log = manager.depict_molecule(self._id, mol_copy)

        if result_log is not None:
            # add conformer
            self._2dmol = result_log.mol
            return result_log
        else:
            return None

    def export_2d_svg(self, file_name, width=500, names=False, atom_highlight={}, bond_highlight=None):
        """
        Save 2d depiction of the component as an SVG file.

        Args:
            file_name (str): path to store 2d depiction.
                width (int, optional): Defaults to 500.
                Width of a frame in pixels.
            names (bool, optional): Defaults to False. Whether or not
                to include atom names in depiction. If atom name is
                not set, element symbol is used instead.
            atomHighlight (dict, optional): Defaults to {}. Atoms names
                to be highlighted along with colors in RGB.
                eg. {'CA': (0.5, 0.5, 0.5)} or {0: (0.5, 0.5, 0.5)}
            bondHighlight (dict, optional): Defaults to None. Bonds to be
                highlighted along with colors in RGB.
                eg. {('CA', 'CB'): (0.5, 0.5, 0.5)} or
                {(0, 1): (0.5, 0.5, 0.5)}
        Raises:
            ValueError: If bond does not exist.
            KeyError: If atom does not exist.

        """
        if self._2dmol is None:
            drawing.save_no_image(file_name, width=width)
            return

        drawer = rdMolDraw2D.MolDraw2DSVG(width, width)
        atom_mapping = {self._get_atom_name(a): i for i, a in enumerate(self._2dmol.GetAtoms())}

        if all(isinstance(i, str) for i in atom_highlight.keys()):
            atom_highlight = {atom_mapping[k]: v for k, v in atom_highlight.items()}

        if bond_highlight is not None:
            if all(isinstance(i[0], str) and isinstance(i[1], str) for i in bond_highlight.keys()):
                temp_highlight = {}
                for k, v in bond_highlight.items():
                    bond = self._2dmol.GetBondBetweenAtoms(atom_mapping[k[0]], atom_mapping[k[1]])
                    if bond is None:
                        raise ValueError('Bond between {} and {} does not exist'.format(k[0], k[1]))
                    temp_highlight[bond.GetIdx()] = v
                bond_highlight = temp_highlight

        if names:
            options = drawer.drawOptions()
            for i, a in enumerate(self._2dmol.GetAtoms()):
                atom_name = self._get_atom_name(a)
                options.atomLabels[i] = atom_name
                a.SetProp('molFileAlias', atom_name)

        self._draw_molecule(drawer, file_name, width, atom_highlight, bond_highlight)

    def compute_3d(self):
        """
        Generate 3D coordinates using ETKDG method from RdKit.

        Returns:
            bool: Result of the structure generation process.
        """
        options = AllChem.ETKDGv2()
        options.clearConfs = False

        try:
            conf_id = AllChem.EmbedMolecule(self.mol, options)
            result = AllChem.UFFOptimizeMolecule(self.mol, confId=conf_id, maxIters=1000)
            self.conformers_mapping[ConformerType.Computed] = conf_id
            return True
        except RuntimeError:
            return False  # Force field issue here
        except ValueError:
            return False  # sanitization issue here

    def sanitize(self, fast=False):
        """
        Attempts to sanitize mol in place. RDKit's standard error can be
        processed in order to find out what went wrong with sanitization
        to fix the molecule.

        Args:
            fast (bool, optional): Defaults to False. If fast option is
                triggered original Oliver's sanitization process is run.


        Returns:
            bool: Result of the sanitization process.
        """
        rwmol = Chem.RWMol(self.mol)
        try:
            success = self._fix_molecule_fast(rwmol) if fast else self._fix_molecule(rwmol)

            if not success:
                return False

            Chem.Kekulize(rwmol)
            AssignAtomChiralTagsFromStructure(rwmol)
            self.mol = rwmol.GetMol()
        except Exception as e:
            print(e, file=sys.stderr)
            return False

        return success

    def has_degenerated_conformer(self, type):
        """
        Determine if given conformer has missing coordinates. This can
        be used to determine, whether or not the coordinates should be
        regenerated.

        Args:
            type (pdbeccdutils.core.ConformerType): type of coformer
                to be inspected.

        Raises:
            ValueError: If given conformer does not exist.

        Returns:
            bool: true if more then 1 atom has coordinates [0, 0, 0]
        """
        conformer = self.mol.GetConformer(self.conformers_mapping[type])
        empty_coords = Chem.rdGeometry.Point3D(0, 0, 0)
        counter = 0

        for i in range(conformer.GetNumAtoms()):
            pos = conformer.GetAtomPosition(i)
            if pos.Distance(empty_coords) == 0.0:
                counter += 1

        if counter > 1:
            return True
        return False

    def locate_fragment(self, mol):
        """
        Identify substructure match in the component.

        Args:
            mol (rdkit.Chem.rdchem.Mol): Fragment to be matched with
                structure

        Returns:
            list(list(rdkit.Chem.rdchem.Atoms)): list of fragments identified in the component as a list of Atoms.
        """
        result = []
        if mol is None:
            return []

        matches = self.mol.GetSubstructMatches(mol)

        for m in matches:
            result.append(list(map(lambda idx: self.mol.GetAtomWithIdx(idx), m)))

        return result

    def library_search(self, fragment_library):
        """Identify fragments from the fragment library in this component

        Args:
            fragment_library (pdbeccdutils.core.FragmentLibrary):
                Fragment library.

        Returns:
            int: number of matches found
        """

        matches_found = 0
        for k, v in fragment_library.library.items():
            try:
                matches = self.mol.GetSubstructMatches(v)
                matches_found += len(matches)

                if len(matches) > 0:
                    self.fragments[k] = matches
            except Exception:
                pass

        return matches_found

    def _fix_molecule(self, rwmol):
        """
        Single molecule sanitization process. Presently, only valence
        errors are taken care are of.

        Args:
            rwmol (rdkit.Chem.rdchem.Mol): rdkit molecule to be
                sanitized

        Returns:
            bool: Whether or not sanitization succeeded
        """
        attempts = 10
        success = False
        while ((not success) and attempts >= 0):
            out = IOGrabber(sys.stderr)
            out.start()
            sanitization_result = Chem.SanitizeMol(rwmol, catchErrors=True)
            if sanitization_result == 0:
                out.stop()
                return True
            out.stop()
            raw_error = out.capturedtext

            sanitization_failure = re.findall('[a-zA-Z]{1,2}, \\d+', raw_error)
            if len(sanitization_failure) == 0:
                return False

            split_object = sanitization_failure[0].split(',')  # [0] element [1] valency
            element = split_object[0]
            valency = int(split_object[1].strip())

            smarts_metal_check = Chem.MolFromSmarts(METALS_SMART + '~[{}]'.format(element))
            metal_atom_bonds = rwmol.GetSubstructMatches(smarts_metal_check)
            Chem.SanitizeMol(rwmol, sanitizeOps=Chem.SanitizeFlags.SANITIZE_CLEANUP)

            for (metal_index, atom_index) in metal_atom_bonds:
                metal_atom = rwmol.GetAtomWithIdx(metal_index)
                erroneous_atom = rwmol.GetAtomWithIdx(atom_index)

                # change the bond type to dative
                bond = rwmol.GetBondBetweenAtoms(metal_atom.GetIdx(), erroneous_atom.GetIdx())
                bond.SetBondType(BondType.DATIVE)

                if erroneous_atom.GetExplicitValence() == valency:
                    erroneous_atom.SetFormalCharge(erroneous_atom.GetFormalCharge() + 1)
                    metal_atom.SetFormalCharge(metal_atom.GetFormalCharge() - 1)

            attempts -= 1
        return False

    def _fix_molecule_fast(self, rwmol):
        """
        Fast sanitization process. Fixes just metal-N valence issues

        Args:
            rwmol (rdkit.Chem.rdchem.Mol): rdkit mol to be sanitized

        Returns:
            bool: Whether or not sanitization succeeded
        """
        smarts_metal_check = Chem.MolFromSmarts(METALS_SMART + '~[N]')
        metal_atom_bonds = rwmol.GetSubstructMatches(smarts_metal_check)
        Chem.SanitizeMol(rwmol, sanitizeOps=Chem.SanitizeFlags.SANITIZE_CLEANUP)
        for (metal_index, atom_index) in metal_atom_bonds:
            metal_atom = rwmol.GetAtomWithIdx(metal_index)
            erroneous_atom = rwmol.GetAtomWithIdx(atom_index)

            # change the bond type to dative
            bond = rwmol.GetBondBetweenAtoms(metal_atom.GetIdx(), erroneous_atom.GetIdx())
            bond.SetBondType(BondType.DATIVE)

            # change the valency
            if erroneous_atom.GetExplicitValence() == 4:
                erroneous_atom.SetFormalCharge(erroneous_atom.GetFormalCharge() + 1)
                metal_atom.SetFormalCharge(metal_atom.GetFormalCharge() - 1)

        sanitization_result = Chem.SanitizeMol(rwmol, catchErrors=True)

        return sanitization_result == 0

    def _draw_molecule(self, drawer, file_name, width, atom_highlight, bond_highlight):
        try:
            copy = rdMolDraw2D.PrepareMolForDrawing(self._2dmol, wedgeBonds=True, kekulize=True,
                                                    addChiralHs=True)
        except RuntimeError:
            copy = rdMolDraw2D.PrepareMolForDrawing(self._2dmol, wedgeBonds=False, kekulize=True,
                                                    addChiralHs=True)

        if bond_highlight is None:
            drawer.DrawMolecule(copy, highlightAtoms=atom_highlight.keys(),
                                highlightAtomColors=atom_highlight)
        else:
            drawer.DrawMolecule(copy, highlightAtoms=atom_highlight.keys(),
                                highlightAtomColors=atom_highlight,
                                highlightBonds=bond_highlight.keys(), highlightBondColors=bond_highlight)
        drawer.FinishDrawing()

        with open(file_name, 'w') as f:
            f.write(drawer.GetDrawingText())

    def _get_atom_name(self, atom):
        """Supplies atom_id obrained from `_chem_comp_atom.atom_id`, see:

        http://mmcif.wwpdb.org/dictionaries/mmcif_pdbx.dic/Categories/chem_comp_atom.html

        If there is no such atom name, it is created from the element
        symbol and atom index.

        Args:
            atom (rdkit.Chem.rdChem.Atom): rdkit atom

        Returns:
            str: atom name
        """
        return atom.GetProp('name') if atom.HasProp('name') else atom.GetSymbol() + str(atom.GetIdx())
