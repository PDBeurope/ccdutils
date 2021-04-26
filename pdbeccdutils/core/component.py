#!/usr/bin/env python
# software from PDBe: Protein Data Bank in Europe; https://pdbe.org
#
# Copyright 2019 EMBL - European Bioinformatics Institute
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

import json
import logging
import re
from datetime import date
from typing import Any, Dict, List, Tuple

import rdkit
from rdkit.Chem import BRICS, Draw
from rdkit.Chem.rdMolDescriptors import Properties
from rdkit.Chem.Scaffolds import MurckoScaffold

from pdbeccdutils.core.depictions import DepictionManager, DepictionResult
from pdbeccdutils.core.exceptions import CCDUtilsError
from pdbeccdutils.helpers import mol_tools
from pdbeccdutils.core.fragment_library import FragmentLibrary
from pdbeccdutils.core.models import (
    CCDProperties,
    ConformerType,
    Descriptor,
    ReleaseStatus,
    ScaffoldingMethod,
    SubstructureMapping,
)
from pdbeccdutils.helpers import conversions, drawing
from pdbeccdutils.utils import web_services


class Component:
    """
    Wrapper for the rdkit.Chem.Mol object enabling some of its
    functionality and handling possible erroneous situations.

    Returns:
        Component: instance object
    """

    def __init__(
        self,
        mol: rdkit.Chem.rdchem.Mol,
        ccd_cif_dict: Dict[str, Any] = None,
        properties: CCDProperties = None,
        descriptors: List[Descriptor] = None,
    ) -> None:

        self.mol = mol
        self._mol_no_h = None
        self.mol2D = None
        self.ccd_cif_dict = ccd_cif_dict
        self._fragments: Dict[str, SubstructureMapping] = {}
        self._scaffolds: Dict[str, SubstructureMapping] = {}
        self._descriptors: List[Descriptor] = []
        self._inchi_from_rdkit = ""
        self._inchikey_from_rdkit = ""
        self._physchem_properties: Dict[str, Any] = {}
        self._external_mapping: List[Tuple[str, str]] = []

        if descriptors is not None:
            self._descriptors = descriptors

        if properties is not None:
            self._cif_properties = properties

    # region properties
    @property
    def id(self) -> str:
        """Supply the unique identifier for the PDB-CCD,
        for example 'ATP'.
        Obtained from CCD's _chem_comp.id:

        http://mmcif.wwpdb.org/dictionaries/mmcif_std.dic/Items/_chem_comp.id.html

        If not defined then the empty string '' will be returned.

        Returns:
            str: the _chem_comp.id or ''.
        """
        return self._cif_properties.id

    @property
    def name(self) -> str:
        """Supply the 'full name' of the PDB-CCD, for example 'ETHANOL'.
        Obtained from PDB-CCD's _chem_comp.name:

        http://mmcif.wwpdb.org/dictionaries/mmcif_std.dic/Items/_chem_comp.name.html

        If not defined then the empty string '' will be returned.

        Returns:
            str: the _chem_comp.name or ''.
        """
        return self._cif_properties.name

    @property
    def formula(self) -> str:
        """Supply the chemical formula for the PDB-CCD,
        for example 'C2 H6 O'.
        Obtained from PDB-CCD's _chem_comp.formula:

        http://mmcif.wwpdb.org/dictionaries/mmcif_std.dic/Items/_chem_comp.formula.html

        If not defined then the empty string '' will be returned.

        Returns:
            str: the _chem_comp.formula or ''.
        """
        return self._cif_properties.formula

    @property
    def pdbx_release_status(self) -> ReleaseStatus:
        """Supply the pdbx_release_status for the PDB-CCD.
        Obtained from PDB-CCD's _chem_comp.pdbx_rel_status:

        http://mmcif.wwpdb.org/dictionaries/mmcif_pdbx.dic/Items/_chem_comp.pdbx_release_status.html

        Returns:
            pdbeccdutils.core.enums.ReleaseStatus: enum of the release
            status (this includes NOT_SET if no value is defined).
        """
        return self._cif_properties.pdbx_release_status

    @property
    def modified_date(self) -> date:
        """Supply the _pdbx_chem_comp_descriptor category for the PDB-CCD
        Obtained from PDB-CCD's _pdbx_chem_comp_descriptor:

        http://mmcif.wwpdb.org/dictionaries/mmcif_pdbx_v50.dic/Items/_chem_comp.pdbx_modified_date.html

        Returns:
            datetime.date: Date of the last entrie's modification.
        """
        return self._cif_properties.modified_date

    @property
    def descriptors(self) -> List[Descriptor]:
        """Supply the pdbx_modified_date for the PDB-CCD
        Obtained from PDB-CCD's _chem_comp.pdbx_modified_date:

        http://mmcif.rcsb.org/dictionaries/mmcif_pdbx.dic/Items/_pdbx_chem_comp_descriptor.program_version.html

        Returns:
            list[Descriptor]: List of descriptors for a given entry.
        """
        return self._descriptors

    @property
    def inchikey(self) -> str:
        """Supply the InChIKey for the PDB-CCD.
        Obtained from `PDB-CCD's _pdbx_chem_comp_descriptor` table line
        with `_pdbx_chem_comp_descriptor.type=InChIKey`, see:

        http://mmcif.wwpdb.org/dictionaries/mmcif_pdbx.dic/Items/_pdbx_chem_comp_descriptor.type.html

        If not defined then the empty string '' will be returned.

        Returns:
            str: the InChIKey or ''.
        """
        return next((x.value for x in self._descriptors if x.type == "InChIKey"), "")

    @property
    def inchi(self) -> str:
        """Supply the InChI for the PDB-CCD.
        Obtained from `PDB-CCD's _pdbx_chem_comp_descriptor` table line
        with `_pdbx_chem_comp_descriptor.type=InChI`, see:

        http://mmcif.wwpdb.org/dictionaries/mmcif_pdbx.dic/Items/_pdbx_chem_comp_descriptor.type.html

        If not defined then the empty string '' will be returned.

        Returns:
            str: the InChI or ''.
        """
        return next((x.value for x in self._descriptors if x.type == "InChI"), "")

    @property
    def inchi_from_rdkit(self) -> str:
        """Provides the InChI worked out by RDKit.

        Returns:
            str: the InChI or empty '' if there was an error finding it.
        """
        if not self._inchi_from_rdkit:
            try:
                self._inchi_from_rdkit = rdkit.Chem.inchi.MolToInchi(self.mol)
            except ValueError:
                self._inchi_from_rdkit = ""
        return self._inchi_from_rdkit

    @property
    def inchikey_from_rdkit(self) -> str:
        """Provides the InChIKey worked out by RDKit.

        Returns:
            str: the InChIKey or '' if there was an error finding it.
        """
        if not self._inchikey_from_rdkit:
            inchi = self.inchi_from_rdkit
            if inchi != "ERROR":
                self._inchikey_from_rdkit = rdkit.Chem.inchi.InchiToInchiKey(inchi)
            else:
                self._inchikey_from_rdkit = ""
            if self._inchikey_from_rdkit is None:
                self._inchikey_from_rdkit = ""
        return self._inchikey_from_rdkit

    @property
    def released(self) -> bool:
        """Tests pdbx_release_status is REL.

        Returns:
            bool: True if PDB-CCD has been released.
        """
        return self._cif_properties.pdbx_release_status == ReleaseStatus.REL

    @property
    def mol_no_h(self) -> rdkit.Chem.rdchem.Mol:
        """RDKit mol object without hydrogens

        Returns:
            rdkit.Chem.rdchem.Mol: RDKit mol object with stripped Hs.
        """
        if self._mol_no_h is None:
            no_h = rdkit.Chem.RemoveHs(self.mol, sanitize=False)
            rdkit.Chem.SanitizeMol(no_h, catchErrors=True)
            self._mol_no_h = no_h

        return self._mol_no_h

    @property
    def number_atoms(self) -> int:
        """Supplies the number of atoms in the _chem_comp_atom table

        Returns:
            int: the number of atoms in the PDB-CCD
        """
        return self.mol.GetNumAtoms()

    @property
    def fragments(self) -> List[SubstructureMapping]:
        """Lists matched fragments and atom names.

        Returns:
            list[SubstructureMapping]: Substructure mapping for
            all discovered fragments.
        """
        return self._id_to_name_mapping(self._fragments)

    @property
    def scaffolds(self) -> List[SubstructureMapping]:
        """Lists matched scaffolds and atom names

        Returns:
            list[SubstructureMapping]: List of substructure mappings.
        """
        return self._id_to_name_mapping(self._scaffolds)

    @property
    def atoms_ids(self) -> Tuple[Any, ...]:
        """Supplies a list of the atom_ids obtained from
        `_chem_comp_atom.atom_id`, see:

        http://mmcif.wwpdb.org/dictionaries/mmcif_pdbx.dic/Categories/chem_comp_atom.html

        The order will reflect the order in the input PDB-CCD.

        The atom_id is also also know as 'atom_name', standard amino
        acids have main chain atom names 'N CA C O'

        Returns:
            tuple[str]: `atom_id's` for the PDB-CCD
        """
        return tuple(atom.GetProp("name") for atom in self.mol.GetAtoms())

    @property
    def physchem_properties(self):
        """RDKit calculated properties related to the CCD compound

        Returns:
            dict[str, float]: A list of RDKit calculated properties
        """
        if not self._physchem_properties:
            try:
                properties = Properties()
                self._physchem_properties = dict(
                    zip(
                        properties.GetPropertyNames(),
                        properties.ComputeProperties(self.mol),
                    )
                )
                self._physchem_properties["NumHeavyAtoms"] = float(
                    self.mol.GetNumHeavyAtoms()
                )
            except (RuntimeError, ValueError):
                return {}

        return self._physchem_properties

    @property
    def external_mappings(self):
        """List external mappings provided by UniChem. fetch_external_mappings()
        was not called before only agreed mapping is retrieved.

        Returns:
            list[tuple[str]]: UniChem mappings
        """

        return self._external_mapping

    @external_mappings.setter
    def external_mappings(self, value):
        """Set mapping for this component obtained with a different mean
        but internal use of UniChem.

        Args:
            list[tuple[str]]: UniChem mappings
        """
        self._external_mapping = value

    # endregion properties

    def fetch_external_mappings(self, all_mappings=False):
        """Retrieve external mapping through UniChem based on the InChi Key.

        Args:
            all_mappings (bool, optional): Get UniChem mappings. Defaults to False.

        Returns:
            dict[str, str]: Return resource ids pairing established by UniChem.
        """
        if all_mappings:
            self._external_mapping = web_services.get_all_unichem_mapping(self.inchikey)
        else:
            self._external_mapping = web_services.get_agreed_unichem_mapping(
                self.inchikey
            )

        return self._external_mapping

    def inchikey_from_rdkit_matches_ccd(self, connectivity_only: bool = False) -> bool:
        """Checks whether inchikey matches between ccd and rdkit

        Args:
            connectivity_only (bool): restrict to the first 14 character - the connectivity information.

        Returns:
            bool: True for match
        """
        if self.inchikey is None or self.inchikey_from_rdkit == "ERROR":
            return False
        if connectivity_only:
            if len(self.inchikey) < 14 or len(self.inchikey_from_rdkit) < 14:
                return False
            if self.inchikey[:14] != self.inchikey_from_rdkit[:14]:
                return False

        if self.inchikey != self.inchikey_from_rdkit:
            return False

        return True

    def compute_2d(
        self, manager: DepictionManager, remove_hs: bool = True
    ) -> DepictionResult:
        """Compute 2d depiction of the component using DepictionManager
        instance.

        Args:
            manager (DepictionManager): Instance of the ligand depiction
                class.
            remove_hs (bool, optional): Defaults to True. Remove
                hydrogens prior to depiction.

        Returns:
            DepictionResult: Object with the details about depiction process.
        """
        mol_copy = rdkit.Chem.RWMol(self.mol)
        if remove_hs:
            mol_copy = rdkit.Chem.RemoveHs(
                mol_copy, updateExplicitCount=True, sanitize=False
            )
            rdkit.Chem.SanitizeMol(mol_copy, catchErrors=True)

        result_log = manager.depict_molecule(self.id, mol_copy)

        self.mol2D = result_log.mol

        return result_log

    def export_2d_svg(
        self,
        file_name: str,
        width: int = 500,
        names: bool = False,
        wedge_bonds: bool = True,
        atom_highlight: Dict[Any, Tuple] = None,
        bond_highlight: Dict[Tuple, Tuple] = None,
    ):
        """Save 2D depiction of the component as an SVG file. Component
        id is generated in case the image cannot be drawn.

        Args:
            file_name (str): path to store 2d depiction
            width (int, optional): Defaults to 500. Width of a frame in pixels.
            names (bool, optional): Defaults to False. Whether or not to
                include atom names in depiction. If atom name is not set, element symbol is used instead.
            wedge_bonds (bool, optional): Defaults to True. Whether or not
                the molecule should be depicted with bond wedging.
            atomHighlight (:obj:`dict` of :obj:`tuple` of :obj:`float`, optional):
                Defaults to None. Atoms names to be highlighted along
                with colors in RGB. e.g. {'CA': (0.5, 0.5, 0.5)} or {0: (0.5, 0.5, 0.5)}
            bondHighlight (:obj:`dict` of :obj:`tuple` of :obj:`float`, optional):
                Defaults to None. Bonds to be highlighted along with
                colors in RGB. e.g. {('CA', 'CB'): (0.5, 0.5, 0.5)} or {(0, 1): (0.5, 0.5, 0.5)}

        Raises:
            CCDUtilsError: If bond or atom does not exist.
        """
        get_atom_name = (
            lambda a: a.GetProp("name")
            if a.HasProp("name")
            else a.GetSymbol() + str(a.GetIdx())
        )

        if self.mol2D is None:
            drawing.save_no_image(file_name, self.id, width)
            return

        drawer = Draw.rdMolDraw2D.MolDraw2DSVG(width, width)
        options = drawer.drawOptions()
        atom_mapping = {
            get_atom_name(a): i for i, a in enumerate(self.mol2D.GetAtoms())
        }

        atom_highlight = {} if atom_highlight is None else atom_highlight
        bond_highlight = {} if bond_highlight is None else bond_highlight

        if width < 201:
            options.bondLineWidth = 1

        if all(isinstance(i, str) for i in atom_highlight.keys()):
            atom_highlight = {atom_mapping[k]: v for k, v in atom_highlight.items()}
        else:
            atom_highlight = {}

        if bond_highlight:
            if all(
                isinstance(i[0], str) and isinstance(i[1], str)
                for i in bond_highlight.keys()
            ):
                temp_highlight = {}
                for k, v in bond_highlight.items():
                    bond = self.mol2D.GetBondBetweenAtoms(
                        atom_mapping[k[0]], atom_mapping[k[1]]
                    )
                    if bond is None:
                        raise CCDUtilsError(
                            "Bond between {} and {} does not exist".format(k[0], k[1])
                        )
                    temp_highlight[bond.GetIdx()] = v
                bond_highlight = temp_highlight

        if names:
            for i, a in enumerate(self.mol2D.GetAtoms()):
                atom_name = get_atom_name(a)
                options.atomLabels[i] = atom_name
                a.SetProp("molFileAlias", atom_name)

        drawing.draw_molecule(
            self.mol2D, drawer, file_name, wedge_bonds, atom_highlight, bond_highlight
        )

    def export_2d_annotation(self, file_name: str, wedge_bonds: bool = True) -> None:
        """Generates 2D depiction in JSON format with annotation of
        bonds and atoms to be redrawn in the interactions component.

        Args:
            file_name (str): Path to the file
        """
        w, h = drawing.get_drawing_scale(self.mol2D)
        drawer = Draw.MolDraw2DSVG(w, h)
        drawer.drawOptions().includeAtomTags = True
        try:
            tmp = rdkit.Chem.Draw.PrepareMolForDrawing(
                self.mol2D, wedgeBonds=wedge_bonds, kekulize=True, addChiralHs=False
            )
        except (RuntimeError, ValueError):
            tmp = rdkit.Chem.Draw.PrepareMolForDrawing(
                self.mol2D, wedgeBonds=False, kekulize=True, addChiralHs=False
            )
        drawer.DrawMolecule(tmp)
        drawer.FinishDrawing()
        svg = drawer.GetDrawingText()
        json_repr = drawing.convert_svg(svg, self.id, self.mol2D)

        with open(file_name, "w") as fp:
            json.dump(json_repr, fp, indent=4, sort_keys=True)

    def compute_3d(self, version="v3") -> bool:
        """Generate 3D coordinates using EKTDG method. Version can be specified

        Args:
            version (str, optional): Version of EKTDG to be used.
                Defaults to "v3".

        Returns:
            bool: Result of the structure generation process.
        """

        if version == "v3":
            options = rdkit.Chem.AllChem.ETKDGv3()
        elif version == "v2":
            options = rdkit.Chem.AllChem.ETKDGv2()
        else:
            options = rdkit.Chem.AllChem.ETKDGv2()

        options.clearConfs = False
        conf_id = -1

        try:
            conf_id = rdkit.Chem.AllChem.EmbedMolecule(self.mol, options)
            rdkit.Chem.AllChem.UFFOptimizeMolecule(
                self.mol, confId=conf_id, maxIters=1000
            )

        except RuntimeError:
            pass  # Force field issue here
        except ValueError:
            pass  # sanitization issue here

        if conf_id != -1:
            conformer = self.mol.GetConformer(conf_id)
            conformer.SetProp("name", ConformerType.Computed.name)
            conformer.SetProp("coord_generation", f"ETKDG{version}")

            return True

        return False

    def has_degenerated_conformer(self, c_type: ConformerType) -> bool:
        """
        Determine if given conformer has missing coordinates or is
        missing completelly from the rdkit.Mol object. This can
        be used to determine, whether or not the coordinates should be
        regenerated.

        Args:
            type (ConformerType): type of conformer
                to be inspected.

        Returns:
            bool: True if more then 1 atom has coordinates [0, 0, 0]
        """
        ok_conformer = False

        for c in self.mol.GetConformers():
            if c.GetProp("name") == c_type.name:
                ok_conformer = mol_tools.is_degenerate_conformer(c)
                break

        return ok_conformer

    def get_conformer(self, c_type) -> rdkit.Chem.rdchem.Conformer:
        """
        Retrieve an rdkit object for a deemed conformer.

        Args:
            c_type (ConformerType): Conformer type to be retrieved.

        Raises:
            ValueError: If conformer does not exist

        Returns:
            rdkit.Chem.rdchem.Conformer: RDKit conformer object
        """
        for c in self.mol.GetConformers():
            try:
                if c.GetProp("name") == c_type.name:
                    return c
            except KeyError:
                pass

        raise ValueError(f"Conformer {c_type.name} does not exist.")

    def locate_fragment(
        self, mol: rdkit.Chem.rdchem.Mol
    ) -> List[List[rdkit.Chem.rdchem.Atom]]:
        """
        Identify substructure match in the component.

        Args:
            mol (rdkit.Chem.rdchem.Mol): Fragment to be matched with
                structure

        Returns:
            list[list[rdkit.Chem.rdchem.Atom]]: List of fragments
            identified in the component as a list of atoms.
        """
        result = []
        if mol is None:
            return []

        matches = self.mol.GetSubstructMatches(mol)

        for m in matches:
            result.append([self.mol.GetAtomWithIdx(idx) for idx in m])

        return result

    def library_search(
        self, fragment_library: FragmentLibrary
    ) -> List[SubstructureMapping]:
        """Identify fragments from the fragment library in this component

        Args:
            fragment_library (FragmentLibrary): Fragment library.

        Returns:
            list[SubstructureMapping]: Matches found in this run
        """
        temp = {}
        for v in fragment_library.library.values():
            try:
                matches = self.mol_no_h.GetSubstructMatches(v.mol)

                if not matches:
                    continue

                key = f"{fragment_library.name}_{v.name}"
                if key not in self._fragments:
                    temp[key] = SubstructureMapping(
                        v.name, rdkit.Chem.MolToSmiles(v.mol), v.source, matches
                    )

            except Exception:
                logging.warning(f"Error mapping fragment {v.name}.")

        self._fragments.update(temp)

        return list(temp.values())

    def get_scaffolds(self, scaffolding_method=ScaffoldingMethod.MurckoScaffold):
        """Compute deemed scaffolds for a given compound.

        Args:
            scaffolding_method (ScaffoldingMethod, optional):
                Defaults to MurckoScaffold. Scaffolding method to use

        Returns:
            list[rdkit.Chem.rdchem.Mol]: Scaffolds found in the component.
        """
        try:
            scaffolds = []

            if scaffolding_method == ScaffoldingMethod.MurckoScaffold:
                scaffolds = [(MurckoScaffold.GetScaffoldForMol(self.mol_no_h))]

            elif scaffolding_method == ScaffoldingMethod.MurckoGeneric:
                scaffolds = [(MurckoScaffold.MakeScaffoldGeneric(self.mol_no_h))]

            elif scaffolding_method == ScaffoldingMethod.Brics:
                scaffolds = BRICS.BRICSDecompose(self.mol_no_h)
                brics_smiles = [
                    re.sub(r"(\[[0-9]*\*\])", "[H]", i) for i in scaffolds
                ]  # replace dummy atoms with H's to get matches https://sourceforge.net/p/rdkit/mailman/message/35261974/
                brics_mols = [rdkit.Chem.MolFromSmiles(x) for x in brics_smiles]

                for mol in brics_mols:
                    rdkit.Chem.RemoveHs(mol)

                brics_hits = [self.mol_no_h.GetSubstructMatches(i) for i in brics_mols]

                for index, brics_hit in enumerate(brics_hits):
                    smiles = rdkit.Chem.MolToSmiles(brics_mols[index])
                    name = scaffolding_method.name
                    source = "RDKit scaffolds"
                    key = f"{name}_{smiles}"
                    brics_hit = conversions.listit(brics_hit)

                    if not smiles:
                        continue

                    if key not in self._scaffolds:
                        self._scaffolds[key] = SubstructureMapping(
                            name, smiles, source, brics_hit
                        )

                return brics_mols

            for s in scaffolds:
                scaffold_atom_names = [atom.GetProp("name") for atom in s.GetAtoms()]
                mapping = []
                for at_name in scaffold_atom_names:
                    idx = [
                        atom.GetIdx()
                        for atom in self.mol.GetAtoms()
                        if atom.GetProp("name") == at_name
                    ][0]
                    mapping.append(idx)

                smiles = rdkit.Chem.MolToSmiles(s)
                name = scaffolding_method.name
                source = "RDKit scaffolds"

                if not smiles:
                    continue

                if name in self._scaffolds:
                    self._scaffolds[name].mappings.append(mapping)
                else:
                    self._scaffolds[name] = SubstructureMapping(
                        name, smiles, source, [mapping]
                    )

            return scaffolds

        except (RuntimeError, ValueError):
            raise CCDUtilsError(
                f"Computing scaffolds using method {scaffolding_method.name} failed."
            )

    def _id_to_name_mapping(self, struct_mapping):
        """Lists matched scaffolds and atom names

        Args:
            struct_mapping (dict[str, SubstructureMapping]): Basic mapping
                for fragments/scaffolds.

        Returns:
            dict[str, SubstructureMapping]: Dictionary with scaffold names
            and matched atoms.
        """
        res = []

        for v in struct_mapping.values():
            mappings = []

            for m in v.mappings:
                atom_names = [self.mol.GetAtomWithIdx(idx).GetProp("name") for idx in m]
                mappings.append(atom_names)
            res.append(SubstructureMapping(v.name, v.smiles, v.source, mappings))

        return res
