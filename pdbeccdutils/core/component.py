import io
import re
import sys
from datetime import date

import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem.rdchem import BondType
from rdkit.Chem.rdmolops import AssignAtomChiralTagsFromStructure
from wurlitzer import pipes

from pdbeccdutils.extensions import drawing
from pdbeccdutils.core.enums import ConformerType

METALS_SMART = '[Li,Na,K,Rb,Cs,F,Be,Mg,Ca,Sr,Ba,Ra,Sc,Ti,V,Cr,Mn,Fe,Co,Ni,Cu,Zn,Al,Ga,Y,Zr,Nb,Mo,Tc,Ru,Rh,Pd,Ag,' \
               'Cd,In,Sn,Hf,Ta,W,Re,Os,Ir,Pt,Au,Hg,Tl,Pb,Bi]'


class Component:
    """Wrapper for the rdkit.Chem.rdchem.Mol object enabling some of its functionality and
    handling possible erroneous situations.

        Presently implemented:
            * sanitization
            * 2D depiction
            * 3D conformation calculation
        TODO:
            *cml generation
            *

    Returns:
        pdbeccdutils.utils.Component -- instance object
    """

    def __init__(self, mol, properties=None, descriptors=None):

        self.mol = mol
        self._2dmol = None
        self._id = ''
        self._name = ''
        self._formula = ''
        self._modified_date = ''
        self._released = False
        self._descriptors = []

        self._conformers_mapping = {ConformerType.Ideal: 0,
                                    ConformerType.Model: 1 if len(mol.GetConformers()) == 2 else 1000,
                                    ConformerType.Computed: 2000}

        if properties is not None:
            mod_date = properties.modified_date.split('-')

            self._id = properties.id
            self._name = properties.name
            self._formula = properties.formula
            self._released = properties.released == 'REL'
            self._modified_date = date(int(mod_date[0]), int(mod_date[1]), int(mod_date[2]))

        if descriptors is not None:
            self._descriptors = descriptors

    # region properties
    @property
    def id(self):
        return self._id

    @property
    def name(self):
        return self._name

    @property
    def formula(self):
        return self._formula

    @property
    def modified_date(self):
        return self._modified_date

    @property
    def descriptors(self):
        return self._descriptors

    @property
    def inchikey(self):
        return next((x.value for x in self._descriptors if x.type == 'InChIKey'), None)

    @property
    def released(self):
        self._released
    # endregion properties

    def export_mol_representation(self, remove_hs=True, str_format='sdf', conf_type=None):
        """Export a component representation in given format. If conf_type is None, all
        the conformers are exported

        Keyword Arguments:
            remove_hs {bool} -- [] (default: {True})
            str_format {str} -- [description] (default: {'sdf'})
            conf_type {pdbeccdutils.utils.ConformerType} -- Which conformer to retrieve. Options:
            Model|Ideal|Computed (default: {ConformerType.Model}). If None is pased all the
            conformers are exported.

        Raises:
            NotImplementedError -- for mmcif option. This will soon to follow.
            ValueError -- for invalid or no exiting conformer and unknown format.
        """
        try:
            conf_id = -1
            mol_to_save = self._2dmol if conf_type == ConformerType.Depiction else self.mol
            mol_to_save = Chem.RemoveHs(mol_to_save) if remove_hs else mol_to_save

            if conf_type is not None:
                conf_id = int(conf_type)

            if str_format == 'sdf':
                return Chem.MolToMolBlock(mol_to_save, confId=conf_id)
            elif str_format == 'pdb':
                return Chem.MolToPDBBlock(mol_to_save, confId=conf_id)
            elif str_format == 'mmcif':
                raise NotImplementedError('Not yet implementd')
            else:
                raise ValueError('Unknown file format: {}'.format(str_format))
        except KeyError:
            raise ValueError('Option {} is not valid for a conformer type.'.format(conf_type))
        except ValueError:
            raise ValueError('Conformer {} does not exist, perhaps it never did.'.format(conf_type))

    def compute_2d(self, manager, remove_hs=True):
        """Compute 2d depiction of the component using Flattening manager instance

        Arguments:
            manager {pdbeccdutils.utils.FlatteningManager} -- instance of the flattening class.

        Keyword Arguments:
            remove_hs {bool} -- Remove hydrogens prior depiction. (default: {True})
        """
        mol_copy = Chem.RWMol(self.mol)
        if remove_hs:
            mol_copy = Chem.RemoveHs(mol_copy, updateExplicitCount=True, sanitize=False)
        result_log = manager.flaten_molecule(self._id, mol_copy)

        if result_log is not None:
            # add conformer
            self._2dmol = result_log.mol
            return result_log
        else:
            return None

    def export_2d_svg(self, file_name, width=500, names=False):
        """Save 2d depiction of the component as an SVG file.

        Arguments:
            file_name {str} -- path to store 2d depiction.

        Keyword Arguments:
            width {int} -- Width of a frame in pixels. (default: {500})
            names {bool} -- Whether or not to include atom names in depiction.
            if atom name is not set. element symbol is used instead. (default: {False})
        """
        drawer = rdMolDraw2D.MolDraw2DSVG(width, width)

        if names:
            options = drawer.drawOptions()
            for atom in self._2dmol.GetAtoms():
                atom_name = atom.GetProp('name') if atom.HasProp('name') == 1 else atom.GetSymbol()
                options.atomLabels[atom.GetIdx()] = atom.GetProp('name')
                atom.SetProp('molFileAlias', atom_name)

        if self._2dmol is None:
            drawing.save_no_image(file_name, width=width)
        else:
            copy = rdMolDraw2D.PrepareMolForDrawing(self._2dmol, wedgeBonds=True, kekulize=True, addChiralHs=True)
            drawer.DrawMolecule(copy)
            drawer.FinishDrawing()
            with open(file_name, 'w') as f:
                f.write(drawer.GetDrawingText())

    def compute_3d(self):
        """Generate 3D coordinates using ETKDG method from RdKit

           Returns:
            bool -- Whether or not the structure generation was succesfull
        """
        options = AllChem.ETKDG()
        options.clearConfs = False

        try:
            conf_id = AllChem.EmbedMolecule(self.mol, options)
            result = AllChem.UFFOptimizeMolecule(self.mol, confId=conf_id, maxIters=1000)
            self._conformers_mapping[ConformerType.Computed] = conf_id
            return True
        except RuntimeError:
            return False  # Force field issue here
        except ValueError:
            return False  # sanitization issue here

    def sanitize(self, fast=True):
        """Attempts to sanitize mol in place. RDKit's standard error can be processed in order to
        find out what went wrong with sanitization to fix the molecule. Presently, that function is
        slow (perhaps due to the stream?) and needs more debuging to find out what's happening there.

        Keyword Arguments:
            fast {bool} -- TEMPORARY until the speed issue is fixed. (default: {true})

        Returns:
            bool -- Whether the sanitization procedure has been succesfull.
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

    def is_degenerated_conformer(self, type):
        """Determine if given conformer has missing coordinates. This can be used to
        determine, whether or not the coordinates should be regenerated.

        Arguments:
            type {ConformerType} -- type of coformer to be inspected

        Raises:
            {ValueError} -- If given conformer does not exist.
        Returns:
            [bool] -- true if more 1 atom has coordinates [0, 0, 0]
        """
        conformer = self.mol.GetConformer(self._conformers_mapping[type])
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
        """Identify substructure match in the component.

        Arguments:
            mol {rdkit.Chem.rdchem.Mol} -- Molecule to be matched with structure

        Returns:
            [list(list(rdkit.Chem.rdchem.Atoms))] -- list of fragments identifiedd in the
            component as a list of rdkit.Chem.rdchem.Atoms.
        """
        result = []
        if mol is None:
            return []

        matches = self.mol.GetSubstructMatches(mol)

        for m in matches:
            result.append(list(map(lambda idx: self.mol.GetAtomWithIdx(idx), m)))

        return result

    def _fix_molecule(self, rwmol):
        """Single molecule sanitization process. Presently only valence errors are taken care are of.
        Capture the C level stream has presently HUGE tradeoff ~ 200ms so should be used with care.

        Arguments:
            rwmol {rdkit.Chem.rdchem.Mol} -- rdkit molecule to be sanitized

        Returns:
            [bool] -- Whether or not sanitization succeeded
        """
        attempts = 10
        success = False

        with io.BytesIO() as stream:
            with pipes(stderr=stream):
                while ((not success) and attempts >= 0):
                    sanitization_result = Chem.SanitizeMol(rwmol, catchErrors=True)

                    if sanitization_result == 0:
                        return True

                    raw_error = stream.getvalue().decode('utf-8')
                    stream.seek(0)
                    stream.truncate(0)

                    sanitization_failure = re.findall('[a-zA-Z]{1,2}, \d+', raw_error)
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

                        # change the valency
                        name = erroneous_atom.GetProp('name')
                        vale = erroneous_atom.GetExplicitValence()
                        total_valence = erroneous_atom.GetTotalValence()
                        explicit = erroneous_atom.GetExplicitValence()
                        degree = erroneous_atom.GetTotalDegree()

                        metal_charge = metal_atom.GetFormalCharge()

                        if erroneous_atom.GetExplicitValence() == valency:
                            erroneous_atom.SetFormalCharge(erroneous_atom.GetFormalCharge() + 1)
                            metal_atom.SetFormalCharge(metal_atom.GetFormalCharge() - 1)

                    attempts -= 1
        return False

    def _fix_molecule_fast(self, rwmol):
        """Fast sanitization process. Fixes just metal-N valence issues

        Arguments:
            rwmol {rdkit.Chem.rdchem.Mol} -- rdkit molecule to be sanitized

        Returns:
            [bool] -- Whether or not sanitization succeeded
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
