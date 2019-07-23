import json
import os
import xml.etree.ElementTree as ET

import pytest
from pdbeccdutils.core import ccd_writer
from pdbeccdutils.core.models import ConformerType
from rdkit import Chem


class TestFileWrites:

    @staticmethod
    def test_write_sdf(component, tmpdir_factory):
        wd = tmpdir_factory.mktemp('sdf_test')
        for ideal in True, False:
            for remove_hs in True, False:
                suffix = f'{("" if remove_hs else "H")}'
                sdf_file = os.path.join(wd, f'{component.id}_{suffix}.sdf')
                conf_type = ConformerType.Ideal if ideal else ConformerType.Model
                ccd_writer.write_molecule(path=sdf_file,
                                          component=component,
                                          conf_type=conf_type,
                                          remove_hs=remove_hs)
                rdkit_mol = component.mol_no_h if remove_hs else component.mol

                assert os.path.isfile(sdf_file)
                assert os.path.getsize(sdf_file) > 0

                mol = Chem.MolFromMolFile(sdf_file, sanitize=False)

                assert isinstance(mol, Chem.rdchem.Mol)
                assert mol.GetNumAtoms() == rdkit_mol.GetNumAtoms()

    @staticmethod
    def test_pdb_write(component, tmpdir_factory):
        wd = tmpdir_factory.mktemp('pdb_test')

        if component.id in ('ASX', 'NA'):
            return

        for ideal in True, False:
            for remove_h in True, False:
                suffix = "" if remove_h else "H"
                pdb_file = os.path.join(wd, f'{component.id}{suffix}.pdb')
                conformer = ConformerType.Ideal if ideal else ConformerType.Model
                ccd_writer.write_molecule(pdb_file, component, remove_hs=remove_h, conf_type=conformer)
                rdkit_mol = component.mol_no_h if remove_h else component.mol

                assert os.path.isfile(pdb_file)
                assert os.path.getsize(pdb_file) > 0

                mol = Chem.MolFromPDBFile(pdb_file, removeHs=False, sanitize=False, )

                assert isinstance(mol, Chem.rdchem.Mol)
                assert mol.GetNumAtoms() == rdkit_mol.GetNumAtoms()

    @staticmethod
    def test_pdb_fallback_write(component):
        if component.id == 'ASX':
            return
        
        for ideal in True, False:
            for remove_h in True, False:
                conformer = ConformerType.Ideal if ideal else ConformerType.Model
                conf_id = component.conformers_mapping[conformer]
                rdkit_mol = component.mol_no_h if remove_h else component.mol

                pdb_repr = ccd_writer._to_pdb_str_fallback(rdkit_mol, component.id, conf_id)

                assert pdb_repr
                mol = Chem.MolFromPDBBlock(pdb_repr, sanitize=False, removeHs=False)

                assert isinstance(mol, Chem.rdchem.Mol)
                assert mol.GetNumAtoms() == rdkit_mol.GetNumAtoms()

    @staticmethod
    def test_write_xyz(component, tmpdir_factory):
        wd = tmpdir_factory.mktemp('xyz_test')
        xyz_file = os.path.join(wd, f'{component.id}.xyz')
        ccd_writer.write_molecule(xyz_file, component)

        assert os.path.isfile(xyz_file)
        assert os.path.getsize(xyz_file) > 0

        with open(xyz_file, 'r') as fp:
            lines = fp.read().splitlines()
            assert lines[1] == component.id
            assert int(lines[0]) == len(lines[2:])

    @staticmethod
    def test_xml_write(component, tmpdir_factory):
        wd = tmpdir_factory.mktemp('xml_test')
        xml_file = os.path.join(wd, f'{component.id}.xml')

        ccd_writer.write_molecule(xml_file, component)

        assert os.path.isfile(xml_file)
        assert os.path.getsize(xml_file) > 0

        xml = ET.parse(xml_file).getroot()

        assert xml.tag == 'chemComp'
        assert xml.find('id').text == component.id
        assert xml.find('name').text == component.name
        assert xml.find('formula').text == component.formula
        assert xml.find('systematicName').text == next((x.value for x in component.descriptors
                                                        if x.type == 'SYSTEMATIC NAME' and x.program == 'ACDLabs'), None)
        assert xml.find('stereoSmiles').text == next((x.value for x in component.descriptors
                                                      if x.type == 'SMILES_CANONICAL' and x.program == 'CACTVS'), None)
        assert xml.find('nonStereoSmiles').text == next((x.value for x in component.descriptors
                                                         if x.type == 'SMILES' and x.program == 'CACTVS'), None)
        assert xml.find('inchi').text == (component.inchi if component.inchi else None)

    @staticmethod
    def test_json_write(component, tmpdir_factory):
        wd = tmpdir_factory.mktemp('json_test')

        for remove_h in True, False:
            suffix = "" if remove_h else "H"
            json_file = os.path.join(wd, f'{component.id}{suffix}.json')
            ccd_writer.write_molecule(json_file, component, remove_hs=remove_h)
            rdkit_mol = component.mol_no_h if remove_h else component.mol

            assert os.path.isfile(json_file)
            assert os.path.getsize(json_file) > 0

            with open(json_file, 'r') as fp:
                js = json.load(fp)
                key = list(js.keys())[0]

                assert key == component.id
                assert len(js[key]['atoms']) == rdkit_mol.GetNumAtoms()
                assert len(js[key]['bonds']) == rdkit_mol.GetNumBonds()
