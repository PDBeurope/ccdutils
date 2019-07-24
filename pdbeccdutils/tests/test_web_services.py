import os

import pytest
from rdkit import Chem

from pdbeccdutils.core import ccd_reader
from pdbeccdutils.tests.tst_utilities import cif_filename
from pdbeccdutils.utils.pubchem_downloader import PubChemDownloader


"""Test of routine for retrieving 2D layouts of the CCDs.
"""

ids_to_test = ['MAN', 'NAG', 'SO4', 'GOL', 'SAC', 'VIA', 'GLU']


class TestWebServices:
    @staticmethod
    @pytest.mark.parametrize('het_code', ids_to_test)
    def test_components_layouts_downloaded(tmpdir, het_code):
        dl = PubChemDownloader(tmpdir)

        to_download = os.path.join(tmpdir, f'{het_code}.sdf')
        comp = ccd_reader.read_pdb_cif_file(cif_filename(het_code)).component
        dl.process_template(comp)

        assert os.path.isfile(to_download)
        assert os.path.getsize(to_download) > 0
        assert comp.mol_no_h.HasSubstructMatch(Chem.MolFromMolFile(to_download, sanitize=True))

    @staticmethod
    @pytest.mark.parametrize('het_code', ids_to_test)
    def test_components_layouts_updated(tmpdir, het_code):
        dl = PubChemDownloader(tmpdir)

        to_download = os.path.join(tmpdir, f'{het_code}.sdf')
        dl.update_ccd_file(cif_filename(het_code))
        mol = Chem.MolFromMolFile(to_download)

        assert os.path.isfile(to_download)
        assert os.path.getsize(to_download) > 0
        assert isinstance(mol, Chem.Mol)
        assert mol.GetNumAtoms() > 0

    @staticmethod
    @pytest.mark.parametrize('het_code', ids_to_test)
    def test_unichem_download(het_code):
        c = ccd_reader.read_pdb_cif_file(cif_filename(het_code)).component

        c.fetch_external_mappings()
        assert len(c.external_mappings) > 0
        
        c.fetch_external_mappings(True)
        assert len(c.external_mappings) > 0
