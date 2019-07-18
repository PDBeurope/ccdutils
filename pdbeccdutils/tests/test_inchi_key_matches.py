sample_ccd_with_inchi_problems = ['7OM', 'ASX', 'CDL', '0OD']


class TestIn:
    @staticmethod
    def test_inchikeys_from_rdkit_and_ccd_match(component):

        if component.id not in sample_ccd_with_inchi_problems:
            assert component.inchikey == component.inchikey_from_rdkit
