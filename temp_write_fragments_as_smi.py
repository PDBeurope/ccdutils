from fragment_library import FragmentLibrary
frag_lib = FragmentLibrary()
for name, smiles in sorted(frag_lib.fragment_name_to_smiles.items()):
    print('{} {}'.format(smiles, name))
