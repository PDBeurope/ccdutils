```{eval-rst}
.. _fragments:
```

# Fragments

`pdbeccdutils` comes with the code to search components using common fragments. Presently, the fragment library contains 2158 fragments which were manually currated at PDBe and collaborating resources (ENAMINE, DSI). Should you wish to check all the fragments which come with code, please check `pdbeccdutils/data/fragment_library.tsv` file.

Alternativelly fragments can be supplied in an external library (*.tsv) provided the library is formatted accordingly:

|   name	|   kind	|   query	|   description |	comment |   url |	source  |
|-----------|-----------|-----------|---------------|-----------|-------|------------|
|   acetylurea  |   SMILES  |	C1C(=O)NC(=O)N1 |       |   unchecked   |		|PDBe |
| phenanthrene | SMARTS | [#6]1:[#6]:[#6]:[#6]2:[#6](:[#6]:1):[#6]:[#6]:[#6]1:[#6]:2:[#6]:[#6]:[#6]:[#6]:1 | | unchecked | | PDBe |


## Identifying fragments of a chemical component

```python
from pdbeccdutils.core import ccd_reader
from pdbeccdutils.core.fragment_library import FragmentLibrary

component = ccd_reader.read_pdb_cif_file('HEL.cif').component
fragment_library = FragmentLibrary()

matches = component.library_search(fragment_library)
print(f'Matches found in the fragment library {matches}.')

fragment_mols = [fragment.mol for fragment in component.fragments]
img = Draw.MolsToGridImage(fragment_mols, legends = [fragment.name for fragment in component.fragments])
img
```
<img src='../_static/fragment_example.svg' style="display:block margin-bottom:5px" />  


## Identifying all chemical components with penicillin fragment

```python
fragment_library = FragmentLibrary()
ccd_dict = ccd_reader.read_pdb_components_file('components.cif')
ccd_with_penicillin_fragment = []
for ccd_id in ccd_dict.keys():
    component = ccd_dict[ccd_id].component
    frag_matches = component.library_search(fragment_library)
    for fragment in component.fragments:
        if fragment.name == 'penicillin':
            ccd_with_penicillin_fragment.append(ccd_id)

ccd_with_penicillin_fragment

['0RN', 'AIC', 'APV', 'CXN', 'HEL', 'IP1', 'MII', 'NFN', 'PN1', 'PNN', 'PNV', 'SOX', 'TAZ', 'WPP', 'X1E']
```
## PDBe supplied fragments

Below you can find actual fragment structures comming with the pdbeccdutil's `FragmentsLibrary` from the PDBe resource:

<img src='../_static/pdbe_fragments.svg' style="display:block"/>
