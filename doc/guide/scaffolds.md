```{eval-rst}
.. _scaffolds:
```
# Scaffolds

`pdbeccdutils` helps to search for scaffolds of small molecules in PDB using RDKit. Presently, the following scaffold identification methods from RDKit are supported: `MurckoScaffold`, `MurckoGeneric`, `Brics`

```python
from pdbeccdutils.core import ccd_reader

component = ccd_reader.read_pdb_cif_file('HEM.cif').component
component.get_scaffolds()
```
