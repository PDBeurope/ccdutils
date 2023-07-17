```{eval-rst}
.. _similarity:
```

## Similarity

`pdbeccdutils` implements [PARITY](https://www.sciencedirect.com/science/article/pii/S0969212618300492?via%3Dihub) method to calcualte similarity between common substructures of small-molecules

```python
from pdbeccdutils.core import ccd_reader
from pdbeccdutils.computations.parity_method import compare_molecules

hem_a = ccd_reader.read_pdb_cif_file('HEA.cif').component
hem_d = ccd_reader.read_pdb_cif_file('DHE.cif').component

similarity = compare_molecules(hem_a.mol_no_h, hem_d.mol_no_h)
similarity
```
