pdbeccdutils.scripts
====================

Protein-ligand interaction pipeline
-----------------------------------
This script is used to generate protein-ligand interactions information
for bound molecules from mmcif files

During this process `_nonpoly_seq` table is used to to identify ligands
and `_struct_conn` table to come up with connectivity among them. In the
next step ChimeraX is used to protonate structures and refactored
version of Arpeggio software [1] is used to come up with protein-ligand
interactions.

[1] Jubb, H. C., Higueruelo, A. P., Ochoa-Montaño, B., Pitt, W. R.,
Ascher, D. B., & Blundell, T. L. (2017). [Arpeggio: A Web Server for Calculating and Visualising Interatomic Interactions in Protein Structures](https://doi.org/10.1016/j.jmb.2016.12.004).
Journal of Molecular Biology, 429(3), 365–371.

The refactored version can be installed from source using:
```bash
git clone -b development https://github.com/lpravda/arpeggio.git 
pip install -e arpeggio
```

```eval_rst
.. automodule:: pdbeccdutils.scripts.interactions_cli
    :members:
````

PDBeChem pipeline
-----------------
```eval_rst
.. automodule:: pdbeccdutils.scripts.process_components_cif_cli
    :members:
```    

Setup pubchem library
---------------------
```eval_rst
.. automodule:: pdbeccdutils.scripts.setup_pubchem_library_cli
    :members:
```    
