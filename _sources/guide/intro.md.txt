```eval_rst
.. _intro:
```

# Introduction

The `pdbeccdutils` is presently under development and new functionality is added regularly as well as its functionality is being revised and updated. When properly installed all the code should have documentation. All the *public* methods do have [static typing](http://mypy-lang.org/) introduced in Python 3.5. All the interfaces should be well documented.

## Installation

The `pdbeccdutils` can be presently obtained from the [EBI Gitlab](https://gitlab.ebi.ac.uk/pdbe/ccdutils.git) using the following command:

```console
pip install git+https://gitlab.ebi.ac.uk/pdbe/ccdutils.git
```

If you want to contribute to the project please fork it first and then do a pull request.

# Getting started

The centerpoint of the `pdbecccdutils` package is a `Component` object, which is a wrapper around the default `rdkit.Chem.rdchem.Mol` object (object property `,mol`) providing most of the functionality and access to its properties. All the conformers are stored in the `rdkit.Chem.rdchem.Mol` with the exception of 2D depiction, as this one does not contain hydrogen atoms. `pdbeccdutils.core.modes.ConformerType` object allows accessing all of them.

Below you can find a few typical use cases.

## Reading CCD mmCIF files

### A single component file

```python
from pdbeccdutils.core import ccd_reader

component = ccd_reader.read_pdb_cif_file('HEM.cif').component
rdkit_mol = component.mol
```

The `rdkit.Chem.rdchem.Mol` object is sanitized already.

### Chemical component dictionary 
Chemical component dictionary can be read in a single command and `rdkit.Chem.rdchem.Mol` representations obtained immediately. Resulting data structure of `ccd_reader.read_pdb_components_file` function is `Dict<str,pdbeccdutils.core.Component>` keyed on component ID as provided by the `data_XXX` element in the respective mmCIF file.

```python
from pdbeccdutils.core import ccd_reader

parsed_components = ccd_reader.read_pdb_components_file('components.cif')
ccd_components = [v.component for k, v in parsed_components.items()]
rdkit_mols = [k for k in components.mol]
```

The `rdkit.Chem.rdchem.Mol` objects are sanitized already!

## Writing CCD files

```python
from pdbeccdutils.core import ccd_writer
from pdbeccdutils.core.models import ConformerType

# write idealized coordinates in the SDF format.
ccd_writer.write_molecule('HEM.sdf', component)

# write model coordinates in the PDB format without hydrogens
ccd_writer.write_molecule('HEM.pdb', component, remove_hs=True, conf_type=ConformerType.Model)

# write model coordinates in the mmCIF format hydrogens
ccd_writer.write_molecule('HEM.cif', component, conf_type=ConformerType.Model)
```
