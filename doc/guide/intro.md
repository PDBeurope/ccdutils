```{eval-rst}
.. _intro:
```

# Introduction

The `pdbeccdutils` is under development and new functionality is added regularly as well as its functionality is being revised and updated. When properly installed all the code should have documentation. All the *public* methods do have [static typing](http://mypy-lang.org/) introduced in Python 3.5. All the interfaces should be well documented.

## Installation

The `pdbeccdutils` can be presently obtained from [PYPI](https://pypi.org/project/pdbeccdutils/) using the following command:

```console
pip install pdbeccdutils
```

Alternativelly, you can install the reposotory from [Github](https://github.com/PDBeurope/ccdutils) using the following command:

```console
pip install git+https://github.com/PDBeurope/ccdutils.git@master#egg=pdbeccdutils
```

If you want to contribute to the project please fork it first and then do a pull request.

# Getting started

The centerpoint of the `pdbecccdutils` package is a `Component` object, which is a wrapper around the default `rdkit.Chem.rdchem.Mol` object (object property `mol`) providing most of the functionality and access to its properties. All the conformers are stored in the `rdkit.Chem.rdchem.Mol` with the exception of 2D depiction, as this one does not contain explicit hydrogens. `pdbeccdutils.core.modes.ConformerType` object allows accessing all of them.

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

### Bound-molecule from PDB model files

The bm_reader module infers all Bound-molecules (covalently-bound chemical components) within a single PDB mmCIF files. The result of
`bm_reader.read_pdb_updated_cif_file` function is a list of instances of `CCDReaderResult`, with each instance representing a single Bound-molecule. The `Component` Object of CCDReaderResult can then be used to probe the properties of each Bound-molecule.

```python
from pdbeccdutils.core import bm_reader

bms = bm_reader.read_pdb_updated_cif_file('/path/to/xxxx_updated.cif',sanitize=True)
bm_components = [bm.component for bm in bms]
rdkit_mols = [k.mol for k in bm_components]
```

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
