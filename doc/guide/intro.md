```{eval-rst}
.. _intro:
```

# Introduction

`pdbeccdutils` is an open-source python package for processing and analyzing small molecules in PDB. Small-molecule data in PDB is available as [Chemical Component Dictionary (CCD)](http://www.wwpdb.org/data/ccd) or [Biologically Interesting Molecule reference Dictioanry (BIRD)](http://www.wwpdb.org/data/bird) in PDBX/mmCIF format. `pdbeccdutils` supports small molecule data in CCD/BIRD dictionaries and provide interoperability to RDKit. It also offers various features, such as generating 2D depictions, matching common fragments to specific molecules, identifying all the covalently attached chemical components in a macromolecule structure, and calcualting similarity between small molecules using Parity method.

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

Structure reading can be done using `ccd_reader.py` module located in the `pdbeccdutils.core` module. By default, the molecules comes sanitized using an augmented RDKit sanitization procedure. However, this option can be turned off by specifying optional parameter `sanitize=False` to the function

```python
from pdbeccdutils.core import ccd_reader

ccd_reader_result = ccd_reader.read_pdb_cif_file('HEM.cif')
ccd_reader_result
```
CCDReaderResult contains a list of possible warnings and errors that were encountered during the structure parsing. There is also a convenience method that allows reading in multiple chemical components provided they are listed in different data blocks in a single mmCIF file at the same time.

### Component

Component is a wrapper around `rdkit.Chem.rdchem.Mol` object providing streamlined access to all metadata information from CCD/BIRD files

```python
component = ccd_reader_result.component
component
```
```python
component.inchikey
```
```python
component.formula
```

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
