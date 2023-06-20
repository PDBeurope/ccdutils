```{eval-rst}
.. _intro:
```

# Introduction

`pdbeccdutils` is an open-source python package for processing and analyzing small molecules in PDB. Small-molecule data in PDB is available as [Chemical Component Dictionary (CCD)](http://www.wwpdb.org/data/ccd) or [Biologically Interesting Molecule reference Dictioanry (BIRD)](http://www.wwpdb.org/data/bird) in PDBX/mmCIF format. `pdbeccdutils` provides streamlined access to all metadata of small molecules in PDB and offers a set of convenient methods to compute various properties of small molecules using RDKIt such as 2D depictions, 3D conformers, physicochemical properties, matching common fragments and scaffolds, mapping to small-molecule databases using UniChem. `pdbeccdutils` also provides methods for identifying all the covalently attached chemical components in a macromolecular structure and calculating similarity among small molecules

**Note**
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

The core structural representation of small-molecules in `pdbecccdutils` package is a `Component` object, which is a wrapper around the default `rdkit.Chem.rdchem.Mol` object (object property `mol`) providing most of the functionality and access to its properties. Both `Ideal` and `Model` conformers are stored in the `mol` attribute and `Computed` and `Depiction` conformers are stroted in `mol3D` and `mol2D` attributes of `Component`. `pdbeccdutils.core.models.ConformerType` object allows accessing all of them.

Below you can find a few typical use cases.

## Reading CCD mmCIF files

CCD structures can be read using `ccd_reader.py` module located in the `pdbeccdutils.core` module. By default, the molecules comes sanitized using an augmented RDKit sanitization procedure. However, this option can be turned off by specifying optional parameter `sanitize=False` to the function

```python
from pdbeccdutils.core import ccd_reader

ccd_reader_result = ccd_reader.read_pdb_cif_file('HEM.cif')
ccd_reader_result
```
CCDReaderResult contains a list of possible warnings and errors that were encountered during the structure parsing. There is also a convenience method that allows reading in multiple chemical components provided they are listed in different data blocks in a single mmCIF file at the same time.

## Reading PRD mmCIF files

PRD structures can be read using `prd_reader.py` module located in `pdbeccdutils.core` module.

```python
from pdbeccdutils.core import prd_reader

prd_reader_result = prd_reader.read_pdb_cif_file('PRDCC_000204.cif')
prd_reader_result
```

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

## Infer Covalently Linked Components (CLC) from PDB model files

Several small and large ligands in PDB are split into individual CCDs. Such splitting of ligands to individual CCDs makes it difficult to correctly identify ligands in PDB and their interactions with macromolecules as well as mapping to other small-molecule databases. `pdbecccdutils` provides `clc_reader` module to infer all covalenty linked components in single PDB model file.

```python
from pdbeccdutils.core import clc_reader

clcs = clc_reader.read_pdb_updated_cif_file('/path/to/xxxx_updated.cif',sanitize=True)
clc_components = [clc.component for clc in clcs]
rdkit_mols = [k.mol for k in clc_components]
```

The result of `clc_reader.read_pdb_updated_cif_file` function is a list of instances of `CLCReaderResult`, with each instance representing a single Covalently Linked Components (CLC). The `Component` Object of `CLCReaderResult` can then be used to probe the properties of each CLC.


## Reading CLC mmCIF files

CLC structures can be read using `clc_reader.py` module located in `pdbeccdutils.core` module.

```python
from pdbeccdutils.core import clc_reader
clc_reader_result = clc_reader.read_pdb_cif_file('CLC_00004.cif')
clc_reader_result
```

## Writing CCD/PRD/CLC files

CCD/PRD/CLC molecules represented as Component objects in `pdbeccdutils` can be exported to different file formats such as mmCIF, SDF, PDB, CML, XML

```python
from pdbeccdutils.core import ccd_writer, prd_writer, clc_writer
from pdbeccdutils.core.models import ConformerType

ccd_component = ccd_reader_result.component
prd_component = prd_reader_result.component
clc_component = clc_reader_result.component

# write idealized coordinates in the SDF format.
ccd_writer.write_molecule('HEM.sdf', ccd_component)
prd_writer.write_molecule('PRDCC_000204.sdf', prd_component)
clc_writer.write_molecule('CLC_00004.sdf', clc_component)

# write model coordinates in the PDB format without hydrogens
ccd_writer.write_molecule('HEM.pdb', ccd_component, remove_hs=True, conf_type=ConformerType.Model)
prd_writer.write_molecule('PRDCC_000204.pdb', prd_component, remove_hs=True, conf_type=ConformerType.Model)
clc_writer.write_molecule('CLC_00004.pdb', clc_component, remove_hs=True, conf_type=ConformerType.Model)

# write model coordinates in the mmCIF format with hydrogens
ccd_writer.write_molecule('HEM.cif', ccd_component, conf_type=ConformerType.Model)
prd_writer.write_molecule('PRDCC_000204.cif', prd_component, conf_type=ConformerType.Model)
clc_writer.write_molecule('CLC_00004.cif', clc_component, conf_type=ConformerType.Model)
```
