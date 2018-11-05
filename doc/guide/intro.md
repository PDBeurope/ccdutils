```eval_rst
.. _intro:
```

# Introduction

The `pdbeccdutils` is presently under development and new functionality is added regularly as well as its functionality is being revised and updated. When properly installed all the code should have documentation. All the *public* methods do have [static typing](http://mypy-lang.org/) introduced in Python 3.5. All the interfaces should be well documented.

The centerpoint of the `pdbecccdutils` package is a `Component` object, which is a wrapper around the default `rdkit.Chem.rdchem.Mol` object providing most of the functionality and access to its properties. All the conformers are stored in the `rdkit.Chem.rdchem.Mol` with the exception of 2D depiction, as this one does not contain hydrogen atoms. `pdbeccdutils.core.modes.ConformerType` object allows accessing all of them.

Below you can find a few typical use cases.

## Reading CCD mmCIF files

```python
from pdbeccdutils.core import ccd_reader

component = ccd_reader.read_pdb_cif_file('HEM.cif').component
rdkit_mol = component.mol
```

The `rdkit_mol` object is sanitized already!

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

## Compute 2D depictions

The depictions rely on two optional external resources. [PubChem layouts](https://www.ncbi.nlm.nih.gov/pccompound) and general templates for substructures. e.g. [porphyrin rings](https://en.wikipedia.org/wiki/Porphyrin). These can be supplied as a paths to folders with 2D layouts for molecules in an SDF format. If neither of the paths is supplied, `pdbeccdutils` uses a few templates provided with the source code.

```python
from pdbeccdutils.core.depictions import DepictionManager

depictions = DepictionManager('pubchem_templates_dir', 'general_templates_dir')
depiction_result = component.compute_2d(depictions)

rdkit_2d_mol = depiction_result.mol
# Usefull to find out whether or not the drawing is planar or contains clashes.
# Ideally should be 0, but everything < 1 is OK.
score = depiction_result.score

component.export_2d_svg('HEM.svg')
component.export_2d_svg('HEM_names.svg', names=True)
```

---
<div align='center'>
    <img src='../_static/HEM_300.svg' />
    <img src='../_static/HEM_300_names.svg' />
</div>
