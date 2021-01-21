[![CodeFactor](https://www.codefactor.io/repository/github/pdbeurope/ccdutils/badge/master)](https://www.codefactor.io/repository/github/pdbeurope/ccdutils/overview/master)  ![PYPi](https://img.shields.io/pypi/v/pdbeccdutils?color=green&style=flat)  ![GitHub](https://img.shields.io/github/license/pdbeurope/ccdutils)   ![ccdutils documentation](https://github.com/PDBeurope/ccdutils/workflows/ccdutils%20documentation/badge.svg) ![ccdutils tests](https://github.com/PDBeurope/ccdutils/workflows/ccdutils%20tests/badge.svg)

# pdbeccdutils

* A set of python tools to deal with PDB chemical components definitions
  for small molecules, taken from the [wwPDB Chemical Component Dictionary](https://www.wwpdb.org/data/ccd)

* The tools use:
  * [RDKit](http://www.rdkit.org/) for chemistry
  * [PDBeCIF](https://github.com/PDBeurope/pdbecif) cif parser.
  * [scipy](https://www.scipy.org/) for depiction quality check.
  * [numpy](https://www.numpy.org/) for molecular scaling.

* Please note that the project is under active development.

## Installation instructions

* `pdbeccdutils` requires RDKit to be installed.
  The official RDKit documentation has [installation instructions for a variety of platforms](http://www.rdkit.org/docs/Install.html).
  For linux/mac OS this is most easily done using the anaconda python with commands similar to:

  ```console
  conda create -c conda-forge -n rdkit-env rdkit python=3.7
  conda activate rdkit-env
  ```

* Once you have installed RDKit, as described above then install pdbeccdutils using pip:

  ```console
  pip install pdbeccdutils
  ```

## Features

* mmCIF CCD read/write.
* Generation of 2D depictions (`No image available` generated if the flattening cannot be done) along with the quality check.
* Generation of 3D conformations.
* Fragment library search.
* Chemical scaffolds (Murcko scaffold, Murcko general, BRICS).
* Lightweight implementation of [parity method](https://doi.org/10.1016/j.str.2018.02.009) by Jon Tyczak.
* RDKit molecular properties per component.
* UniChem mapping.

## TODO list

* Add more unit/regression tests to get higher code coverage.
* Further improvements of the documentation.

## Notes

* Protein-ligand interaction has been extracted [here](https://gitlab.ebi.ac.uk/pdbe/release/interactions). This was because of the fact that at the end of the day it was not using any of the pdbeccdutils functionality and introduced additional dependencies on the package.

## Documentation

The documentation depends on the following packages:

* `sphinx`
* `sphinx_rtd_theme`
* `recommonmark`
* `sphinx-autodoc-typehints`

Note that `sphinx` needs to be a part of the virtual environment, if you want to generate documentation by yourself.
Otherwise it cannot pick `rdkit` module. `sphinx_rtd_theme` is a theme providing nice `ReadtheDocs` mobile friendly style.

* Generate *.rst* files to be included as a part of the documentation. Inside the directory `pdbeccdutils/doc` run the following commands to generate documentation.
* Alternatively, use the `recommonmark` package along with the proper configuration to get the Markdown working.

 Use the following to generate initial markup files to be used by sphinx.  This needs to be used when adding another sub-packages.

```console
sphinx-apidoc -f -o /path/to/output/dir ../pdbeccdutils/
```

Use this to re-generate the documentation from the doc/ directory:

```console
make html
```
