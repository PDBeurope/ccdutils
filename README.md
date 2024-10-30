[![CodeFactor](https://www.codefactor.io/repository/github/PDBeurope/ccdutils/badge/master)](https://www.codefactor.io/repository/github/PDBeurope/ccdutils/overview/master)  ![PYPi](https://img.shields.io/pypi/v/pdbeccdutils?color=green&style=flat)  ![GitHub](https://img.shields.io/github/license/PDBeurope/ccdutils)   ![ccdutils documentation](https://github.com/PDBeurope/ccdutils/workflows/ccdutils%20documentation/badge.svg) ![ccdutils tests](https://github.com/PDBeurope/ccdutils/workflows/ccdutils%20tests/badge.svg) ![PyPI Downloads](https://img.shields.io/pypi/dm/pdbeccdutils)


# pdbeccdutils

An RDKit-based python toolkit for parsing and processing small molecule definitions in [wwPDB Chemical Component Dictionary](https://www.wwpdb.org/data/ccd) and [wwPDB The Biologically Interesting Molecule Reference Dictionary](https://www.wwpdb.org/data/bird).`pdbeccdutils` provides streamlined access to all metadata of small molecules in the PDB and offers a set of convenient methods to compute various properties of small molecules using RDKIt such as 2D depictions, 3D conformers, physicochemical properties, matching common fragments and scaffolds, mapping to small-molecule databases using UniChem.

## Features

* `gemmi` CCD read/write.
* Generation of 2D depictions (`No image available` generated if the flattening cannot be done) along with the quality check.
* Generation of 3D conformations.
* Fragment library search (PDBe hand-curated library, ENAMINE, DSI).
* Chemical scaffolds (Murcko scaffold, Murcko general, BRICS).
* Lightweight implementation of [parity method](https://doi.org/10.1016/j.str.2018.02.009) by Jon Tyzack.
* RDKit molecular properties per component.
* UniChem mapping.
* Generating complete representation of multiple [Covalently Linked Components (CLC)](https://www.ebi.ac.uk/pdbe/news/introducing-covalently-linked-components)

## Dependencies

  * [RDKit](http://www.rdkit.org/) for small molecule representation. Presently tested with `2023.9.6`
  * [GEMMI](https://gemmi.readthedocs.io/en/latest/index.html) for parsing mmCIF files.
  * [scipy](https://www.scipy.org/) for depiction quality check.
  * [numpy](https://www.numpy.org/) for molecular scaling.
  * [networkx](https://networkx.org/) for bound-molecules.


## Installation

create a [virtual environment](https://packaging.python.org/en/latest/guides/installing-using-pip-and-virtual-environments/#create-and-use-virtual-environments) and install using pip

  ```bash
  pip install pdbeccdutils
  ```

## Contribution
We encourage you to contribute to this project. The package uses [poetry](https://python-poetry.org/) for packaging and dependency management. You can develop locally using:

```bash
git clone https://github.com/PDBeurope/ccdutils.git
cd ccdutils
pip install poetry
poetry install --with tests,docs
pre-commit install
```

The pre-commit hook will run linting, formatting and update `poetry.lock`. The `poetry.lock` file will lock all dependencies and ensure that they match pyproject.toml versions.

To add a new dependency

```bash
# Latest resolvable version
poetry add <package>

# Optionally fix a version
poetry add <package>@<version>
```

To change a version of a dependency, either edit pyproject.toml and run:

```bash
poetry sync --with dev
```

or

```bash
poetry add <package>@<version>
```


## Documentation

The documentation is generated using `sphinx` in `sphinx_rtd_theme` and hosted on GitHub Pages. To generate the documentation locally,

```bash
cd doc
poetry run sphinx-build -b html . _build/html

# See the documentation at http://localhost:8080.
python -m http.server 8080 -d _build/html
```
