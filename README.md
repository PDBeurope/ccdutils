[![pipeline status](https://gitlab.ebi.ac.uk/pdbe/ccdutils/badges/master/pipeline.svg)](https://gitlab.ebi.ac.uk/pdbe/ccdutils/commits/master)
[![coverage report](https://gitlab.ebi.ac.uk/pdbe/ccdutils/badges/master/coverage.svg)](https://gitlab.ebi.ac.uk/pdbe/ccdutils/commits/master)


# pdbeccdutils
* A set of python tools to deal with PDB chemical components definitions
  for small molecules, taken from the 
  [wwPDB Chemical Component Dictionary](https://www.wwpdb.org/data/ccd)
* Written to replace [PDBeChem](http://pdbe.org/chemistry/) back end
  processing.
  * the [wwPDB validation pipeline](https://www.wwpdb.org/validation/validation-reports)
* The tools use:
  * [RDKit](http://www.rdkit.org/) for chemistry
  * [PDBeCIF](https://gitlab.com/glenveegee/PDBeCIF.git) cif parser.
* Please note that the project is under active development includes some rough 
  preliminary scripts and is and not yet ready for wider use!


## Installation instructions.
* `pdbeccdutils` requires RDKit to be installed.
  The official RDKit documentation has
  [installation instructions for a variety of platforms](http://www.rdkit.org/docs/Install.html).
  For linux/mac OS this is most easily done using the anaconda python with
  commands similar to:

  ```console
  conda create -c rdkit -n rdkit-env rdkit python=3
  source activate rdkit-env
  ```
* Once you have installed RDKit, as described above then install pdbeccdutils using pip:

  ```console
  pip install git+https://gitlab.ebi.ac.uk/pdbe/ccdutils.git
  ```

## Documentation
The documentation depends on the following packages:
  * `sphinx`
  * `sphinx_rtd_theme`
  * `recommonmark`
  * `sphinx-autodoc-typehints`

Note that `sphinx` needs to be a part of the virtual environmnent, if you want to generate documentation by yourself. 
Otherwise it cannot pick `rdkit` module. `sphinx_rtd_theme` is a theme providing nice `ReadtheDocs` mobile friendly style.

  * Generate *.rst files to be included as a part of the documentation. Inside the directory `pdbeccdutils/doc` run the following commands to generate documentation.
  * Alternativelly, use the `recommonmark` package along with the proper configuration to get the Markdown working.
  
  Use the following to generate initial markup files to be used by sphinx.
  This needs to be used when adding another subpackages.
  
  ```
  sphinx-apidoc -f -o /path/to/output/dir ../pdbeccdutils/
  ```

  Use this to re-generate the documentation from the doc/ directory:
  ```
  make html
  ```


## Features
  * Generation of 2D depictions (`No image available` generated if the flattening cannot be done) along with the quality check.
  * Generation of 3D conformations.
  * Fragment library search.
  * Lightweight implementation of parity method by Jon Tyczak.

## TODO list
  * Port rest of the important functionality implemented by Oliver
  * Add more unit/regression tests to get at least 100% code coverage.
  * Add protein-ligand environment generation  
  * Further improvement of the documentation
  * Scafolding (Abhik)








