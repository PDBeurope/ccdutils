# pdbeccdutils

* A set of python tools to deal with PDB chemical components definitions
  for small molecules, taken from the 
  [wwPDB Chemical Component Dictionary](https://www.wwpdb.org/data/ccd)
* Written to replace [PDBeChem](http://pdbe.org/chemistry/) back end
  processing.
  * the [wwPDB validation pipeline](https://www.wwpdb.org/validation/validation-reports)
* The tools use:
  * [RDKit](http://www.rdkit.org/) for chemistry
  * [PDBeCIF](https://github.com/glenveegee/PDBeCIF.git) cif parser.
* Please note that the project is under active development includes some rough 
  preliminary scripts and is and not yet ready for wider use!

## Installation instructions.

* `pdbeccdutils` requires RDKit to be installed.
  The official RDKit documentation has
  [installation instructions for a variety of platforms](http://www.rdkit.org/docs/Install.html).
  For linux/mac OS this is most easily done using the anaconda python with
  commands similar to:

  ```
  conda create -c rdkit -n rdkit-env rdkit python=3
  source activate rdkit-env
  ```
* `pdbeccdutils` also requires the [PDBeCIF](https://github.com/glenveegee/PDBeCIF.git) cif parser.
  Once you have installed RDKit, as described above then install PDBeCIF and pdbeccdutils using pip:

  ```
  pip install git+https://github.com/glenveegee/PDBeCIF.git 
  pip install https://gitlab.com/pdbe/ccd_utils/repository/latest_release/archive.zip
  ```
* Alternatively if you want to contribute to the project fork the repository and clone
  it. Then after getting rdkit setup, in the top `pdbeccdutils` directory:
  
  ```
  pip install -e ccdutils
  ```

## Documentation
 Note that `sphinx` needs to be a part of the virtual environmnent, 
 otherwise it cannot pick `rdkit` module. `sphinx_rtd_theme` is a theme providing
 nice `ReadtheDocs` mobile friendly style.
  
  * Generate *.rst files to be included as a part of the documentation. Inside the directory
  `pdbeccdutils/doc` run the following commands to generate documentation.
  
  Use the following to generate initial markup files to be used by sphinx.
  This needs to be used when other package but `core`, `utils` and `helpers` is implemented
  
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

## TODO list
  * Port rest of the important functionality implemented by Oliver
  * Port cofactors to this solution
  * Improve documentation








