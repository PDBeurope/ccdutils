# ccd_utils

* A set of python tools to deal with PDB chemical components definitions
  for small molecules, taken from the 
  [wwPDB Chemical Component Dictionary](https://www.wwpdb.org/data/ccd)
* Written to replace [PDBeChem](http://pdbe.org/chemistry/) back end
  processing.
  * the [wwPDB validation pipeline](https://www.wwpdb.org/validation/validation-reports)
* The tools use:
  * [RDKit](http://www.rdkit.org/) for chemistry
  * for cif parsing:
    * either the [PDBeCIF](https://github.com/glenveegee/PDBeCIF.git) cif parser.
    * or the CifFile CIF parser used in wwPDB OneDep (not extensively tested).
* Please note that the project is under active development includes some rough 
  preliminary scripts and is and not yet ready for wider use!

## installation instructions.

* `ccd_utils` requires RDKit to be installed.
  The official RDKit documentation has
  [installation instructions for a variety of platforms](http://www.rdkit.org/docs/Install.html).
  For linux this is most easily done using the anaconda python with
  commands similar to:

  ```
  conda create -c rdkit -n my-rdkit-env rdkit python=3
  source activate my-rdkit-env
  ```
* `ccd_utils` also requires the [PDBeCIF](https://github.com/glenveegee/PDBeCIF.git) cif parser.
  Once you have installed RDKit, as described above then install PDBeCIF and ccd_utils using pip:

  ```
  pip install https://github.com/glenveegee/PDBeCIF/zipball/master 
  pip install https://gitlab.com/pdbe/ccd_utils/repository/latest_release/archive.zip
  ```
* Alternatively if you want to contribute to the project fork the repository and clone
  it. Then after getting rdkit setup, in the top `ccd_utils` directory:
  
  ```
  pip install https://github.com/glenveegee/PDBeCIF/zipball/master
  pip install -e . 
  ```
  






