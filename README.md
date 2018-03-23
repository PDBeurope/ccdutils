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
  conda create -c rdkit -n my-rdkit-env rdkit python=3
  source activate my-rdkit-env
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
  pip install https://github.com/glenveegee/PDBeCIF/zipball/master
  pip install -e . 
  ```
  
  ## Features
  * Generation of 2D depictions (`No image available` generated if the flattening cannot be done) along with the quality check.
  * Generation of 3D conformations.
  * Fragment search

  ## TODO list
  * Port rest of the important functionality implemented by Oliver
  * write mmcif Compound exporter (Lukas)
  * Once pip 9.1 is released change the way how pdbecif parser is installed (https://stackoverflow.com/questions/15221473/how-do-i-update-pip-itself-from-inside-my-virtual-environment)







