# ccd_utils

* A set of python tools to deal with PDB chemical components definitions
  for small molecules, taken from the 
  [wwPDB Chemical Component Dictionary](https://www.wwpdb.org/data/ccd)
* Written for use in:
  * replacing [PDBeChem](http://www.ebi.ac.uk/pdbe-srv/pdbechem/) back end 
  processing.
  * the [wwPDB validation pipeline](https://www.wwpdb.org/validation/validation-reports)
* The tools use:
  * [RDKit](http://www.rdkit.org/) for chemistry
  * either the [PDBeCIF](https://github.com/glenveegee/PDBeCIF.git) 
  or the CifFile CIF parser used in wwPDB OneDep.
* Please note that the project is under active development includes some rough 
  preliminary scripts and is and not yet ready for wider use!

# current rather poor installation instructions:
* Set up rdkit with anaconda following http://www.rdkit.org/docs/Install.html
* Then issue commands like:
```
git clone https://github.com/glenveegee/PDBeCIF.git
git clone https://gitlab.com/pdbe/ccd_utils.git
pip install nose
cd ccd_utils
```








