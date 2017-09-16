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









