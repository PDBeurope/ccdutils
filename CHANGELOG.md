# Changelog

## RELEASE 1.0.0 - Sep 4, 2025

### Features
* Support for RDKit 2025 and NumPy 2.0
* Uses both coordgen library and the default RDKit functionality to generate 2D depictions

## RELEASE 0.8.6 - Oct 28, 2024

### Features
* Enable parsing of a subset of CCDs from the Chemical Component Dictionary
* Added CCDC to UniChem resources


## RELEASE 0.8.5 - May 26, 2024

### Features
* Updated for RDKit 2023.09.6
* Removed the formal charge adjustment for valence issues, instead uses Dative bonds
* Uses Zero bond type (dotted lines) for depiction of coordinate bonds
* Changes Dative bonds to single bonds for inchi calculation
* Corrected the format of PDB files
* Corrected the header of SDF file

### Breaking changes
* Removed smile from SubstructureMapping and added mol

## RELEASE 0.8.0 - July 26, 2023

### Features
* Added clc_reader and clc_writer to read and write Covalently Linked Components (CLCs)
* Added prd_reader and prd_writer to read and write Protein Reference Dictionarys (PRDs)
* Added boundmolecule_cli, a pipeline to infer CLCs from PDB model files

## RELEASE 0.7.0 - February 26, 2023

### Features
* Replaced pdbecif with gemmi for mmcif parsing
* Support for rdkit 2022.09.x
* Component representation of Covalently Linked Components (complex multi-component ligands typically represented as individual components represented by individual CCDs)

### Breaking changes
* Removed ccd_cif_dict (dict) property of Component. Please use ccd_cif_block (gemmi.cif.Block) to access data from mmcif file

## RELEASE 0.6 - April 26, 2021

### Features

* A lot of minor bug fixes and code improvements
* New templates for 2D layouts
* Support for rdkit 2021.03.x

## RELEASE 0.5 - May 15, 2019

### Features

* Add support for UniChem mapping.
* Add bond information to the SVG decomposition.
* Allow PARITY method to be atom/bond specific.
* Improve and extend physchem properties.
* Enhanced CIF export (physchem, scaffolds, fragments, 2D, mapping).
* Few improvements to match newest RDKit version (**breaking changes**).

## RELEASE 0.4 - January 12, 2019

### Features

* Add SVG decomposition in the SVG format.
* Protein-ligand interaction pipeline moved to separate [repository](https://gitlab.ebi.ac.uk/pdbe/release/interactions).
* Add basic properties calculation (Abhik).
* Extension and improvements of the Fragment library.

## RELEASE 0.3 - October 12, 2018

### Features

* Introduce CoordGen from RDKit.
* Refactored pdbeccdutils.core (**breaking changes**).
* Scaffolding (Abhik).
* Add protein-ligand interaction pipeline.

## RELEASE 0.2 - June 14, 2018

### Features

* PDBeChem pipeline.
* Support for PARITY method.
* Introduce EKTGv2() method for 3D conformer generation.
