```{eval-rst}
.. _pipelines:
```

# Pipelines

`pdbeccdutils` can be used to process standard ligand definitions defined by [CCD](http://www.wwpdb.org/data/ccd) or [BIRD](https://www.wwpdb.org/data/bird) dictionaries, identify boundmolecules from PDB model files and infer Covalently Linked Components (CLC). The packages comes with two pipelines, one to process CCD/PRDs and an another one to process boundmolecules from PDB model files

## Processing of CCD/PRD

`process_components_cif` entry point can be used to run the pipeline to process a CCD or PRD. The pipeline enriches (data-enrichment-process) the information in standard wwPDB CCD/PRD files, generates 2D images and exports to multiple file formats.

To find all the input options, use:

```bash
process_components_cif -h
```
For example, following files are generated for ATP

* `ATP.cif` - Standard wwPDB CCD file with the [data enrichments](data-enrichment-process) | `mmCIF`.
* `ATP_ideal.pdb` - ideal cooordinates | `PDB`.
* `ATP_ideal_alt.pdb` - ideal cooordinates with atom alternate names | `PDB`.
* `ATP_model.pdb` - model cooordinates | `PDB`.
* `ATP_model_alt.pdb` - model cooordinates with atom alternate names | `PDB`.
* `ATP_N.svg` - 2D depiction in N x N resolution. Where N is (100,200,300,400,500) pixels.
* `ATP_N_names.svg` - 2D depiction in N x N resolution with atom names. Where N in (100,200,300,400,500) pixels.
* `ATP_model.sdf` - model coordinates | `MOL`
* `ATP_ideal.sdf` - ideal cooordinates | `MOL`.
* `ATP.cml` - model coordinates | `CML`.
* `ATP_annotation.json` - 2D depiction in JSON format with some additional annotation.


(data-enrichment-process)=
### Data enrichment process

Aside from the standard content of the wwPDB CCD files, `process_components_cif` pipeline adds the following information in the `mmCIF` files.

* 2D coordinations of the ligand.
* RDKit regenerated 3D conformer.
* Information about scaffold and fragments found in the entry.


## Processing of boundmolecules from PDB model files

Most of the structures in PDB has atleast one small molecule bound to it. Many of these small molecules are ligands of biological significance, such as cofactors, metabolites, carbohydrates, or lipids. `read_boundmolecule` entry point can be used to run the pipeline to identify these boundmolecules from PDB entries and assign them as Covalently Linked Components (CLC) if they are composed of individual components represented by individual CCDs.

To find all the input options, use:

```bash
read_boundmolecule -h
```
For example, the pipeline generates the following files in PDB entry 1d83

* `1d83_processed.cif` - The input structure of 1d83 after removing alternate conformers
* `bound_molecules.json` - Details of boundmolecules found in the input structure
* `CLC_1` - Folder containing details of CLC identified from input structure
* `CLC_2` - Folder containing details of CLC identified from input structure

In each folder of CLCs identified from an input structure, the following files are stored:

* `CLC_K.cif` - CIF files of CLC after [data enrichments](data-enrichment-process) | `mmCIF`.
* `CLC_K_model.pdb` - model cooordinates | `PDB`.
* `CLC_K_N.svg` - 2D depiction in N x N resolution. Where N is (100,200,300,400,500) pixels.
* `CLC_K_N_names.svg` - 2D depiction in N x N resolution with atom names. Where N in (100,200,300,400,500) pixels.
* `CLC_K_model.sdf` - model coordinates | `MOL`
* `CLC_K.cml` - component representation | `CML`.
* `CLC_K_annotation.json` - 2D depiction in 'natural format' (i.e. 50px per 1Å) with some additional annotation.

Where K is a positive interger representing the K<sup>th</sup> CLC identified from an input structure
