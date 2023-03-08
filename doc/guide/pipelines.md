```{eval-rst}
.. _pipelines:
```

# Pipelines

Presently there are four pipelines in the production based on the `pdbeccdutils` code.

## PDBeChem

This pipeline generates all the chemistry data consumed on the PDBe pages. It can process standard ligand definitions defined by [CCD](http://www.wwpdb.org/data/ccd) or [BIRD](https://www.wwpdb.org/data/bird) dictionaries.

### FTP area content

Following files are generated in the tree-like structure (A/ATP/...) in our [FTP area](http://ftp.ebi.ac.uk/pub/databases/msd/pdbechem_v2/).

* `ATP.cif` - Standard wwPDB CCD file with the [data enrichments](data-enrichment-process) | `mmCIF`.
* `ATP_ideal.pdb` - ideal cooordinates | `PDB`.
* `ATP_ideal_alt.pdb` - ideal cooordinates with atom alternate names | `PDB`.
* `ATP_model.pdb` - model cooordinates | `PDB`.
* `ATP_model_alt.pdb` - model cooordinates with atom alternate names | `PDB`.
* `ATP_N.svg` - 2D depiction in N x N resolution. Where N is (100,200,300,400,500) pixels.
* `ATP_N_names.svg` - 2D depiction in N x N resolution with atom names. Where N in (100,200,300,400,500) pixels.
* `ATP_model.sdf` - model coordinates | `MOL`
* `ATP_ideal.sdf` - ideal cooordinates | `MOL`.
* `ATP.cml` - component representation | `CML`.
* `ATP_annotation.json` - 2D depiction in 'natural format' (i.e. 50px per 1Å) with some additional annotation. This file is consumed by the protein-ligand interaction viewer.

On the top of that these files are generated:

* `chem_comp.list` - List of processed ids
* `chem_comp_list.xml` - XML file with some additional metadata per entry.
* `components.cif` - Aggregated files with all the CCD information.
* `pdbechem.tar.gz` - The archive in the `*.tar.gz` format:

(data-enrichment-process)=
### Data enrichment process

Aside from the standard content of the wwPDB CCD files. PDBeChem adds the following information in the `mmCIF` files. See example for [ATP](http://ftp.ebi.ac.uk/pub/databases/msd/pdbechem_v2/A/ATP/ATP.cif):

* ID mapping to some popular resources like ChEMBL, ChEBI, DrugBank, etc. The curated list is available in the pdbeccdutils package.(`pdbeccdutils.utils.web_services.agreed_resources`).
* ChEMBL synonyms.
* DrugBank details - summary, synonyms, international brands, taxonomy, and known targets.
* 2D coordinations of the ligand.
* RDKit regenerated 3D conformer.
* Information about scaffold and fragments found in the entry.


## [Cofactors](https://pdbe.gitdocs.ebi.ac.uk/release/relic/documentation/cofactors.html)


## [Reactants](https://pdbe.gitdocs.ebi.ac.uk/release/relic/documentation/reactants.html)

## [Similarities](https://pdbe.gitdocs.ebi.ac.uk/release/relic/documentation/similarities.html)
