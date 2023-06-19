# virtualEmbryo
This repository contains the scripts used during the virtualEmbryo project. The virtualEmbryo refers to [...]

A visual representation of the virtualEmbryo can be found under the Drosophila Virtual MUlti-modal eXplorer : [DVMuX](https://shiny.mdc-berlin.de/DVMuX/)

## Overview scripts
| Name    | Description | 
| -------- | ------- | 
|  annotate_ATAC_peaks.Rmd | annotating ATAC-peaks as putative enhancer,promoters or CDS  |
|  create_archetypes.Rmd |  apply multivariate factor analysis (NMF) on spatially reconstructed gene expression patterns and chromatin accessibilites  | 
|  enrich6.Rmd |  motif enrichment analysis | 
|  functions_regulation.R |  collection of functions used throughout the analysis  |
|  p2_stage6_withmarkers.ipynb | spatial reconstruction of the sequenced cells with [novoSpaRc](https://github.com/rajewsky-lab/novosparc)   |
|  predictions.Rmd |  using the glmnet package to predict enhancers |
|  sample_AB2.Rmd |  from sequencing data to count matrices. Clustering of the cells.  |
|  virtual_tfbs.Rmd | computing the response accessibilities for a specific TFBS  | 

