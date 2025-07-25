# SCARLET: automatic identification of cell types from single-cell RNA sequencing data
SCARLET (Single Cell Annotation Likelihood Tool) is an innovative method which introduces a paradigm-shift for the analysis of scRNAseq data. In this approach, cells are annotated to a specific type at individual level, by using a simple but elegant method based on maximum likelihood, without the need for clustering, dimensionality reduction or manual annotation. SCARLET leverages a collection of 293 lists of cell-type specific genes, constructed by extensive re-analysis of comprehensive and expert curated catalogues (HPA and DISCO).

For further details, we advise to read the paper available at DOI

As mentioned in the paper, SCARLET is composed of other parallel utilities that allow to perform more complex analysis such as:

1. build cell type specific lists of genes in a deterministic fashion starting from a counts matrix and the correspoding annotation for each cell;
2. build cell type specific lists of genes starting from a count matrix and a gathering of user-defined cell type specific lists of genes making the use of an hypergeometric test.


## SCARLET manual
The manual of SCARLET is available at [SCARLET documentation](https://SCALT.readthedocs.io/en/latest/)
