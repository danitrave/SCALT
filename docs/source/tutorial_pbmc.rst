PBMC 3K tutorial
================

The tutorial consist in analyzing a sample of perifieral blood mononuclear cells (PBMCs) that can be downloaded from the Seurat vignette available at the following `link <https://satijalab.org/seurat/articles/pbmc3k_tutorial>`_. 
Make sure to convert counts from a 10X to a tsv format.

Once obtained a file in tsv format, the only command needed to run SCALT is the following:

::

   python3 SCALT.py countsPBMC3K.tsv -Notation ensembl_id  

The output consists in a series of files and metadata table that should provide a complete picture of the composition of cells present in the sample according to SCALT.
The first file is a **report** in html format comprehensive of series of pictures reported below.

.. figure:: pictures/countsPBMC_adj_UMAP_2D.html
   :align: center
   :scale: 50%
