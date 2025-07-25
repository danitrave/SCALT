Inputs and Outputs
================

The input of SCARLET is a matrix of raw counts from single cell RNA sequencing. The matrix must be in **.tsv** extension reporting:

1. genes on the rows;
2. cells on the columns.

The first row of the matrix must contain the ids of each cell while the first column must provide the gene ids written either as **gene symbol** or **ensembl id** (not together). 

.. Note::

   In the case of enesembl gene id, the usage of the version of the id is not allowed. For example, **ENSG00000235430.1** cannot be used in the counts matrix. Use **ENSG00000235430** instead. 

An example of the input file is reported below.

.. list-table::  
   :widths: 50 50 50 50 50 50 50 50 50 50 50 50 50 50 50 50 50 50 50 50 50
   :header-rows: 1

   * - genes/cells
     - 1 
     - 2
     - 3
     - 4
     - 5
     - 6
     - 7
     - 8
     - 9
     - 10
     - 11
     - 12
     - 13
     - 14
     - 15
     - 16
     - 17
     - 18
     - 19
     - 20
   * - ENSG00000160072
     - 0
     - 0
     - 1
     - 5
     - 7
     - 0
     - 0
     - 1
     - 9
     - 0 
     - 10
     - 1
     - 0
     - 11
     - 2
     - 0
     - 0
     - 0
     - 1
     - 1
   * - ENSG00000142611
     - 4
     - 10
     - 0
     - 0
     - 0
     - 0
     - 1
     - 3
     - 1
     - 12
     - 10
     - 1
     - 0
     - 9
     - 7
     - 0
     - 0
     - 0
     - 0
     - 4
   * - ENSG00000157911
     - 0
     - 0
     - 0
     - 0
     - 0
     - 0
     - 0
     - 0
     - 0
     - 0 
     - 0
     - 0
     - 0
     - 0
     - 0
     - 0
     - 0
     - 0
     - 0
     - 0
   * - ENSG00000142655
     - 0
     - 0
     - 0
     - 0
     - 16
     - 14
     - 0
     - 0
     - 0
     - 0
     - 0
     - 0
     - 0
     - 20
     - 0
     - 0
     - 0
     - 0
     - 0
     - 0
   * - ENSG00000149527
     - 1
     - 1
     - 3
     - 2
     - 1
     - 2
     - 0
     - 1
     - 0
     - 5
     - 1
     - 1
     - 0
     - 2
     - 3
     - 1
     - 0
     - 5
     - 1
     - 1

The application returns three outputs:

1. a report file in **.html** format;
2. a tabular file in **.tsv** format summarizing the outcome of the classification per each cell;
3. a directory (named by default as **output**) hosting a collection results and metadata produced upon classification.

More details about the report and tabular outputs can be found in the **PBMC 3K tutorial**

The directory contains a series of files which are produced automatically during the classification step and are required for the generation of the report. Among them, we find:

1. **_adj_p_values.tsv**. Tabular file reporting a collection of p-values per each cell. Each p-value indicates to a likehood test. Therefore, the number corresponds to the number of cell types tested;
2. **_adj_deltas.tsv**. Tabular file reporting a series of likelihood differences between the cell type specific model and the mean cell type. The number corresponds to the number of cell types tested;
3. **_adj.tsv** file. Counts table after input set-up performed by SCARLET.py by default;
4. **_adj_genesExpressed_filter.tsv** file. Tabular file reporting either **PASS** or **EXCLUDE** if the cell expresses at least the minimum number of genes set in the -Min (or --Min) parameter or not;
5. **_barplot_cellTypesAboundance.html**, **_barplot_survivedCells.html**, **_UMAP_2D.html**, **_UMAP_2D_ONTO.html** and **_UMAP_3D.html**. Collection of plots visualized in the report file.
6. **_umap_2d_coords.tsv** and **_umap_3d_coords.tsv**. Tabular files reporting respectively the coordinates of the 2-dimensional and 3-dimensional UMAPs present in the report file.


SCARLET parameters
================

SCARLET.py makes usage of a collection of both positional arguments and parameters that can visualized typing the following command:

:: 

  python3 SCARLET.py -h

Or:

::

  python3 SCARLET --help

The documentation should appear as follows:

::

  usage: SCARLET.py [-h] [-Min --Min] [-Notation --Notation]
             [-Types --Types] [-CPUs --CPUs] [-pvalue --pvalue]
               [-out --out] Sample

1. **Sample** is the only positional argument of the tool. It represents the scRNA seq counts matrix file;
2. **-h** or **--help** shows the documentation.
3. **-Min** or **--Min** is the minimum number of genes that a cell must express to be classified. The **default** value is **250**;
4. **-Notation** or **--Notation** is the type of gene notation present in the counts. The defaul is **ensembl id**. Instead, write **gene_symbol** to switch to the gene symbol nomenclature;
5. **-Types** or **--Types** is the name of the directory containing the lists of the cell types to use in the likelihood test. By default, only the 293 pre-compiled lists (DISCO, HPA) are used. To use only the custom lists generated from annotation, insert **custom**;
6. **-CPUs** or **--CPUs** is number of processors employed. The default is **1**;
7. **-pvalue** or **--pvalue** indicates the significance level corresponding to the likelihood difference that there must be between the most signficant classification and the other pluasible ones in order to unequivocally annotate a cell to a type. If the likelihood difference is not reached, the cell will be classified as **multiassigned**. By default, the p-value threshold is **0.05** and corresponds to a likelihood difference of 6. Set the threshold to **0.01** to increase the stringency of the likelihood difference up to 9;
8. **-out** or **--out** is the name referring to the final output files. By default is **output**.

Run SCARLET.py
=========

SCARLET.py is quite straightforward since it requires just the counts table as positional input. 

Leaving default parameters, the basic command appears as follows:

::

   python3 SCARLET.py read_counts.tsv

By default, the ensembl ids are used. 

If the **gene symbol** is used in the counts matrix, the notation must be specified as follows:

::

   python3 SCARLET.py read_counts.tsv -Notation gene_symbol

Or:

::

   python3 SCARLET.py read_counts.tsv --Notation gene_symbol

By default, a cell is classified if it expresses at least **250** genes. Managing the SCARLET.py parameters, this threshold che be modified with any number as follows:

::

   python3 SCARLET.py read_counts.tsv -Min 500

Or:

::

   python3 SCARLET.py read_counts.tsv --Min 500

In addition, the computational time can be reduced if the number of processors is increased as reported:

::

   python3 SCARLET.py read_counts.tsv -CPUs 4

Or:

::

   python3 SCARLET.py read_counts.tsv --CPUs 4

Make sure to have available the number of desidered processors on your machine.

The significance threshold can be modified in the following way:

::

   python3 SCARLET.py read_counts.tsv -pvalue 0.01

Or:

::

   python3 SCARLET.py read_counts.tsv --pvalue 0.01

Finally, the name present in the output files can be changed as follows:

::

   python3 SCARLET.py read_counts.tsv -out my_output

Or:

::

   python3 SCARLET.py read_counts.tsv --out my_output

Adjusting the parameters in a unique call, the final command should appear as follows:

::

   python3 SCARLET.py read_counts.tsv -Notation gene_symbol -Min 500 -CPUs 4 -pvalue 0.01 -out my_output

Or:

::

   python3 SCARLET.py read_counts.tsv --Notation gene_symbol --Threshold 500 --CPUs 4 --pvalue 0.01 --out my_output

The order of parameters is irrelevant.


Report
======

The report is a file in html format composed of a collection of plots summarizing the general statistics and classification results of the analysis. The file reports four different plots:

1. a bar plot showing how many cells express or not the minimum number of genes for classification;
2. a second barplot counting how many cells were classified to each cell type cathegory;
3. a 2D UMAP;
4. a 3D UMAP.
5. a 2D UMAP where cells are colored based on the cell ontology.

.. note::
   The genes used for the creation of the UMAPs coordinates are the union of genes deriving from the 293 cell types that managed to annotate at least 50 cells without repetitions.

Workflow 
========

Running SCARLET.py, the following workflow is performed:

.. figure:: pictures/SCARLET_workflow.png
   :align: center
   :scale: 50%

1. **inputPreparation.py** is a python script that adjustes the input counts table in order to be properly analyzed by SCARLET.py;
2. **likelihood_ratio_test.py** is the python script that performs the actual likelihood test;
3. **reportGenerator.py** is the python script that creates the final report.


