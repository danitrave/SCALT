Inputs formats
==============

SCALT_AnnotaionListsBuilder.py demands two inputs:

1. a scRNA seq row counts matrix. The matrix must be in .**tsv** extension;
2. a table having **one column** in which the annotation of each cell present in the counts is reported.

The counts matrix must present genes on the rows and cells on the columns. The first row of the matrix must contain the ids of each cell while the first column must provide the gene ids written either as **gene symbol** or **ensembl id**. 

To see an example of the table, see the section **SCALT: CLASSIFICATION - Inputs & Outputs** of this manual.

An example of the annotation table is reported here:

.. list-table:: 
   :align: center
   :widths: 80 
   :header-rows: 1

   * - annotation
   * - acinar
   * - beta
   * - beta
   * - alpha
   * - alpha
   * - alpha

Parameters
==========

SCALT_AnnotaionListsBuilder.py takes the advantage of a collection of both positional arguments and parameters that can visualized typing the following command:

:: 

  python3 SCALT_AnnotaionListsBuilder.py -h

Or:

:: 

  python3 SCALT_AnnotaionListsBuilder.py --help

The documentation should appear as follows:

::

   usage: SCALT_AnnotaionListsBuilder.py [-h] [-Boo --Boostraps] [-Cells --Cells]
                                      [-Genes --Genes] [-Notation --Notation]
                                      [-CPUs --CPUs]
                                      Sample Annotation


1. **Sample** is the first positional argument and addresses the name of the counts matrix;
2. **Annotation** is the second positional argument and refers to the name of the annotation table;
3. **-Boo** or **--Boostraps** indicates the number of boostrap samples to generate. The default number is **100**;
4. **-Cells** or **--Cells** specifies the number of cells to pick randomly per each cell type during the probability inference process. By default, the number is set to **100**;
5. **-Genes** or **--Genes** refers to the number of genes that the final cell type lists of genes must contain at the end. The default number is **100**;
6. **-Notation** is used to underline the kind of gene notation present in the counts. The user can choose between **gene_symbol** or **ensembl_id**. By default, ensembl_id is set;
7. **-CPUs** or **--CPUs** indicates the number of processors to use. By default, **1** is used.


Run SCALT_AnnotaionListsBuilder.py
==================================

SCALT_AnnotaionListsBuilder.py is quite straightforward since it requires the counts table and the annotation table as inputs. 

Leaving default parameters, the basic comand appears as follows:

::

   python3 SCALT_AnnotaionListsBuilder.py read_counts.tsv annotation.tsv

By default, it is assumed that the counts have the ensembl id for genes notation. 100 boostrap samples are generated in which 100 cells will be randomly sampled per each cell type. At the end, the final lists will contain 100 genes each. The number of processors used is 1.

If the **gene symbol** is used in the counts matrix, the notation must be specified as follows:

::

   python3 SCALT_AnnotaionListsBuilder.py read_counts.tsv annotation.tsv -Notation gene_symbol

Or:

::

   python3 SCALT_AnnotaionListsBuilder.py read_counts.tsv annotation.tsv --Notation gene_symbol

The number of boostrap samples can be changed in any moment specifying the number in the proper parameter:

::

   python3 SCALT_AnnotaionListsBuilder.py read_counts.tsv annotation.tsv -Boo 80

Or:

::

   python3 SCALT_AnnotaionListsBuilder.py read_counts.tsv annotation.tsv --Boostraps 80

The number of cell per cell type included in each boostrap sample can be modified making usage of the proper parameter:

::

   python3 SCALT_AnnotaionListsBuilder.py read_counts.tsv annotation.tsv -Cells 50

Or:

::

   python3 SCALT_AnnotaionListsBuilder.py read_counts.tsv annotation.tsv --Cells 50

The same logic for the number of genes that will be included in the final lists of genes:

::

   python3 SCALT_AnnotaionListsBuilder.py read_counts.tsv annotation.tsv -Genes 75

Or:

::

   python3 SCALT_AnnotaionListsBuilder.py read_counts.tsv annotation.tsv --Genes 75

In addition, the computational time can be reduced if the number of processors is increased as reported:

::

   python3 SCALT_AnnotaionListsBuilder.py read_counts.tsv annotation.tsv -CPUs 4

Or:

::

   python3 SCALT_AnnotaionListsBuilder.py read_counts.tsv annotation.tsv --CPUs 4

Make sure to have available the number of desidered processors in your machine.

To conclude, the different parameters can be modified in a unique call:

::

   python3 SCALT_AnnotaionListsBuilder.py read_counts.tsv annotation.tsv -Notation gene_symbol -Boo 80 -Cells 50 -Genes 75 -CPUs 4

Or:

::

   python3 SCALT_AnnotaionListsBuilder.py read_counts.tsv annotation.tsv --Notation gene_symbol --Boostraps 80 --Cells 50 --Genes 75 --CPUs 4

The order of parameters is irrelevant.

Outputs
=======

The tool returns two output:

1. a directory called **custom** containing the final lists of genes;
2. a directory named **AnnolistsBuilder_results** hosting a collection of metadata.

The metadata consists in a series of files and directories which are produced automatically during the process and were utilized for the generation of the final lists:
  
1. **originalTables_zipped.zip** is a zipped repository containing the original input data;
2. **groupped_cell_types** is the directory that contains the counts matrix groupped by cell type. Each tsv file groups the cells annotated with same cell type;
3. **boostraps_samples** is the folder in which all the boostrap samples are saved;
4. **genesGeneral_probabilities.tsv** is a tabular file that reports the probability of each gene to be expressed in a generical cell estimated from the boostrap samples;
5. **genesCellTypes_probabilities.tsv** is a table that provides the the probability of each gene to be expressed in any cell type from the annotation. As already mentioned, the probability is estimated from the boostrap samples;
6. **genesProbabilities_ratios.tsv** is a tab separated file reporting the ratios between the two previously mentioned probabilities;
7. **genesRanking.tsv** show the ranking of the genes on the basis of the ratios reported in the genesProbabilities_ratios.tsv file;
8. **genes_entropy.tsv** gives the entropy of each gene calculated over the probabilites of a gene to be expressed in any cell type;
9. **genes2remove.tsv** lists the genes to remove from the final lists;
10. **cellTypes_fromAnnotationHeatmap.png** is an heatmap showing the percentage of overlap among each couple of final cell type specific list of genes;
11. **TABLE_OF_GENES.tsv** is a simple tabular file reporting the genes from the counts in the proper order.

