Inputs formats
==============

SCALT_NaiveListsBuilder.py demands two inputs:

1. a scRNA seq row counts matrix. The matrix must be in **.tsv** extension;
2. a text file (**.txt**) having **one column** where the user-defined lists must be inserted.

The counts matrix must present genes on the rows and cells on the columns. The first row of the matrix must contain the ids of each cell while the first column must provide the gene ids written either as **gene symbol** or **ensembl id**. 

To see an example of the table, see the section **SCALT: CLASSIFICATION - Inputs & Outputs** of this manual.

An example of the text file is reported here:

.. list-table:: 
   :align: center
   :widths: 80 

   * - -B-cell
   * - LCN8
   * - TRBC2
   * - CCR9
   * - -Langerhans cell
   * - GFI1B
   * - MMP9
   * - ETV3L
   * - -Macrophage cell
   * - RNASE2
   * - CEBPE
   * - SPIC

.. note::
   Avoid the usage of too small lists of genes since they might decrease the accuracy of the tool. Larger lists with small overlaps are recomended.

Parameters
==========

SCALT_NaiveListsBuilder.py takes the advantage of a collection of both positional arguments and parameters that can visualized typing the following command:

:: 

  python3 SCALT_NaiveListsBuilder.py -h

The documentation should appear as follows:

::

   usage: SCALT_NaiveListsBuilder.py [-h] [-Boo --Boostraps] [-Cells --Cells]
                                      [-Genes --Genes] [-Notation --Notation]
                                      [-CPUs --CPUs]
                                      Sample CustomLists


1. **Sample** is the first positional argument and addresses the name of the counts matrix;
2. **CustomLists** is the second positional argument and refers to the **.txt** file having the user-defined custom lists of genes;
3. **-Boo** or **--Boostraps** indicates the number of boostrap samples to generate. The default number is **100**;
4. **-Cells** or **--Cells** specifies the number of cells to pick randomly per each cell type during the probability inference process. By default, the number is set to **100**;
5. **-Genes** or **--Genes** refers to the number of genes that the final cell type lists of genes must contain at the end. The default number is **100**;
6. **-Notation** is used to underline the kind of gene notation present in the counts. The user can choose between **gene_symbol** or **ensembl_id**. By default, ensembl_id is set;
7. **-CPUs** or **--CPUs** indicates the number of processors to use. By default, **1** is used.


Run SCALT_NaiveListsBuilder.py
==================================

SCALT_NaiveListsBuilder.py is quite straightforward sine it requires the counts table and the text file as inputs. 

Leaving default parameters, the basic comand appears as follows:

::

   python3 SCALT_NaiveListsBuilder.py read_counts.tsv mylists.txt

By default, it is assumed that the counts have the ensembl id for genes notation. 100 boostrap samples are generated in which 100 cells will be randomly sampled per each cell type. At the end, the final lists will contain 100 genes each. The number of processors used is 1.

If the **gene symbol** is used in the counts matrix, the notation must be specified as follows:

::

   python3 SCALT_NaiveListsBuilder.py read_counts.tsv mylists.txt -Notation gene_symbol

Or:

::

   python3 SCALT_NaiveListsBuilder.py read_counts.tsv mylists.txt --Notation gene_symbol

The number of boostrap samples can be changed in any moment specifying the number in the proper parameter:

::

   python3 SCALT_NaiveListsBuilder.py read_counts.tsv mylists.txt -Boo 80

Or:

::

   python3 SCALT_NaiveListsBuilder.py read_counts.tsv mylists.txt --Boostraps 80

The number of cell per cell type included in each boostrap sample can be modified making usage of the proper parameter:

::

   python3 SCALT_NaiveListsBuilder.py read_counts.tsv mylists.txt -Cells 50

Or:

::

   python3 SCALT_NaiveListsBuilder.py read_counts.tsv mylists.txt --Cells 50

The same logic for the number of genes that will be included in the final lists of genes:

::

   python3 SCALT_NaiveListsBuilder.py read_counts.tsv mylists.txt -Genes 75

Or:

::

   python3 SCALT_NaiveListsBuilder.py read_counts.tsv mylists.txt --Genes 75

In addition, the computational time can be reduced if the number of processors is increased as reported:

::

   python3 SCALT_NaiveListsBuilder.py read_counts.tsv mylists.txt -CPUs 4

Or:

::

   python3 SCALT_NaiveListsBuilder.py read_counts.tsv mylists.txt --CPUs 4

Make sure to have available the number of desidered processors in your machine.

To conclude, the different parameters can be modified in a unique call:

::

   python3 SCALT_NaiveListsBuilder.py read_counts.tsv mylists.txt -Notation gene_symbol -Boo 80 -Cells 50 -Genes 75 -CPUs 4

Or:

::

   python3 SCALT_NaiveListsBuilder.py read_counts.tsv mylists.txt --Notation gene_symbol --Boostraps 80 --Cells 50 --Genes 75 --CPUs 4

The order of parameters is irrelevant.

Outputs
=======

The tool returns two outputs:

1. a directory called **naive** containing the final lists of genes;
2. a directory named **NaivelistsBuilder_results** hosting a collection of metadata.

The metadata consists in a series of files and directories which are produced automatically during the process and were utilized for the generation of the final lists:
  
1. **originalTables_zipped.zip** is a zipped repository containing the original input data;
2. **FDR_table.tsv** file that is the tabular file containing the **False Discovery Rate** of each hypergeometric done;
3. **naive_annotation.tsv** file is the table reporting the naive annotation file required for the subsequent steps;
4. **groupped_cell_types** is the directory that containg the counts matrix clustered per cell type. Each tsv file groups the cells annotated with same cell type;
5. **boostraps_samples** is the folder in which all the boostrap samples are saved;
6. **genesGeneral_probabilities.tsv** is a tabular file that reports the probability of each gene to be expressed in a generical cell estimated from the boostrap samples;
7. **genesCellTypes_probabilities.tsv** is a table that provides the the probability of each gene to be expressed in any cell type from the naive annotation. As already mentioned, the probability is estimated from the boostrap samples;
8. **genesProbabilities_ratios.tsv** is a tab separated file reporting the ratios between the two previously mentioned probabilities;
9. **genesRanking.tsv** show the ranking of the genes on the basis of the ratios reported in the genesProbabilities_ratios.tsv file;
10. **genes_entropy.tsv** gives the entropy of each gene calculated over the probabilites of a gene to be expressed in any cell type;
11. **genes2remove.tsv** lists the genes to remove from the final lists;
12. **newCellTypes_fromNaiveHeatmap.png** is an heatmap showing the percentage of overlap among each couple of final cell type specific list of genes;
13. **TABLE_OF_GENES.tsv** is a simple tabular file reporting the genes from the counts in the proper order.

