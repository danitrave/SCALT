Inputs formats
==============

The required input of SCALT_merge.py is the path to the **custom** directory containing the cell type specific lists of genes deriving from SCALT_AnnotationListBuilder.py.
This input must be specified as positional argument of the tool.

Parameters
==========

SCALT_merge.py takes the advantage of a collection of parameters that can visualized typing the following command:

:: 

  python3 SCALT_merge.py -h

The documentation should appear as follows:

::

   usage: SCALT_merge.py [-h] [-Notation --Notation] [-Genes --Genes] 
                         [-Boo --Bootstrap] [-out --out] Lists


1. **Lists** is the only positional argument and refers to the path where the **custom** directory is located;
2. **-h** or **--help** show the documentation and the parameters;
3. **-Notation** is used to underline the kind of gene notation present in the counts. The user can choose between **gene_symbol** or **ensembl_id**. By default, ensembl_id is set;
4. **-Genes** or **--Genes** refers to the number of genes contained in each cell type specific lists of genes. The default value is 100;
5. **-Boo** or **--Boostraps** indicates the number of boostrap samples available The default number is **100**;
6. **-out** or **--outs** specifies the name of the name of the final report. By default the name is **merging_diagnose.txt**


Run SCALT_merge.py
==================================

SCALT_merge.py is quite straightforward given that only the path to the **custom** directory must be specified. 

Leaving default parameters, the basic comand appears as follows:

::

   python3 SCALT_merge.py custom

By default, it is assumed that 100 bootstrap samples are available and that each cell type specific list of genes in the custom directory has 100 genes. Moreover, ensembl ids are assumed and the final name of the report withb be **mutual_information_diagnose.txt**.

If the **gene symbol** is present in the cell type specific lists of genes, the notation must be specified as follows:

::

   python3 SCALT_merge.py custom -Notation gene_symbol

Or:

::

   python3 SCALT_merge.py custom --Notation gene_symbol

The number of boostrap samples available can be changed in any moment specifying the number in the proper parameter:

::

   python3 SCALT_merge.py custom -Boo 80

Or:

::

   python3 SCALT_merge.py custom --Boostraps 80


The same process can be translated for the number of genes available in the cell type defyning lists:

::

   python3 SCALT_merge.py custom -Genes 75

Or:

::

   python3 SCALT_merge.py custom --Genes 75

The name of the output report can be changed adjusting the proper parameter:

::

   python3 SCALT_merge.py custom -out report_one.txt

Or:

::

   python3 SCALT_merge.py custom --out report_one.txt


To conclude, the different parameters can be modified in a unique call:

::

   python3 SCALT_merge.py custom -Notation gene_symbol -Boo 80 -Genes 75 -out report_one.txt

Or:

::

   python3 SCALT_merge.py custom --Notation gene_symbol --Boostraps 80 --Genes 75 --out report_one.txt

The order of parameters is irrelevant.

Outputs
=======

The tool returns two outputs:

1. a report (default name merging_diagnose.txt) in **.txt** format highlighting which cell types should be merged in a unique one with the corresponding similarity values expressed in terms of mutual information;
2. a directory (default name mi_merge_results) grouping all the additional results and metadata deriving from SCALT_merge.py

The directory contains the following files:
1. **boos_mutual_information.tsv** reports the mutual information values calculated per each pair of bootstrap lists of the same cell type; 
2. **thresholdsMI.tsv** is a tabular file reporting a series of metrics summarizing the distribution of mutual information values deriving from lists of the same cell types, including the cell type specific mutual information threshold; 
3. **miProbTables/** is a directory containing a series of **.tsv** files each reporting the probability of finding each gene across a number of bootstrap lists of the same cell type;
4. **mi/** is a directory containing the bootstrap lists deriving from each bootstrap sample;
5. **experiment_MIs.tsv** is a tabular file reporting the mutual information values calculated between each pair of lists from the **custom** directory, excluding self-comparison.


