For impatient people
====================

The tool runs on a Conda envirnoment which can be installed following the steps reported in the **PREREQUISITES** section of this manual.

The input of SCALT is a matrix of **read counts** from single cell RNA sequencing, having the genes on the rows and the cells ids on the columns.
The gene **gene notation** reported in the table must be either **gene_symbol** or **ensembl_id** and the matrix must be in **.tsv** extension. 
The program for classification is called **SCALT.py** 

By default, SCALT classifies each cell to one of the **293** pre-defined cell types available using the followinmg command:

::

   python3 SCALT.py read_counts.tsv -Notation ensembl_id  

The outpus are:

1. a **REPORT.html** file reporting a collection of plots that summarize the results of the experiment;
2. a directory called **results_directory** which contains all the metadata originated from the analysis. 
