For impatient people
====================

The tool must be executed on a Conda envirnoment which can be installed following the steps reported in the **PREREQUISITES** section of this manual.

SCALT requires a scRNA **read counts** matrix, having the genes on the rows and the cells ids on the columns, and the **gene notation** reported in the table, either **gene_symbol** or **ensembl_id**. The program for classification is called **SCALT.py** 

SCALT classifies each cell to one of the **471 cell types** available using the followinmg command:

::

   python3 SCALT.py read_counts.py -Notation ensembl_id  

The outpus are:

1. a **REPORT.html** file reporting a collection of plots that summarize the results of the experiment;
2. a directory called **results_directory** which contains all the metadata originated from the analysis. 
