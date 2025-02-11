For impatient people
====================

The tool runs on UNIX systems and requires a Conda envirnoment which can be installed following the steps reported in the **PREREQUISITES** section of this manual.

The input of SCALT is a matrix of **read counts** from single cell RNA sequencing, having the genes on the rows and the cells ids on the columns.
The **gene notation** reported in the table must be either **gene_symbol** or **ensembl_id** and the matrix must be in **.tsv** extension. 
The program for classification is called **SCALT.py** 

By default, SCALT classifies each cell to one of the **293** pre-defined cell types available using the followinmg command:

::

   python3 SCALT.py read_counts.tsv -Notation ensembl_id  

The outpus are:

1. a **report** file in HTML format reporting a collection of plots that summarize the results of the experiment;
2. a table in TSV format highlighting the classification outcome for each singular cell;
3. a directory containing additional results and metadata originated from the analysis. By default, the name of the directory is **results_directory**. 
