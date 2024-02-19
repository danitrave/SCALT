Formats converter
=================

SCALT posseses an auxiliary utility which is used to extract the scRNA seq counts matrix from files in **.rds**, **.RData** and **.h5ad** and convert it to a **.tsv** format.

The programm is called **formatConverter2tsv.py**. The script requires just one positional argument which is the file to convert. It automatically understands the kind of format.

Its usage is quite straightforward. Example of the command for the three kind of files able to convert are reported below:

::

   python3 formatConverter2tsv.py data.rds

::

   python3 formatConverter2tsv.py data.RData

::

   python3 formatConverter2tsv.py data.h5ad

