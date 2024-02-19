Conda Installation
=====
 
SCALT is written in python, R and bash programming language. To ensure the adaptability accross different operating systems, SCALT requires the installation of a **Conda** environment. Conda is a powerful command line tool for package and environment management that runs on Windows, macOS, and Linux.

To install Conda, visit the `Conda Installation Page <https://conda.io/projects/conda/en/latest/user-guide/install/index.html>`_ and follow the instructions step by step for the regular installation on your operating system.


SCALT download
==============

Subsequent to the installation of conda, make sure that conda is activated with the following command:

::

  conda activate

Then, download SCALT either from the following `Github repository <https://github.com/danitrave/SCALT>`_ as a zipped file or directly using the following command on the command line:

::

  git clone https://github.com/danitrave/SCALT.git


Move to the SCALT directory making use of **cd** command as follows:

::

  cd path/to/directory/SCALT

Inside the directory, the configuration file **SCALT_conda_envSetup.yml** enables the installation of the proper SCALT environment in which all the packages and programs required for the tool are provided. Install the environment running the following command:

::

  conda env create -f SCALT_conda_envSetup.yml

Finally, make sure to activate the environment if not already activated using the following commad:

::

  conda activate SCALT

Now, you are ready to play with SCALT!

